// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


/// \file ImageComposte.h
///
/// A view class that represents a composite mosaic of images.
///
#ifndef __VW_MOSAIC_IMAGECOMPOSITE_H__
#define __VW_MOSAIC_IMAGECOMPOSITE_H__

#include <iostream>
#include <vector>
#include <list>

#include <vw/Core/Cache.h>
#include <vw/Core/ProgressCallback.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/Grassfire.h>
#include <vw/Image/Transform.h>
#include <vw/Image/Filter.h>
#include <vw/Image/SparseImageCheck.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/ImageChannels.h>

namespace vw {
namespace mosaic {

  // *******************************************************************
  // PositionedImage
  // *******************************************************************

  /// ?
  template <class PixelT>
  class PositionedImage : public ImageViewBase<PositionedImage<PixelT> > {
  public:
    int m_cols, m_rows;
    ImageView<PixelT> image;
    BBox2i bbox;

    typedef PixelT pixel_type;
    typedef PixelT result_type;

    template <class ImageT>
    PositionedImage( int cols, int rows, ImageT const& image, BBox2i const& bbox ) 
      : m_cols(cols), m_rows(rows), image(image), bbox(bbox) {}

    /// ?
    PositionedImage reduce() const {
      const int border = 1;

      int left   = std::min( border + (bbox.min().x()   +border)%2, bbox.min().x()                      );
      int top    = std::min( border + (bbox.min().y()   +border)%2, bbox.min().y()                      );
      int right  = std::min( border + (bbox.width()+left+border)%2, m_cols-bbox.min().x()-bbox.width()  );
      int bottom = std::min( border + (bbox.height()+top+border)%2, m_rows-bbox.min().y()-bbox.height() );

      std::vector<float> kernel(3); kernel[0]=0.25; kernel[1]=0.5; kernel[2]=0.25;
      // I don't quite yet understand why (if?) this is the correct bounding box,
      // but bad things happen without the final "+1"s:
      BBox2i new_bbox( Vector2i( (bbox.min().x()-left)/2, 
                                 (bbox.min().y()-top )/2 ),
                       Vector2i( (bbox.min().x()-left)/2 + (bbox.width()+left+right +1)/2+1,
                                 (bbox.min().y()-top )/2 + (bbox.height()+top+bottom+1)/2+1 ) );
      // We use vw::rasterize() here rather than ordinary assignment because it is
      // faster for this particular combination of filtering and subsampling.
      ImageView<PixelT> new_image( new_bbox.width(), new_bbox.height() );
      vw::rasterize( edge_extend( subsample( separable_convolution_filter( edge_extend( image, -left, -top, image.cols()+left+right, image.rows()+top+bottom, ZeroEdgeExtension() ),
                                                                           kernel, kernel, ZeroEdgeExtension() ), 2 ), 0, 0, new_bbox.width(), new_bbox.height(), vw::ConstantEdgeExtension() ),
                     new_image );
      return PositionedImage( (m_cols+1)/2, (m_rows+1)/2, new_image, new_bbox );
    }

    void unpremultiply() {
      image /= select_alpha_channel( image );
    }

    void addto( ImageView<PixelT> const& dest ) const {
      crop( dest, bbox ) += image;
    }

    // Performs additive composition when overlay==false (the default).
    // When overlay==true, it overlays the image on top of the destination,
    // respecting any alpha channel.
    void addto( ImageView<PixelT> const& dest, int ox, int oy, bool overlay = false ) {
      BBox2i sum_bbox = bbox;
      sum_bbox.crop( BBox2i( Vector2i(ox,oy), Vector2i(ox+dest.cols(),oy+dest.rows()) ) );
      if( sum_bbox.empty() ) return;
      if( overlay ) {
        if( PixelHasAlpha<PixelT>::value ) {
          crop( dest, sum_bbox-Vector2i(ox,oy) ) *= 1.0 - select_alpha_channel( crop( image, sum_bbox-bbox.min() ) ) / (double)ChannelRange<PixelT>::max();
          crop( dest, sum_bbox-Vector2i(ox,oy) ) += crop( image, sum_bbox-bbox.min() );
        }
        else {
          crop( dest, sum_bbox-Vector2i(ox,oy) ) = crop( image, sum_bbox-bbox.min() );
        }
      }
      else {
        crop( dest, sum_bbox-Vector2i(ox,oy) ) += crop( image, sum_bbox-bbox.min() );
      }
    }

    void subtract_expanded( PositionedImage const& other ) {
      vw::rasterize( image - edge_extend( resample( other.image, 2 ), bbox-2*other.bbox.min(), ZeroEdgeExtension() ), image );
    }

    template <class OtherPixT>
    PositionedImage& operator*=( PositionedImage<OtherPixT> const& other ) {
      image *= edge_extend( other.image, bbox-other.bbox.min(), ZeroEdgeExtension() );
      return *this;
    }

    int32 cols  () const { return m_cols; }
    int32 rows  () const { return m_rows; }
    int32 planes() const { return 1;      }

    typedef PositionedImage prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& /*bbox*/ ) const { return *this; }
    template <class DestT> inline void rasterize( DestT const& /*dest*/, BBox2i /*bbox*/ ) const {
      vw_throw( NoImplErr() << "PositionedImage does not support rasterize!" );
    }
  };


  // *******************************************************************
  // ImageComposite
  // *******************************************************************

  /// ?
  template <class PixelT>
  class ImageComposite : public ImageViewBase<ImageComposite<PixelT> > {
  public:
    typedef PixelT pixel_type;
    typedef typename PixelChannelType<PixelT>::type channel_type;

  private:
    struct Pyramid {
      std::vector<PositionedImage<pixel_type  > > images;
      std::vector<PositionedImage<channel_type> > masks;
    };

    /// ?
    class SourceGenerator {
      ImageViewRef<pixel_type> m_source;
    public:
      typedef ImageView<pixel_type> value_type;
      SourceGenerator( ImageViewRef<pixel_type> const& source ) : m_source(source) {}
      size_t size() const {
        return m_source.cols() * m_source.rows() * sizeof(pixel_type);
      }
      boost::shared_ptr<value_type> generate() const {
        return boost::shared_ptr<value_type>( new value_type(m_source) );
      }
    };

    class GrassfireGenerator {
      ImageViewRef<pixel_type> m_source;
    public:
      typedef ImageView<float32> value_type;
      GrassfireGenerator( ImageViewRef<pixel_type> const& source ) : m_source(source) {}
      size_t size() const {
        return m_source.cols() * m_source.rows() * sizeof(float32);
      }
      boost::shared_ptr<value_type> generate() const {
        ImageView<float32> alpha = channel_cast<float32>( select_alpha_channel( copy( m_source ) ) );
        return boost::shared_ptr<value_type>( new value_type( grassfire( alpha ) ) );
      }
    };

    class AlphaGenerator {
    public:
      ImageComposite& m_composite;
      size_t m_index;
    public:
      typedef ImageView<channel_type> value_type;
      AlphaGenerator( ImageComposite& composite, size_t index ) : m_composite(composite), m_index(index) {}
      size_t size() const {
        return m_composite.sources[m_index].size() / PixelNumChannels<pixel_type>::value;
      }
      boost::shared_ptr<value_type> generate() const {
        ImageView<pixel_type> source = *m_composite.sources[m_index];
        m_composite.sources[m_index].release();
        m_composite.sources[m_index].deprioritize();
        return boost::shared_ptr<value_type>( new value_type( select_alpha_channel( source ) ) );
      }
    };

    class PyramidGenerator {
    public:
      ImageComposite& m_composite;
      size_t m_index;
    public:
      typedef Pyramid value_type;
      PyramidGenerator( ImageComposite& composite, size_t index ) : m_composite(composite), m_index(index) {}
      size_t size() const {
        return size_t( double(m_composite.sources[m_index].size()) * 1.66 ); // 1.66 = (5/4)*(4/3)
      }
      boost::shared_ptr<value_type> generate() const;
    };

    friend class PyramidGenerator;

    std::vector<BBox2i > bboxes;
    BBox2i view_bbox, data_bbox;
    int    mindim, levels;
    bool   m_draft_mode;
    bool   m_fill_holes;
    bool   m_reuse_masks;
    Cache& m_cache;
    std::vector<ImageViewRef<pixel_type> >        sourcerefs;
    std::vector<Cache::Handle<SourceGenerator > > sources;
    std::vector<Cache::Handle<AlphaGenerator  > > alphas;
    std::vector<Cache::Handle<PyramidGenerator> > pyramids;

    void generate_masks( ProgressCallback const& progress_callback ) const;

    /// Generates a full-resolution patch of the mosaic corresponding
    /// to the given bounding box.
    ImageView<pixel_type> blend_patch( BBox2i const& patch_bbox ) const;

    // Generates a full-resolution patch of the mosaic corresponding
    // to the given bounding box WITHOUT blending.
    ImageView<pixel_type> draft_patch( BBox2i const& patch_bbox ) const;

  public:
    typedef pixel_type result_type;

    ImageComposite() : m_draft_mode (false), m_fill_holes(false),
                       m_reuse_masks(false), m_cache(vw_system_cache()) {}

    void insert( ImageViewRef<pixel_type> const& image, int x, int y );

    void prepare( const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() );
    void prepare( BBox2i const& total_bbox, const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() );

    /// Generate a section of the output image.
    ImageView<pixel_type> generate_patch( BBox2i const& patch_bbox ) const {
      if( m_draft_mode ) return draft_patch( patch_bbox );
      else return blend_patch( patch_bbox );
    }

    /// If draft mode is on no image blending is performed.
    void set_draft_mode (bool draft_mode ) { m_draft_mode = draft_mode; }

    void set_fill_holes (bool fill_holes ) { m_fill_holes = fill_holes; }

    void set_reuse_masks(bool reuse_masks) { m_reuse_masks = reuse_masks; }

    int32 cols  () const { return view_bbox.width();  }
    int32 rows  () const { return view_bbox.height(); }
    int32 planes() const { return 1;                  }
        
    BBox2i const& bbox            () const { return data_bbox; }
    BBox2i const& source_data_bbox() const { return view_bbox; }

    pixel_type operator()( int x, int y, int p=0 ) const {
      // FIXME: This is horribly slow, and totally untested for
      // multi-band blending.  We should do something faster for draft
      // mode, and possibly cache output blocks in multi-band mode?
      return generate_patch(BBox2i(x,y,1,1))(0,0,p);
    }

    typedef ProceduralPixelAccessor<ImageComposite> pixel_accessor;
    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

    typedef CropView<ImageView<PixelT> > prerasterize_type;

    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      ImageView<PixelT> buf = generate_patch(bbox);
      return CropView<ImageView<PixelT> >( buf, BBox2i(-bbox.min().x(),-bbox.min().y(),cols(),rows()) );
    }

    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
      vw::rasterize( prerasterize(bbox), dest, bbox );
    }

    bool sparse_check( BBox2i const& bbox ) const {
      for (unsigned int i = 0; i < bboxes.size(); ++i) {
        BBox2i src_bbox = bboxes[i];
        src_bbox.crop(bbox);
        if( ! src_bbox.empty() ) {
          if( vw::sparse_check( sourcerefs[i], src_bbox-bboxes[i].min() ) ) {
            return true;
          }
        }
      }
      return false;
    }

  };

} // namespace mosaic

  // This specializes the SparseImageCheck template for the
  // ImageComposite type of ImageView.  The ImageComposite might be
  // very sparsely covered by images, so it is sometimes handy to be
  // able to check to see whether an arbitrary bbox intersects with
  // any of the the ImageComposite's source images.
  template <class PixelT>
  class SparseImageCheck<mosaic::ImageComposite<PixelT> > {
    mosaic::ImageComposite<PixelT> const& composite;
  public:
    SparseImageCheck(mosaic::ImageComposite<PixelT> const& source) : composite(source) {}
    bool operator() (BBox2i const& bbox) {
      return composite.sparse_check(bbox);
    }
  }; // End class ImageComposite

} // namespace vw


template <class PixelT>
void vw::mosaic::ImageComposite<PixelT>::generate_masks( vw::ProgressCallback const& progress_callback ) const {
  vw_out(DebugMessage, "mosaic") << "Generating masks..." << std::endl;
  std::vector<Cache::Handle<GrassfireGenerator> > grassfires;
  for( unsigned i=0; i<sources.size(); ++i )
    grassfires.push_back( m_cache.insert( GrassfireGenerator( sourcerefs[i] ) ) );
  for( unsigned p1=0; p1<sources.size(); ++p1 ) {
    ImageView<float> mask = copy( *(grassfires[p1]) );
    for( unsigned p2=0; p2<sources.size(); ++p2 ) {
      if( p1 == p2 ) continue;
      int ox = bboxes[p2].min().x() - bboxes[p1].min().x();
      int oy = bboxes[p2].min().y() - bboxes[p1].min().y();
      if( ! ( ox >= bboxes[p1].width() ||
              oy >= bboxes[p1].height() ||
              -ox >= bboxes[p2].width() ||
              -oy >= bboxes[p2].height() ) )
      {
        ImageView<float> other = *grassfires[p2];
        int left = std::max( ox, 0 );
        int top = std::max( oy, 0 );
        int right = std::min( bboxes[p2].width()+ox, bboxes[p1].width() );
        int bottom = std::min( bboxes[p2].height()+oy, bboxes[p1].height() );
        for( int j=top; j<bottom; ++j ) {
          for( int i=left; i<right; ++i ) {
            if( ( other(i-ox,j-oy) > mask(i,j) ) ||
                ( other(i-ox,j-oy) == mask(i,j) && p2 > p1 ) )
              mask(i,j) = 0;
          }
        }
      }
      progress_callback.report_fractional_progress( double(p1*(sources.size()+1)+p2+1), double((sources.size()+1)*sources.size()) );
    }
    mask = threshold( mask );
    std::ostringstream filename;
    filename << "mask." << p1 << ".png";
    write_image( filename.str(), mask );
    progress_callback.report_fractional_progress( double((p1+1)*(sources.size()+1)), double((sources.size()+1)*sources.size()) );
  }
  // report_finished() called by prepare(), so don't call it here
}


template <class PixelT>
boost::shared_ptr<typename vw::mosaic::ImageComposite<PixelT>::Pyramid> vw::mosaic::ImageComposite<PixelT>::PyramidGenerator::generate() const {
  vw_out(DebugMessage, "mosaic") << "ImageComposite generating pyramid " << m_index << std::endl;
  boost::shared_ptr<Pyramid> ptr( new Pyramid );
  ImageView<pixel_type> source = copy(*m_composite.sources[m_index]);
  m_composite.sources[m_index].release();
  m_composite.sources[m_index].deprioritize();

  // This is sort of a kluge: the hole-filling algorithm currently
  // doesn't cope well with partially-transparent source pixels.
  if( m_composite.m_fill_holes ) source /= select_alpha_channel(source);

  PositionedImage<pixel_type> image_high( m_composite.view_bbox.width(), m_composite.view_bbox.height(), source, m_composite.bboxes[m_index] );
  PositionedImage<pixel_type> image_low = image_high.reduce();
  ImageView<channel_type> mask_image;

  std::ostringstream mask_filename;
  mask_filename << "mask." << m_index << ".png";
  read_image( mask_image, mask_filename.str() );
  PositionedImage<channel_type> mask( m_composite.view_bbox.width(), m_composite.view_bbox.height(), mask_image, m_composite.bboxes[m_index] );

  for( int l=0; l<m_composite.levels; ++l ) {
    PositionedImage<pixel_type> diff = image_high;
    if( l > 0 ) mask = mask.reduce();
    if( l < m_composite.levels-1 ) {
      PositionedImage<pixel_type> next_image_low = image_low.reduce();
      image_low.unpremultiply();
      diff.subtract_expanded( image_low );
      image_high = image_low;
      image_low = next_image_low;
    }
    diff *= mask;
    ptr->images.push_back( diff );
    ptr->masks.push_back( mask );
  }
  return ptr;
}


template <class PixelT>
void vw::mosaic::ImageComposite<PixelT>::insert( ImageViewRef<pixel_type> const& image, int x, int y ) {
  sourcerefs.push_back( image );
  sources.push_back( m_cache.insert( SourceGenerator( image ) ) );
  alphas.push_back( m_cache.insert( AlphaGenerator( *this, pyramids.size() ) ) );
  pyramids.push_back( m_cache.insert( PyramidGenerator( *this, pyramids.size() ) ) );

  int cols = image.cols(), rows = image.rows();
  BBox2i image_bbox( Vector2i(x, y), Vector2i(x+cols, y+rows) );
  bboxes.push_back( image_bbox );
  if( bboxes.size() == 1 ) {
    view_bbox = bboxes.back();
    data_bbox = bboxes.back();
    mindim = std::min(cols,rows);
  }
  else {
    view_bbox.grow( image_bbox );
    data_bbox.grow( image_bbox );
    mindim = std::min( mindim, std::min(cols,rows) );
  }
}


template <class PixelT>
void vw::mosaic::ImageComposite<PixelT>::prepare( vw::ProgressCallback const& progress_callback ) {
  // Translate bboxes to origin
  for( unsigned i=0; i<sources.size(); ++i )
    bboxes[i] -= view_bbox.min();
  data_bbox -= view_bbox.min();

  levels = (int) floorf( logf( float(mindim)/2.0f ) / logf(2.0f) ) - 1;
  if( levels < 1 ) levels = 1;

  if( !m_draft_mode && !m_reuse_masks ) {
    generate_masks( progress_callback );
  }
  progress_callback.report_finished();
}

template <class PixelT>
void vw::mosaic::ImageComposite<PixelT>::prepare( BBox2i const& total_bbox,
                                                  vw::ProgressCallback const& progress_callback ) {
  view_bbox = total_bbox;
  prepare( progress_callback );
}

// Suppose a destination image patch at a given level of the pyramid
// has a bounding box that begins at offset x and has width w.  It
// is affected by a range of pixels at the next level of the pyramid
// starting at x/2 with width (x+w)/2-x/2+1 = (w+x%2)/2+1.  This in
// turn is affected by source image pixels at the current level in
// the range starting at 2*(x/2)-1 = x-x%2-1 with width
// (2*(x+w)/2+1)-(2*(x/2)-1)+1 = w-(x+w)%2+x%2+3.

// Generates a full-resolution patch of the mosaic corresponding
// to the given bounding box.
template <class PixelT>
vw::ImageView<PixelT> vw::mosaic::ImageComposite<PixelT>::blend_patch( BBox2i const& patch_bbox ) const {
#if VW_DEBUG_LEVEL > 1
  vw_out(DebugMessage, "mosaic") << "ImageComposite compositing patch " << patch_bbox << "..." << std::endl;
#endif
  // Compute bboxes and allocate the pyramids
  std::vector<BBox2i> bbox_pyr;
  std::vector<ImageView<pixel_type> > sum_pyr(levels);
  std::vector<ImageView<channel_type> > msum_pyr(levels);
  for( int l=0; l<levels; ++l ) {
    if( l==0 ) bbox_pyr.push_back( patch_bbox );
    else bbox_pyr.push_back( BBox2i( Vector2i( bbox_pyr[l-1].min().x() / 2,
                                               bbox_pyr[l-1].min().y() / 2 ),
                                     Vector2i( bbox_pyr[l-1].min().x() / 2 + ( bbox_pyr[l-1].width() + bbox_pyr[l-1].min().x() % 2 ) / 2 + 1,
                                               bbox_pyr[l-1].min().y() / 2 + ( bbox_pyr[l-1].height() + bbox_pyr[l-1].min().y() % 2 ) / 2 + 1) ) );
    sum_pyr [l] = ImageView<pixel_type  >( bbox_pyr[l].width(), bbox_pyr[l].height() );
    msum_pyr[l] = ImageView<channel_type>( bbox_pyr[l].width(), bbox_pyr[l].height() );
  }

  // Compute the bounding box for source pixels that could
  // impact the patch.
  BBox2i padded_bbox = patch_bbox;
  for( int l=0; l<levels-1; ++l ) {
    padded_bbox.min().x() = padded_bbox.min().x()/2;
    padded_bbox.min().y() = padded_bbox.min().y()/2;
    padded_bbox.max().x() = padded_bbox.max().x()/2+1;
    padded_bbox.max().y() = padded_bbox.max().y()/2+1;
  }
  for( int l=0; l<levels-1; ++l ) {
    padded_bbox.min().x() = 2*padded_bbox.min().x()-1;
    padded_bbox.min().y() = 2*padded_bbox.min().y()-1;
    padded_bbox.max().x() = 2*padded_bbox.max().x();
    padded_bbox.max().y() = 2*padded_bbox.max().y();
  }

  // Make a list of the images whose bounding boxes permit them to
  // impact the patch, prioritizing ones that are already in memory.
  std::list<unsigned> image_list;
  for( unsigned p=0; p<sources.size(); ++p ) {
    if( ! padded_bbox.intersects( bboxes[p] ) ) continue;
    if( ! pyramids[p].valid() ) image_list.push_back( p );
    else image_list.push_front( p );
  }

  // Add each source image pyramid to the blend pyramid.
  std::list<unsigned>::iterator ili=image_list.begin(), ilend=image_list.end();
  for( ; ili!=ilend; ++ili ) {
    unsigned p = *ili;
    boost::shared_ptr<Pyramid> pyr = pyramids[p];
    for( int l=0; l<levels; ++l ) {
      pyr->images[l].addto( sum_pyr [l], bbox_pyr[l].min().x(), bbox_pyr[l].min().y() );
      pyr->masks [l].addto( msum_pyr[l], bbox_pyr[l].min().x(), bbox_pyr[l].min().y() );
    }
    pyramids[p].release();
  }

  // Collapse the pyramid
  ImageView<pixel_type> composite( sum_pyr[levels-1].cols(), sum_pyr[levels-1].rows() );
  for( int l=levels; l; --l ) {
    if( l < levels ) {
      composite = ImageView<pixel_type>( crop( resample( composite, 2 ),
                                               bbox_pyr[l-1].min().x()-2*bbox_pyr[l].min().x(),
                                               bbox_pyr[l-1].min().y()-2*bbox_pyr[l].min().y(),
                                               sum_pyr[l-1].cols(), sum_pyr[l-1].rows() ) );
    }
    composite += sum_pyr[l-1] / msum_pyr[l-1];
    sum_pyr.pop_back();
    msum_pyr.pop_back();
  }

  if( m_fill_holes ) {
    composite /= select_alpha_channel( composite );
  }
  else {

    // Trim to the maximal source alpha, reloading images if needed
    ImageView<channel_type> alpha( patch_bbox.width(), patch_bbox.height() );
    for( unsigned p=0; p<sources.size(); ++p ) {
      if( ! patch_bbox.intersects( bboxes[p] ) ) continue;

      ImageView<channel_type> source_alpha = *alphas[p];
      alphas[p].release();

      BBox2i overlap = patch_bbox;
      overlap.crop( bboxes[p] );
      for( int j=0; j<overlap.height(); ++j ) {
        for( int i=0; i<overlap.width(); ++i ) {
          if( source_alpha( overlap.min().x()+i-bboxes[p].min().x(), overlap.min().y()+j-bboxes[p].min().y() ) >
              alpha( overlap.min().x()+i-patch_bbox.min().x(), overlap.min().y()+j-patch_bbox.min().y() ) ) {
            alpha( overlap.min().x()+i-patch_bbox.min().x(), overlap.min().y()+j-patch_bbox.min().y() ) =
              source_alpha( overlap.min().x()+i-bboxes[p].min().x(), overlap.min().y()+j-bboxes[p].min().y() );
          }
        }
      }
    }

    composite *= alpha / select_alpha_channel( composite );
  }

  return composite;
}


// Generates a full-resolution patch of the mosaic corresponding
// to the given bounding box WITHOUT blending.
template <class PixelT>
vw::ImageView<PixelT> vw::mosaic::ImageComposite<PixelT>::draft_patch( BBox2i const& patch_bbox ) const {
#if VW_DEBUG_LEVEL > 1
  vw_out(DebugMessage, "mosaic") << "ImageComposite compositing patch " << patch_bbox << "..." << std::endl;
#endif
  ImageView<pixel_type> composite(patch_bbox.width(),patch_bbox.height());

  // Add each image to the composite.
  for( unsigned p=0; p<sources.size(); ++p ) {
    if( ! patch_bbox.intersects( bboxes[p] ) ) continue;
    BBox2i bbox = patch_bbox;
    bbox.crop( bboxes[p] );
    PositionedImage<pixel_type> image( view_bbox.width(), view_bbox.height(),
                                       crop(sourcerefs[p],bbox-bboxes[p].min()), bbox );
    image.addto( composite, patch_bbox.min().x(), patch_bbox.min().y(), true );
  }

  return composite;
}

#endif // __VW_MOSAIC_IMAGECOMPOSITE_H__
