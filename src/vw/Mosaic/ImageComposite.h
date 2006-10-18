#ifndef __VW_MOSAIC_IMAGECOMPOSITE_H__
#define __VW_MOSAIC_IMAGECOMPOSITE_H__

#include <vector>
#include <list>

#include <vw/Core/Cache.h>
#include <vw/Math/BBox.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/Filter.h>
#include <vw/Image/EdgeExtend.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/Transform.h>

namespace vw {
namespace mosaic {

  // *******************************************************************
  // PositionedImage
  // *******************************************************************

  template <class PixelT>
  class PositionedImage : public ImageViewBase<PositionedImage<PixelT> > {
  public:
    int cols_, rows_;
    ImageView<PixelT> image;
    BBox<int,2> bbox;

    template <class ImageT>
    PositionedImage( int cols, int rows, ImageT const& image, BBox<int,2> const& bbox ) : cols_(cols), rows_(rows), image(image), bbox(bbox) {}

    PositionedImage reduce() const {
      const int border = 1;
      int left = std::min( border + (bbox.min().x()+border)%2, bbox.min().x() );
      int top = std::min( border + (bbox.min().y()+border)%2, bbox.min().y() );
      int right = std::min( border + (bbox.width()+left+border)%2, cols_-bbox.min().x()-bbox.width() );
      int bottom = std::min( border + (bbox.height()+top+border)%2, rows_-bbox.min().y()-bbox.height() );
      std::vector<float> kernel(3); kernel[0]=0.25; kernel[1]=0.5; kernel[2]=0.25;
      // I don't quite yet understand why (if?) this is the correct bounding box,
      // but bad things happen without the final "+1"s:
      BBox<int,2> new_bbox( Vector<int,2>( (bbox.min().x()-left)/2, (bbox.min().y()-top)/2 ),
                            Vector<int,2>( (bbox.min().x()-left)/2 + (bbox.width()+left+right+1)/2+1,
                                                (bbox.min().y()-top)/2 + (bbox.height()+top+bottom+1)/2+1 ) );
      // We use vw::rasterize() here rather than ordinary assignment because it is 
      // faster for this particular combination of filtering and subsampling.
      ImageView<PixelT> new_image( new_bbox.width(), new_bbox.height() );
      rasterize( edge_extend( subsample( separable_convolution_filter( edge_extend( image, -left, -top, image.cols()+left+right, image.rows()+top+bottom, ZeroEdgeExtend() ),
                                                                     kernel, kernel, ZeroEdgeExtend() ), 2 ), 0, 0, new_bbox.width(), new_bbox.height(), vw::ConstantEdgeExtend() ),
                 new_image );
      return PositionedImage( (cols_+1)/2, (rows_+1)/2, new_image, new_bbox );
    }
    
    void unpremultiply() {
      image /= select_channel( image, 3 );
    }

    void addto( ImageView<PixelT> const& dest ) const {
      crop( dest, bbox ) += image;
    }

    // Performs additive composition when overlay==false (the default).
    // When overlay==true, it overlays the image on top of the destination,
    // respecting the alpha channel, but that currently only works for 
    // floating-point RGBA pixels.
    void addto( ImageView<PixelT> const& dest, int ox, int oy, bool overlay = false ) {
      BBox<int,2> sum_bbox = bbox;
      sum_bbox.crop( BBox<int,2>( Vector<int,2>(ox,oy), Vector<int,2>(ox+dest.cols(),oy+dest.rows()) ) );
      if( sum_bbox.empty() ) return;
      if( overlay ) {
        ImageView<float> alpha = select_channel(crop( image, sum_bbox-bbox.min() ),3);
        crop( dest, sum_bbox-Vector<int,2>(ox,oy) ) *= 1.0-alpha;
        crop( dest, sum_bbox-Vector<int,2>(ox,oy) ) += crop( image, sum_bbox-bbox.min() );
      }
      else {
        crop( dest, sum_bbox-Vector<int,2>(ox,oy) ) += crop( image, sum_bbox-bbox.min() );
      }
    }

    void subtract_expanded( PositionedImage const& other ) {
      rasterize( image - edge_extend( resample( other.image, 2 ), bbox-2*other.bbox.min(), ZeroEdgeExtend() ), image );
    }

    template <class OtherPixT>
    PositionedImage& operator*=( PositionedImage<OtherPixT> const& other ) {
      image *= edge_extend( other.image, bbox-other.bbox.min(), ZeroEdgeExtend() );
      return *this;
    }
    
    int cols() const { return cols_; }
    int rows() const { return rows_; }
    int planes() const { return 1; }
  };


  // *******************************************************************
  // ImageComposite
  // *******************************************************************

  class ImageComposite : public ImageViewBase<ImageComposite> {

    struct Pyramid {
      std::vector<PositionedImage<PixelRGBA<float> > > images;
      std::vector<PositionedImage<float> > masks;
    };

    class SourceGenerator {
      ImageViewRef<PixelRGBA<float> > m_source;
    public:
      typedef ImageView<PixelRGBA<float> > value_type;
      SourceGenerator( ImageViewRef<PixelRGBA<float> > const& source ) : m_source(source) {}
      size_t size() const {
        return m_source.cols() * m_source.rows() * sizeof(PixelRGBA<float>);
      }
      boost::shared_ptr<value_type> generate() const {
        return boost::shared_ptr<value_type>( new value_type(m_source) );
      }
    };
    
    class AlphaGenerator {
    public:
      ImageComposite& m_composite;
      int m_index;
    public:
      typedef ImageView<float> value_type;
      AlphaGenerator( ImageComposite& composite, int index ) : m_composite(composite), m_index(index) {}
      size_t size() const {
        return m_composite.sources[m_index].size() / 4;
      }
      boost::shared_ptr<value_type> generate() const { 
        ImageView<PixelRGBA<float> > source = *m_composite.sources[m_index];
        m_composite.sources[m_index].deprioritize();
        return boost::shared_ptr<value_type>( new value_type( select_channel( source, 3 ) ) );
      }
    };

    class PyramidGenerator {
    public:
      ImageComposite& m_composite;
      int m_index;
    public:
      typedef Pyramid value_type;
      PyramidGenerator( ImageComposite& composite, int index ) : m_composite(composite), m_index(index) {}
      size_t size() const {
        return size_t( m_composite.sources[m_index].size() * 1.66 ); // 1.66 = (5/4)*(4/3)
      }
      boost::shared_ptr<value_type> generate() const;
    };

    friend class PyramidGenerator;

    std::vector<BBox<int,2> > bboxes;
    BBox<int,2> bbox;
    int mindim, levels;
    bool m_draft_mode;
    Cache& m_cache;
    std::vector<Cache::Handle<SourceGenerator> > sources;
    std::vector<Cache::Handle<AlphaGenerator> > alphas;
    std::vector<Cache::Handle<PyramidGenerator> > pyramids;

    void generate_masks() const;

    ImageView<PixelRGBA<float> > blend_patch( BBox<int,2> const& patch_bbox ) const;
    ImageView<PixelRGBA<float> > draft_patch( BBox<int,2> const& patch_bbox ) const;

  public:
    typedef PixelRGBA<float> pixel_type;
    typedef PixelRGBA<float> result_type;
    
    ImageComposite() : m_draft_mode(false), m_cache(Cache::system_cache()) {}

    void insert( ImageViewRef<PixelRGBA<float> > const& image, int x, int y );

    void prepare();

    ImageView<PixelRGBA<float> > generate_patch( BBox<int,2> const& patch_bbox ) const {
      if( m_draft_mode ) return draft_patch( patch_bbox );
      else return blend_patch( patch_bbox );
    }

    void set_draft_mode( bool draft_mode ) { m_draft_mode = draft_mode; }

    int cols() const {
      return bbox.width();
    }

    int rows() const {
      return bbox.height();
    }

    int planes() const {
      return 1;
    }

    PixelRGBA<float> operator()( int x, int y, int p=0 ) const {
      throw NoImplErr() << "ImageComposite does not support individual pixel access!";
    }

    typedef ProceduralPixelAccessor<ImageComposite> pixel_accessor;
    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }    

    typedef ImageComposite prerasterize_type;
    inline prerasterize_type const& prerasterize( BBox2i const& bbox ) const { return *this; }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const { dest = generate_patch(bbox); }

  };

} // namespace mosaic
} // namespace vw

#endif // __VW_MOSAIC_IMAGECOMPOSITE_H__
