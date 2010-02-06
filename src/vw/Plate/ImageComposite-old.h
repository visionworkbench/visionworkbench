// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_IMAGE_COMPOSITE_H__
#define __VW_PLATE_IMAGE_COMPOSITE_H__

// TOOD: Once this code has finished maturing, it should probably be
// moved into the mosaic module.

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/EdgeExtension.h>

namespace vw {
namespace platefile {

  template <class PixelT>
  class CompositeView : public ImageViewBase<CompositeView<PixelT> > {

    std::vector<ImageViewRef<PixelT> > m_images;
    std::vector<BBox2i> m_bboxes;
    BBox2i m_view_bbox, m_data_bbox;

    ImageView<PixelT> generate_patch( BBox2i const& patch_bbox ) const;

  public:
    typedef PixelT pixel_type;
    typedef PixelT result_type;
    typedef ProceduralPixelAccessor<CompositeView> pixel_accessor;

    CompositeView() {}
    CompositeView(std::vector<ImageViewRef<PixelT> > images,
                  std::vector<BBox2i> bboxes) : m_images(images),
                                                m_bboxes(bboxes) {}

    void insert( ImageViewRef<PixelT> const& image, int x, int y );

    void prepare( const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() );
    void prepare( BBox2i const& total_bbox, 
                  const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() );

    /// Returns the number of images that intersect with the supplied
    /// bbox.  Useful if you are deciding whether it is worthwhile to
    /// rasterize this bounding box.
    int intersects(BBox2i bbox) const {
      int num_intersections = 0;
      for (unsigned i = 0; i < m_bboxes.size() ; ++i) {
        if(bbox.intersects( m_bboxes[i] ) ) 
          ++num_intersections;
      }
      return num_intersections;
    }

    inline int32 cols() const { return m_data_bbox.width(); }
    inline int32 rows() const { return m_data_bbox.height(); }
    inline int32 planes() const { return 1; }

    inline pixel_accessor origin() const { return pixel_accessor(*this); }
  
    inline result_type operator()( int32 i, int32 j, int32 /*p=0*/ ) const {
      for (int p = m_images.size()-1; p >= 0; --p) {
        if( m_bboxes[p].contains( Vector2i(i,j) ) ) { 
          int ii = i - m_bboxes[p].min().x();
          int jj = j - m_bboxes[p].min().y();
          result_type px = m_images[p](ii,jj);
          if ( ! is_transparent(px) ) 
            return px;
        }
      }
      return result_type();
    }

    /// \cond INTERNAL
    typedef CompositeView<PixelT> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      //      std::cout << "CompositeView()::prerasterize called with bbox=" << bbox << "\n";
      std::vector<ImageViewRef<PixelT> > cropped_images(m_images.size());
      std::vector<BBox2i> cropped_bboxes(m_images.size());
      for (unsigned int i = 0; i < m_images.size(); ++i) {
        if( ! bbox.intersects( m_bboxes[i] ) ) { 

          // We don't bother to rasterize images that don't overlap with
          // the current bounding box.
          cropped_images[i] = m_images[i];
          cropped_bboxes[i] = m_bboxes[i];
         
        } else {

          // For the rest of the images, we create cropped, edge
          // extended versions of the image and modify the bounding
          // box.
          ImageView<PixelT> patch = crop( edge_extend(m_images[i], ZeroEdgeExtension()), 
                                          bbox-m_bboxes[i].min());
          cropped_images[i] = patch;
          cropped_bboxes[i] = bbox;
        }
      }
      return prerasterize_type(cropped_images, cropped_bboxes);
    }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const {
      //      std::cout << "CompositeView()::rasterize called with bbox=" << bbox << "\n";
      vw::rasterize( prerasterize(bbox), dest, bbox );
    }
    /// \endcond
  };

}} // wv::plate


template <class PixelT>
void vw::platefile::CompositeView<PixelT>::insert( ImageViewRef<PixelT> const& image, int x, int y ) {
  m_images.push_back( image );
  int cols = image.cols(), rows = image.rows();
  BBox2i image_bbox( Vector2i(x, y), Vector2i(x+cols, y+rows) );
  m_bboxes.push_back( image_bbox );
  if( m_bboxes.size() == 1 ) {
    m_view_bbox = m_bboxes.back();
    m_data_bbox = m_bboxes.back();
  }
  else {
    m_view_bbox.grow( image_bbox );
    m_data_bbox.grow( image_bbox );
  }
}


template <class PixelT>
void vw::platefile::CompositeView<PixelT>::prepare( vw::ProgressCallback const& progress_callback ) {
  // Translate bboxes to origin
  for( unsigned i=0; i<m_images.size(); ++i )
    m_bboxes[i] -= m_view_bbox.min();
  m_data_bbox -= m_view_bbox.min();
}

template <class PixelT>
void vw::platefile::CompositeView<PixelT>::prepare( BBox2i const& total_bbox,
                                                    vw::ProgressCallback const& progress_callback ) {
  m_view_bbox = total_bbox;
  prepare( progress_callback );
}

#endif
