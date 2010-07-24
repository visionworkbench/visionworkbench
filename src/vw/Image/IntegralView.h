// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

// IntegralView.h
//
// Defines useful images views that create and use Integral Image Views.

#ifndef __VW_IMAGE_INTEGRAL_VIEW_H__
#define __VW_IMAGE_INTEGRAL_VIEW_H__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Image/Manipulation.h>
#include <boost/foreach.hpp>

namespace {
  template <typename T>
  int cmp_range(const T& value, const T& min_, const T& max_) {
    if (value < min_)
      return -1;
    else if (value > max_)
      return 1;
    return 0;
  }
}

namespace vw {
  // An edge extention type that obeys the integral image expression:
  //    sum(x,y) = image(x,y) + sum(x-1,y) + sum(x,y-1) - sum(x-1,y-1)
  struct IntegralEdgeExtension : EdgeExtensionBase {
      template <class ViewT>
      inline typename ViewT::pixel_type operator()( const ViewT &view, int32 i, int32 j, int32 p ) const {

        const int xr = cmp_range(i, 0, view.cols()-1);
        const int yr = cmp_range(j, 0, view.rows()-1);

        if(xr == 0 && yr == 0)
          return view(i,j,p);
        else if (xr == -1 || yr == -1)
          return typename ViewT::pixel_type();
        else if (xr == 1 && yr == 0 )
          return view(view.cols()-1,j,p);
        else if (xr == 0 && yr == 1 )
          return view(i,view.rows()-1,p);
        else
          return view(view.cols()-1,view.rows()-1,p);
      }

      template <class ViewT>
      inline BBox2i source_bbox( ViewT const& view, BBox2i const& bbox ) const {
        BBox2i result = bbox;
        result.crop( BBox2i( 0, 0, view.cols(), view.rows() ) );
        return result;
      }
  };

  template <class ImageT>
  inline ImageView<typename CompoundChannelCast<typename ImageT::pixel_type,typename AccumulatorType<typename CompoundChannelType<typename ImageT::pixel_type>::type>::type>::type>
  integral_image( ImageViewBase<ImageT> const& src ) {

    typedef typename ImageT::pixel_type input_pixel_type;
    typedef typename CompoundChannelCast<input_pixel_type,typename AccumulatorType<typename CompoundChannelType<input_pixel_type>::type>::type>::type pixel_type;
    typedef ImageT SrcT;
    typedef ImageView<pixel_type> DstT;

    ImageView<pixel_type> dst( src.impl().cols(), src.impl().rows() );

    typedef typename ImageView<pixel_type>::pixel_accessor DstAccessT;
    typedef typename ImageT::pixel_accessor SrcAccessT;

    SrcAccessT splane  = src.impl().origin();
    DstAccessT dplane = dst.impl().origin();

    // seed the first row and column
    for( int32 p=0; p<dst.planes(); ++p ) {
      *dplane = *splane;
      SrcAccessT srow = splane.next_row_copy();
      DstAccessT drow = dplane.next_row_copy();
      for( int32 y=1; y<dst.rows(); ++y ) {
        *drow = *srow + *drow.prev_row_copy();
        drow.next_row();
        srow.next_row();
      }
      SrcAccessT scol = splane.next_col_copy();
      DstAccessT dcol = dplane.next_col_copy();
      for( int32 y=1; y<dst.cols(); ++y ) {
        *dcol = *scol + *dcol.prev_col_copy();
        dcol.next_col();
        scol.next_col();
      }
    }

    for( int32 p=0; p<dst.planes(); ++p ) {
      SrcAccessT srow = splane.advance_copy(1,1);
      DstAccessT drow = dplane.advance_copy(1,1);
      for( int32 y=1; y<dst.rows(); ++y ) {
        SrcAccessT scol = srow;
        DstAccessT dcol = drow;
        for( int32 x=1; x<dst.cols(); ++x ) {
          *dcol = *scol + *dcol.advance_copy(-1,0) + *dcol.advance_copy(0,-1) - *dcol.advance_copy(-1,-1);
          scol.next_col();
          dcol.next_col();
        }
        srow.next_row();
        drow.next_row();
      }
      splane.next_plane();
      dplane.next_plane();
    }

    return dst;
  }

  // BlockSumView
  //
  // Applies a view where every pixel of the result represents the
  // block sum of the input. This is equivalent to convolution but much quicker
  // but also only for constant value square kernels.
  template <class ImageT>
  class BlockSumView : public ImageViewBase<BlockSumView<ImageT> > {
    ImageT m_child;
    uint32 m_range;

    // This type can hold image pixels and sum them with minimal risk of overflow
    typedef typename CompoundChannelCast<typename ImageT::pixel_type,typename AccumulatorType<typename CompoundChannelType<typename ImageT::pixel_type>::type>::type>::type sum_t;

    template <typename PixelAccessorT>
    inline sum_t sum_block(const PixelAccessorT& upper_left) const {
      return
          *upper_left
        + *upper_left.advance_copy(2*m_range, 2*m_range)
        - *upper_left.advance_copy(0, 2*m_range)
        - *upper_left.advance_copy(2*m_range, 0);
    }

  public:

    typedef typename ImageT::pixel_type pixel_type;
    typedef pixel_type result_type;
    typedef typename CompoundChannelType<pixel_type>::type channel_type;
    typedef ProceduralPixelAccessor<BlockSumView> pixel_accessor;

    BlockSumView( ImageT const& image, uint32 window_size = 3)
      : m_child(image), m_range((window_size-1)/2) {
        VW_ASSERT(window_size % 2 == 1, vw::LogicErr() << "BlockSumView's window_size must be odd");
    }

    inline int32 cols() const { return m_child.cols(); }
    inline int32 rows() const { return m_child.rows(); }
    inline int32 planes() const { return m_child.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(*this); }

    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const {
      ImageView<result_type> img(1,1,m_child.planes());
      this->rasterize(img, BBox2i(i,j,1,1));
      return img(0,0,p);
    }

    typedef CropView<ImageView<pixel_type> > prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      ImageView<pixel_type> dst(bbox.width(), bbox.height(), m_child.planes());
      rasterize( dst, bbox );
      return crop(dst, BBox2i(-bbox.min().x(),-bbox.min().y(),
                              m_child.cols(), m_child.rows()) );
    }

    template <class DstT>
    inline void rasterize( DstT const& dst, BBox2i const& bbox ) const {

      BBox2i child_bbox(bbox);
      child_bbox.expand(m_range);

      typedef ImageView<typename CompoundChannelCast<pixel_type,typename AccumulatorType<typename CompoundChannelType<pixel_type>::type>::type>::type> SrcT;

      SrcT src = edge_extend(integral_image(m_child), child_bbox, IntegralEdgeExtension());

      typedef typename SrcT::pixel_accessor SrcAccessT;
      typedef typename DstT::pixel_accessor DstAccessT;

      SrcAccessT splane = src.origin();
      DstAccessT dplane = dst.origin();
      for( int32 p=0; p<dst.planes(); ++p ) {
        SrcAccessT srow = splane;
        DstAccessT drow = dplane;
        for( int32 y=0; y<dst.rows(); ++y ) {
          SrcAccessT scol = srow;
          DstAccessT dcol = drow;
          for( int32 x=0; x<dst.cols(); ++x ) {
            *dcol = sum_block(scol);
            scol.next_col();
            dcol.next_col();
          }
          srow.next_row();
          drow.next_row();
        }
        splane.next_plane();
        dplane.next_plane();
      }

    }
  };

  template <class ImageT>
  BlockSumView<ImageT> block_sum( const ImageViewBase<ImageT>& img, uint32 window_size) {
    return BlockSumView<ImageT>( img.impl(), window_size );
  }

} // end namespace vw

#endif // __VW_IMAGE_INTEGRAL_VIEW_H__
