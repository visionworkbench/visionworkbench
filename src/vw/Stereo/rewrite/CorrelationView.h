// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __VW_STEREO_REWRITE_CORRELATION_VIEW_H__
#define __VW_STEREO_REWRITE_CORRELATION_VIEW_H__

#include <vw/Image/Algorithms.h>
#include <vw/Stereo/rewrite/Correlation.h>

namespace vw {
namespace stereo {
namespace rewrite {

  /// An image view for performing image correlation
  template <class Image1T, class Image2T, class PreFilterT>
  class CorrelationView : public ImageViewBase<CorrelationView<Image1T, Image2T, PreFilterT> > {

    Image1T m_left_image;
    Image2T m_right_image;
    PreFilterT m_prefilter;
    BBox2i m_search_region;
    Vector2i m_kernel_size;
    rewrite::CostFunctionType m_cost_type;
    float m_consistency_threshold; // 0 = means don't do a consistency check

  public:
    typedef PixelMask<Vector2i> pixel_type;
    typedef PixelMask<Vector2i> result_type;
    typedef ProceduralPixelAccessor<CorrelationView> pixel_accessor;

    CorrelationView( ImageViewBase<Image1T> const& left,
                     ImageViewBase<Image2T> const& right,
                     PreFilterT const& prefilter,
                     BBox2i const& search_region, Vector2i const& kernel_size,
                     CostFunctionType cost_type = ABSOLUTE_DIFFERENCE,
                     float consistency_threshold = -1 ) :
      m_left_image(left.impl()), m_right_image(right.impl()),
      m_prefilter(prefilter), m_search_region(search_region), m_kernel_size(kernel_size),
      m_cost_type(cost_type), m_consistency_threshold(consistency_threshold) {}

    // Standard required ImageView interfaces
    inline int32 cols() const { return m_left_image.cols(); }
    inline int32 rows() const { return m_left_image.rows(); }
    inline int32 planes() const { return 1; }

    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }
    inline pixel_type operator()( int32 /*i*/, int32 /*j*/, int32 /*p*/ = 0) const {
      vw_throw( NoImplErr() << "CorrelationView::operator()(....) has not been implemented." );
      return pixel_type();
    }

    // Block rasterization section that does actual work
    typedef CropView<ImageView<pixel_type> > prerasterize_type;
    inline prerasterize_type prerasterize(BBox2i const& bbox) const {

      // 1.) Expand the left raster region by the kernel size.
      Vector2i half_kernel = m_kernel_size/2;
      BBox2i left_region = bbox;
      left_region.min() -= half_kernel;
      left_region.max() += half_kernel;

      // 2.) Calculate the region of the right image that we're using.
      BBox2i right_region = left_region + m_search_region.min();
      right_region.max() += m_search_region.size();

      // 3.) Correlate with options that they requested
      ImageView<pixel_type> result;
      // Shutting off the consistency check
      switch ( m_cost_type ) {
      case CROSS_CORRELATION:
        result =
          best_of_search_convolution<NCCCost>( m_prefilter.filter(crop(edge_extend(m_left_image),left_region)),
                                               m_prefilter.filter(crop(edge_extend(m_right_image),right_region)),
                                               left_region - left_region.min(),
                                               m_search_region.size() + Vector2i(1,1),
                                               m_kernel_size );
        break;
      case SQUARED_DIFFERENCE:
        result =
          best_of_search_convolution<SquaredCost>( m_prefilter.filter(crop(edge_extend(m_left_image),left_region)),
                                                   m_prefilter.filter(crop(edge_extend(m_right_image),right_region)),
                                                   left_region - left_region.min(),
                                                   m_search_region.size() + Vector2i(1,1),
                                                   m_kernel_size );
        break;
      case ABSOLUTE_DIFFERENCE:
      default:
        result =
          best_of_search_convolution<AbsoluteCost>( m_prefilter.filter(crop(edge_extend(m_left_image),left_region)),
                                                    m_prefilter.filter(crop(edge_extend(m_right_image),right_region)),
                                                    left_region - left_region.min(),
                                                    m_search_region.size() + Vector2i(1,1),
                                                    m_kernel_size );
      }

      // 4.0 ) Do a consistency check if they asked for it
      if ( m_consistency_threshold >= 0 ) {
        ImageView<pixel_type> rl_result;

        switch ( m_cost_type ) {
        case CROSS_CORRELATION:
          // Getting the crops correctly here is not important as best
          // of search convolution will recrop. The important bit is
          // just aligning up the origins.
          rl_result =
            best_of_search_convolution<NCCCost>( m_prefilter.filter(crop(edge_extend(m_right_image),right_region)),
                                                 m_prefilter.filter(crop(edge_extend(m_left_image),left_region-(m_search_region.size()+Vector2i(1,1)))),
                                                 right_region - right_region.min(),
                                                 m_search_region.size() + Vector2i(1,1),
                                                 m_kernel_size ) -
            PixelMask<Vector2i>(m_search_region.size()+Vector2i(1,1));
          break;
        case SQUARED_DIFFERENCE:
          // Getting the crops correctly here is not important as best
          // of search convolution will recrop. The important bit is
          // just aligning up the origins.
          rl_result =
            best_of_search_convolution<SquaredCost>( m_prefilter.filter(crop(edge_extend(m_right_image),right_region)),
                                                     m_prefilter.filter(crop(edge_extend(m_left_image),left_region-(m_search_region.size()+Vector2i(1,1)))),
                                                     right_region - right_region.min(),
                                                     m_search_region.size() + Vector2i(1,1),
                                                     m_kernel_size ) -
            PixelMask<Vector2i>(m_search_region.size()+Vector2i(1,1));
          break;
        case ABSOLUTE_DIFFERENCE:
        default:
          // Getting the crops correctly here is not important as best
          // of search convolution will recrop. The important bit is
          // just aligning up the origins.
          rl_result =
            best_of_search_convolution<AbsoluteCost>( m_prefilter.filter(crop(edge_extend(m_right_image),right_region)),
                                                      m_prefilter.filter(crop(edge_extend(m_left_image),left_region-(m_search_region.size()+Vector2i(1,1)))),
                                                      right_region - right_region.min(),
                                                      m_search_region.size() + Vector2i(1,1),
                                                      m_kernel_size ) -
            PixelMask<Vector2i>(m_search_region.size()+Vector2i(1,1));
        }
        stereo::cross_corr_consistency_check( result, rl_result,
                                              m_consistency_threshold, false );
      }
      VW_DEBUG_ASSERT( bbox.size() == bounding_box(result).size(),
                       MathErr() << "CorrelationView::prerasterize got a bad return from best_of_search_convolution." );

      // 5.) Convert back to original coordinates
      result += pixel_type(m_search_region.min());
      return prerasterize_type( result, -bbox.min().x(), -bbox.min().y(), cols(), rows() );
    }

    template <class DestT>
    inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
      vw::rasterize(prerasterize(bbox), dest, bbox);
    }
  };

  template <class Image1T, class Image2T, class PreFilterT>
  CorrelationView<Image1T,Image2T,PreFilterT>
  correlate( ImageViewBase<Image1T> const& left,
             ImageViewBase<Image2T> const& right,
             PreFilterT const& filter,
             BBox2i const& search_region, Vector2i const& kernel_size,
             CostFunctionType cost_type = ABSOLUTE_DIFFERENCE,
             float consistency_threshold = -1 ) {
    typedef CorrelationView<Image1T,Image2T,PreFilterT> result_type;
    return result_type( left.impl(), right.impl(), filter, search_region,
                        kernel_size, cost_type, consistency_threshold );
  }
}}}

#endif//__VW_STEREO_REWRITE_CORRELATION_VIEW_H__
