// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __VW_STEREO_REWRITE_SUBPIXEL_VIEW_H__
#define __VW_STEREO_REWRITE_SUBPIXEL_VIEW_H__

#include <vw/Image/Algorithms.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Stereo/rewrite/PreFilter.h>
#include <vw/Stereo/rewrite/CostFunctions.h>
#include <vw/Stereo/rewrite/Correlation.h>

#include <boost/foreach.hpp>

namespace vw {
namespace stereo {
namespace rewrite {

  // is_equal_elements, tells if all the elements are equal to each other
  template <class T>
  struct IsEqualElements {
    bool result;

    IsEqualElements() : result(true), m_comparing(false) {}

    void operator()( T const& i ) {
      if ( m_comparing ) {
        result = result && (*m_first_element == i);
      } else {
        m_comparing = true;
        m_first_element = &i;
      }
    }
  protected:
    bool m_comparing;
    T const* m_first_element;
  };

  template <class DImageT, class Image1T, class Image2T, class PreFilterT>
  class ParabolaSubpixelView : public ImageViewBase<ParabolaSubpixelView<DImageT,Image1T,Image2T,PreFilterT> > {
    DImageT m_disparity;
    Image1T m_left_image;
    Image2T m_right_image;
    PreFilterT m_prefilter;
    Vector2i m_kernel_size;

    Matrix<float,6,9> m_p_A_matrix;

    template <class FImage1T, class FImage2T>
    ImageView<PixelMask<Vector2f> >
    evaluate( ImageView<PixelMask<Vector2i> > const& integer_disparity,
              ImageViewBase<FImage1T> const& left_filter,
              ImageViewBase<FImage2T> const& right_filter,
              BBox2i const& left_region, BBox2i const& right_region,
              BBox2i const& disparity_region, BBox2i const& search_range ) const {
      ImageView<typename FImage1T::pixel_type> left_raster =
        crop(left_filter,left_region);
      ImageView<typename FImage2T::pixel_type> right_raster =
        crop(right_filter,right_region);

      // This will use 2.25 MB for a 256^2 pixel region
      ImageView<Vector<float,9> > cost_patch( integer_disparity.cols(),
                                              integer_disparity.rows() );
      typedef AbsoluteCost<FImage1T,boost::is_integral<typename PixelChannelType<typename FImage1T::pixel_type>::type>::value> CostType;
      typedef typename CostType::accumulator_type AccumChannelT;
      typedef typename PixelChannelCast<typename FImage1T::pixel_type,AccumChannelT>::type AccumT;
      CostType cost_function( left_raster, right_raster, m_kernel_size );
      typedef typename ImageView<Vector<float,9> >::pixel_accessor PatchAcc;
      typedef typename ImageView<AccumT>::pixel_accessor MetricAcc;
      typedef typename ImageView<PixelMask<Vector2i> >::pixel_accessor IDispAcc;

      // Subdivide the input disparity into smaller sections that have
      // the same disparity.
      std::list<SearchParam> zones;
      subdivide_regions( integer_disparity,
                         bounding_box(integer_disparity),
                         zones, m_kernel_size );
      BOOST_FOREACH( SearchParam& zone, zones ) {
        zone.second.expand(1);

        ImageView<AccumT> cost_metric( zone.first.width(), zone.first.height() );

        BBox2i left_zone = zone.first;
        left_zone.max() += m_kernel_size - Vector2i(1,1);

        for ( int32 dx = 0; dx < zone.second.width(); ++dx ) {
          for ( int32 dy = 0; dy < zone.second.height(); ++dy ) {
            Vector2i disparity( dx, dy );

            cost_metric =
              fast_box_sum<AccumChannelT>(cost_function( crop(left_raster,left_zone),
                                                         crop(right_raster,left_zone+disparity+zone.second.min()-search_range.min())),
                                          m_kernel_size );
            cost_function.cost_modification( cost_metric, disparity );

            VW_DEBUG_ASSERT( zone.first.width() == cost_metric.cols() &&
                             zone.first.height() == cost_metric.rows(),
                             MathErr() << "Cost Metric seems to have wrong size" );

            Vector2i disparity_adj = disparity + zone.second.min();

            // Now take the fast evaluate cost function and see where
            // the cost patch they apply
            PatchAcc patch_row = cost_patch.origin();
            patch_row.advance( zone.first.min().x(), zone.first.min().y() );
            MetricAcc metric_row = cost_metric.origin();
            IDispAcc idisp_row = integer_disparity.origin();
            idisp_row.advance( zone.first.min().x(), zone.first.min().y() );
            for ( int32 j = zone.first.height(); j; --j ) {
              PatchAcc patch_col = patch_row;
              MetricAcc metric_col = metric_row;
              IDispAcc idisp_col = idisp_row;
              for ( int32 i = zone.first.width(); i; --i ) {
                Vector2i delta = disparity_adj - (*idisp_col).child();
                if ( delta == Vector2i(-1,-1) ) {
                  (*patch_col)[0] = *metric_col;
                } else if ( delta == Vector2i( 0, -1 ) ) {
                  (*patch_col)[3] = *metric_col;
                } else if ( delta == Vector2i( 1, -1 ) ) {
                  (*patch_col)[6] = *metric_col;
                } else if ( delta == Vector2i( -1, 0 ) ) {
                  (*patch_col)[1] = *metric_col;
                } else if ( delta == Vector2i( 0, 0 ) ) {
                  (*patch_col)[4] = *metric_col;
                } else if ( delta == Vector2i( 1, 0 ) ) {
                  (*patch_col)[7] = *metric_col;
                } else if ( delta == Vector2i( -1, 1 ) ) {
                  (*patch_col)[2] = *metric_col;
                } else if ( delta == Vector2i( 0, 1 ) ) {
                  (*patch_col)[5] = *metric_col;
                } else if ( delta == Vector2i( 1, 1 ) ) {
                  (*patch_col)[8] = *metric_col;
                }
                patch_col.next_col();
                metric_col.next_col();
                idisp_col.next_col();
              }
              patch_row.next_row();
              metric_row.next_row();
              idisp_row.next_row();
            }
          }
        } // end disparity search
      } // end zone search

      // Calculate the new floating point location
      ImageView<PixelMask<Vector2f> > result( disparity_region.width(),
                                              disparity_region.height() );
      typedef typename ImageView<PixelMask<Vector2f> >::pixel_accessor ResultAcc;
      ResultAcc result_row = result.origin();
      IDispAcc idisp_row = integer_disparity.origin();
      PatchAcc patch_row = cost_patch.origin();
      for ( int32 j = disparity_region.height(); j; --j ) {
        ResultAcc result_col = result_row;
        IDispAcc idisp_col = idisp_row;
        PatchAcc patch_col = patch_row;
        for ( int32 i = disparity_region.width(); i; --i ) {
          // Given our 9 points of cost around our disparity pixel
          // (patch_col), find the 2D minimum and that will be our
          // floating point update.
          if ( is_valid( *idisp_col ) ) {
            IsEqualElements<float> is_same =
              std::for_each( (*patch_col).begin(), (*patch_col).end(),
                             IsEqualElements<float>() );
            if ( !is_same.result ) {
              // First, compute the parameters of the hyperbolic surface by
              // fitting the nine points in 'points' using a linear least
              // squares fit.  This process is fairly fast, since we have
              // already pre-computed the inverse of the A matrix in Ax = b.
              Vector<float,6> x( m_p_A_matrix * (*patch_col) );
              // With these parameters, we have a closed form expression for
              // the surface.  We compute the derivative, and find the point
              // where the slope is zero.  This is our maximum.
              //
              //  Max is at [x,y] where:
              //   dz/dx = 2ax + cy + d = 0
              //   dz/dy = 2by + cx + e = 0
              float denom = 4 * x[0] * x[1] - ( x[2] * x[2] );
              Vector2f offset( ( x[2] * x[4] - 2 * x[1] * x[3] ) / denom,
                               ( x[2] * x[3] - 2 * x[0] * x[4] ) / denom );
              if ( norm_2(offset) < 5.0 )
                *result_col = PixelMask<Vector2f>( remove_mask(*idisp_col) + offset );
              else
                *result_col = PixelMask<Vector2f>(*idisp_col);
            } else {
              *result_col = PixelMask<Vector2f>(*idisp_col);
            }
          } else {
            *result_col = PixelMask<Vector2f>();
          }
          result_col.next_col();
          idisp_col.next_col();
          patch_col.next_col();
        }
        result_row.next_row();
        idisp_row.next_row();
        patch_row.next_row();
      }

      return result;
    }

  public:
    typedef PixelMask<Vector2f> pixel_type;
    typedef pixel_type result_type;
    typedef ProceduralPixelAccessor<ParabolaSubpixelView> pixel_accessor;

    ParabolaSubpixelView(ImageViewBase<DImageT> const& disparity,
                         ImageViewBase<Image1T> const& left_image,
                         ImageViewBase<Image2T> const& right_image,
                         PreFilterT const& prefilter,
                         Vector2i const& kernel_size ) :
      m_disparity( disparity.impl() ), m_left_image( left_image.impl() ),
      m_right_image( right_image.impl() ), m_prefilter( prefilter ),
      m_kernel_size( kernel_size ) {
      VW_ASSERT( m_disparity.cols() == m_left_image.cols() &&
                 m_disparity.rows() == m_left_image.rows(),
                 ArgumentErr() << "SubpixelView: Disparity image must match left image." );

      // We get a considerable speedup in our 2d subpixel correlation if
      // we go ahead and compute the pseudoinverse of the A matrix (where
      // each row in A is [ x^2 y^2 xy x y 1] (our 2d parabolic surface)
      // for the range of x = [-1:1] and y = [-1:1].
      static float pinvA_data[] =
        { 1.0/6,  1.0/6,  1.0/6, -1.0/3, -1.0/3, -1.0/3,  1.0/6,  1.0/6,  1.0/6,
          1.0/6, -1.0/3,  1.0/6,  1.0/6, -1.0/3,  1.0/6,  1.0/6, -1.0/3,  1.0/6,
          1.0/4,    0.0, -1.0/4,    0.0,    0.0,    0.0, -1.0/4,    0.0,  1.0/4,
          -1.0/6, -1.0/6, -1.0/6,    0.0,    0.0,   0.0,  1.0/6,  1.0/6,  1.0/6,
          -1.0/6,    0.0,  1.0/6, -1.0/6,    0.0, 1.0/6, -1.0/6,    0.0,  1.0/6,
          -1.0/9,  2.0/9, -1.0/9,  2.0/9,  5.0/9, 2.0/9, -1.0/9,  2.0/9, -1.0/9 };
      m_p_A_matrix = Matrix<float,6,9>( pinvA_data );
    }

    inline int32 cols() const { return m_disparity.cols(); }
    inline int32 rows() const { return m_disparity.rows(); }
    inline int32 planes() const { return 1; }

    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }
    inline pixel_type operator() ( int32 /*i*/, int32 /*j*/, int32 /*p*/ = 0 ) const {
      vw_throw( NoImplErr() << "SubpixelView:operator() has not been implemented." );
      return pixel_type();
    }

    // Block rasterization section does the actual work
    typedef CropView<ImageView<pixel_type> > prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {

      // Prerasterize the disparity for this crop and work out it's
      // range. This could be faster if we subdivided this down by the
      // areas with equal values.
      ImageView<PixelMask<Vector2i> > disparity_subsection =
        crop(m_disparity,bbox);

      // We calculate the entire search range only so that we can
      // prerasterize the sections we'll need for left and right
      // images.
      BBox2i entire_search_range = stereo::get_disparity_range(disparity_subsection);
      entire_search_range.max() += Vector2i(1,1);
      entire_search_range.expand(1); // Because we are going to check
                                     // the neighboring pixels

      // Calculate our left an right regions that we need
      Vector2i half_kernel = m_kernel_size/2;
      BBox2i left_region = bbox;
      left_region.min() -= half_kernel;
      left_region.max() += half_kernel;
      BBox2i right_region = left_region + entire_search_range.min();
      right_region.max() += entire_search_range.size();

      // Rasterize the result with the prefilter
      return prerasterize_type(evaluate(disparity_subsection,
                                        m_prefilter.filter(m_left_image),
                                        m_prefilter.filter(m_right_image),
                                        left_region, right_region, bbox,
                                        entire_search_range ),
                               -bbox.min().x(), -bbox.min().y(), cols(), rows());
    }

    template <class DestT>
    inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
      vw::rasterize(prerasterize(bbox), dest, bbox );
    }
  };

  template <class DImageT, class Image1T, class Image2T, class PreFilterT>
  ParabolaSubpixelView<DImageT, Image1T, Image2T, PreFilterT>
  parabola_subpixel( ImageViewBase<DImageT> const& disparity,
                     ImageViewBase<Image1T> const& left_image,
                     ImageViewBase<Image2T> const& right_image,
                     PreFilterT const& prefilter,
                     Vector2i const& kernel_size ) {
    typedef ParabolaSubpixelView<DImageT, Image1T, Image2T, PreFilterT> result_type;
    return result_type( disparity.impl(), left_image.impl(), right_image.impl(),
                        prefilter, kernel_size );
  }

}}} // end namespace vw::stereo::rewrite

#endif//__VW_STEREO_REWRITE_SUBPIXEL_VIEW_H__
