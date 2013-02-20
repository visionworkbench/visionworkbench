// __BEGIN_LICENSE__
//  Copyright (c) 2006-2012, United States Government as represented by the
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


#ifndef __VW_STEREO_SUBPIXEL_VIEW__
#define __VW_STEREO_SUBPIXEL_VIEW__

#include <vw/Image.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Stereo/Correlate.h>
#include <vw/Stereo/PreFilter.h>
#include <vw/Stereo/CostFunctions.h>
#include <vw/Stereo/Correlation.h>

#include <boost/foreach.hpp>

namespace vw {
namespace stereo {

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
            if ( !std::equal( (*patch_col).begin()+1, (*patch_col).end(),
                              (*patch_col).begin() ) ) {
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
                         PreFilterBase<PreFilterT> const& prefilter,
                         Vector2i const& kernel_size ) :
      m_disparity( disparity.impl() ), m_left_image( left_image.impl() ),
      m_right_image( right_image.impl() ), m_prefilter( prefilter.impl() ),
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
                     PreFilterBase<PreFilterT> const& prefilter,
                     Vector2i const& kernel_size ) {
    typedef ParabolaSubpixelView<DImageT, Image1T, Image2T, PreFilterT> result_type;
    return result_type( disparity.impl(), left_image.impl(), right_image.impl(),
                        prefilter.impl(), kernel_size );
  }

  /// An image view for performing image correlation
  template <class PreprocFilterT, class ImageT1, class ImageT2, class ImageTD>
  class BayesEMSubpixelView : public ImageViewBase<BayesEMSubpixelView<PreprocFilterT, ImageT1, ImageT2, ImageTD> > {

    ImageTD m_disparity_map;
    ImageT1 m_left_image;
    ImageT2 m_right_image;

    // General Settings
    Vector2i m_kernel_size;
    PreprocFilterT m_preproc_filter;
    int32 m_max_pyramid_levels;

  public:
    typedef PixelMask<Vector2f> pixel_type;
    typedef pixel_type result_type;
    typedef ProceduralPixelAccessor<BayesEMSubpixelView> pixel_accessor;

    BayesEMSubpixelView(ImageViewBase<ImageTD> const& disparity_map,
                        ImageViewBase<ImageT1> const& left_image,
                        ImageViewBase<ImageT2> const& right_image,
                        PreFilterBase<PreprocFilterT> const& preproc_filter,
                        Vector2i const& kernel_size,
                        int32 max_pyramid_levels ) :
      m_disparity_map(disparity_map.impl()),
      m_left_image(left_image.impl()), m_right_image(right_image.impl()),
      m_kernel_size(kernel_size), m_preproc_filter(preproc_filter.impl()),
      m_max_pyramid_levels(max_pyramid_levels) {

      // Basic assertions
      VW_ASSERT( m_disparity_map.cols() == m_left_image.cols() &&
                 m_disparity_map.rows() == m_left_image.rows(),
                 ArgumentErr() << "BayesEMSubpixelView::BayesEMSubpixelView():  Disparity image must match left image.\n");

      VW_ASSERT((m_left_image.channels() == 1) && (m_left_image.planes() == 1) &&
                (m_right_image.channels() == 1) && (m_right_image.planes() == 1),
                ArgumentErr() << "BayesEMSubpixelView::BayesEMSubpixelView(): multi-channel, multi-plane images not supported.\n");

      // Max pyramid can't go below 0 ... or bayes em won't process anything
      if ( m_max_pyramid_levels < 0 ) m_max_pyramid_levels = 0;
    }

    // Standard ImageView interface methods
    inline int32 cols() const { return m_left_image.cols(); }
    inline int32 rows() const { return m_left_image.rows(); }
    inline int32 planes() const { return 1; }

    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

    inline pixel_type operator()(float /*x*/, float /*y*/, int32 /*p*/ = 0) const {
      vw_throw(NoImplErr() << "BayesEMSubpixelView::operator() is not yet implemented.");
      return PixelMask<Vector2f>(); // Never reached
    }


    /// \cond INTERNAL
    typedef CropView<ImageView<pixel_type> > prerasterize_type;
    inline prerasterize_type prerasterize(BBox2i const& bbox) const {

#if VW_DEBUG_LEVEL > 0
      Stopwatch watch;
      watch.start();
#endif
      // Find the range of disparity values for this patch.
      ImageView<pixel_type > disparity_map_patch =
        crop(m_disparity_map, bbox);
      BBox2i search_range = get_disparity_range(disparity_map_patch);

      // The area in the right image that we'll be searching is
      // determined by the bbox of the left image plus the search
      // range.
      BBox2i left_crop_bbox(bbox);
      BBox2i right_crop_bbox(bbox.min() + search_range.min(),
                             bbox.max() + search_range.max());

      // The correlator requires the images to be the same size. The
      // search bbox will always be larger than the given left image
      // bbox, so we just make the left bbox the same size as the
      // right bbox.
      left_crop_bbox.max() = left_crop_bbox.min() + right_crop_bbox.size();

      // Finally, we must adjust both bounding boxes to account for
      // the size of the kernel itself.
      right_crop_bbox.min() -= m_kernel_size;
      right_crop_bbox.max() += m_kernel_size;
      left_crop_bbox.min() -=  m_kernel_size;
      left_crop_bbox.max() +=  m_kernel_size;

      // We crop the images to the expanded bounding box and edge
      // extend in case the new bbox extends past the image bounds.
      ImageView<float> left_image_patch, right_image_patch;
      // parabola subpixel does the same preprocessing as the pyramid correlator
      left_image_patch =
        crop(m_preproc_filter.filter(m_left_image), left_crop_bbox);
      right_image_patch =
        crop(m_preproc_filter.filter(m_right_image), right_crop_bbox);
      disparity_map_patch =
        crop(edge_extend(m_disparity_map, ZeroEdgeExtension()),
             left_crop_bbox);

      // Adjust the disparities to be relative to the cropped
      // image pixel locations
      PixelMask<Vector2f> disparity_patch_translation( search_range.min() );
      disparity_map_patch -= disparity_patch_translation;

      std::vector< ImageView<float> > l_patches, r_patches;
      std::vector<BBox2i> rois;
      ImageView<pixel_type > d_subpatch;

      // This initialization makes a difference only if the number of
      // levels is 0, which is a valid situation.
      d_subpatch = disparity_map_patch;

      // I'd like for image subsampling to use a gaussian when
      // downsampling however it was introducing some edge effects
      // that I couldn't figure out within a reasonable time frame.
      for ( int32 i = 0; i < m_max_pyramid_levels; i++ ) {
        if ( i > 0 ) {
          // Building all other levels
          l_patches.push_back( subsample( l_patches.back(), 2 ) );
          r_patches.push_back( subsample( r_patches.back(), 2 ) );
          ImageView<PixelMask<Vector2f> > d_subpatch_buf =
            disparity_subsample( d_subpatch );
          d_subpatch = d_subpatch_buf;
          rois.push_back( rois.back()/2 );
        } else {
          // First step down from native resolution
          l_patches.push_back( subsample( left_image_patch, 2 ) );
          r_patches.push_back( subsample( right_image_patch, 2 ) );
          d_subpatch = disparity_subsample( disparity_map_patch );
          rois.push_back( BBox2i( m_kernel_size[0], m_kernel_size[1],
                                  bbox.width(), bbox.height() ) );
        }
      }

      for ( int32 i = m_max_pyramid_levels-1; i >= 0; i-- ) {
        subpixel_optimized_affine_2d_EM(d_subpatch,
                                        l_patches[i], r_patches[i],
                                        m_kernel_size[0], m_kernel_size[1],
                                        rois[i], true, true, false );
        BBox2i crop_bbox;
        if ( i > 0 )
          crop_bbox = BBox2i(0,0,l_patches[i-1].cols(),
                             l_patches[i-1].rows());
        else
          crop_bbox = BBox2i(0,0,left_image_patch.cols(),
                             left_image_patch.rows());
        ImageView<pixel_type > d_subpatch_buf =
          crop(disparity_upsample(edge_extend(d_subpatch)), crop_bbox);
        d_subpatch = d_subpatch_buf;
      }

      disparity_map_patch = d_subpatch;

      // Perfrom final pass at native resolution
      subpixel_optimized_affine_2d_EM(disparity_map_patch,
                                      left_image_patch,
                                      right_image_patch,
                                      m_kernel_size[0], m_kernel_size[1],
                                      BBox2i(m_kernel_size[0],m_kernel_size[1],
                                             bbox.width(), bbox.height()),
                                      true, true, false );

      // Undo the above adjustment
      disparity_map_patch += disparity_patch_translation;

#if VW_DEBUG_LEVEL > 0
      watch.stop();
      vw_out(DebugMessage,"stereo") << "Tile " << bbox << " processed in " << watch.elapsed_seconds() << " s\n";
#endif

      // This may seem confusing, but we must crop here so that the
      // good pixel data is placed into the coordinates specified by
      // the bbox.  This allows rasterize to touch those pixels
      // using the coordinates inside the bbox.  The pixels outside
      // those coordinates are invalid, and they never get accessed.
      return crop(disparity_map_patch, BBox2i(m_kernel_size[0]-bbox.min()[0],
                                              m_kernel_size[1]-bbox.min()[1],
                                              m_left_image.cols(),
                                              m_left_image.rows() ));
    }

    template <class DestT> inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
      vw::rasterize(prerasterize(bbox), dest, bbox);
    }
    /// \endcond
  };

  template <class PreprocFilterT, class ImageT1, class ImageT2, class DisparityT>
  BayesEMSubpixelView<PreprocFilterT, ImageT1, ImageT2, DisparityT>
  bayes_em_subpixel( ImageViewBase<DisparityT> const& disparity_map,
                     ImageViewBase<ImageT1> const& left_image,
                     ImageViewBase<ImageT2> const& right_image,
                     PreFilterBase<PreprocFilterT> const& filter,
                     Vector2i const& kernel_size,
                     int max_pyramid_levels = 2 ) {
    typedef BayesEMSubpixelView<PreprocFilterT,ImageT1,ImageT2,DisparityT> result_type;
    return result_type( disparity_map.impl(), left_image.impl(),
                        right_image.impl(), filter.impl(), kernel_size,
                        max_pyramid_levels );
  }

}} // namespace vw::stereo

#endif // __VW_STEREO_CORRELATOR_VIEW__
