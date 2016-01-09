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


#ifndef __VW_STEREO_CORRELATION_VIEW_H__
#define __VW_STEREO_CORRELATION_VIEW_H__

#include <vw/Core/Exception.h>
#include <vw/Core/Stopwatch.h>
#include <vw/Core/Thread.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/PerPixelAccessorViews.h>
#include <vw/FileIO.h>
#include <vw/Stereo/Correlation.h>
#include <vw/Stereo/Correlate.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Stereo/PreFilter.h>
#include <boost/foreach.hpp>
#include <ctime>

namespace vw {
namespace stereo {

  /// An image view for performing image correlation
  /// - For each left image pixel, compute disparity vector to matching pixel in the right image.
  template <class Image1T, class Image2T, class PreFilterT>
  class CorrelationView : public ImageViewBase<CorrelationView<Image1T, Image2T, PreFilterT> > {

    Image1T          m_left_image;
    Image2T          m_right_image;
    PreFilterT       m_prefilter;
    BBox2i           m_search_region;
    Vector2i         m_kernel_size;
    CostFunctionType m_cost_type;
    float            m_consistency_threshold; // 0 = means don't do a consistency check

  public:
    typedef PixelMask<Vector2i> pixel_type;
    typedef PixelMask<Vector2i> result_type;
    typedef ProceduralPixelAccessor<CorrelationView> pixel_accessor;

    CorrelationView( ImageViewBase<Image1T>    const& left,
                     ImageViewBase<Image2T>    const& right,
                     PreFilterBase<PreFilterT> const& prefilter,
                     BBox2i const& search_region, Vector2i const& kernel_size,
                     CostFunctionType cost_type = ABSOLUTE_DIFFERENCE,
                     float            consistency_threshold = -1 ) :
      m_left_image(left.impl()), m_right_image(right.impl()),
      m_prefilter(prefilter.impl()), m_search_region(search_region), m_kernel_size(kernel_size),
      m_cost_type(cost_type), m_consistency_threshold(consistency_threshold) {}

    // Standard required ImageView interfaces
    inline int32 cols  () const { return m_left_image.cols(); }
    inline int32 rows  () const { return m_left_image.rows(); }
    inline int32 planes() const { return 1; }

    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }
    inline pixel_type operator()( int32 /*i*/, int32 /*j*/, int32 /*p*/ = 0) const {
      vw_throw( NoImplErr() << "CorrelationView::operator()(....) has not been implemented." );
      return pixel_type();
    }

    // Block rasterization section that does actual work
    typedef CropView<ImageView<pixel_type> > prerasterize_type;
    inline prerasterize_type prerasterize(BBox2i const& bbox) const {

#if VW_DEBUG_LEVEL > 0
      Stopwatch watch;
      watch.start();
#endif

      // 1.) Expand the left raster region by the kernel size.
      Vector2i half_kernel = m_kernel_size/2;
      BBox2i  left_region  = bbox;
      left_region.min() -= half_kernel;
      left_region.max() += half_kernel;

      // 2.) Calculate the region of the right image that we're using.
      BBox2i right_region = left_region + m_search_region.min();
      right_region.max() += m_search_region.size();

      // 3.) Calculate the disparity
      ImageView<pixel_type> result
        = calc_disparity(m_cost_type,
                         crop(m_prefilter.filter(m_left_image),left_region),
                         crop(m_prefilter.filter(m_right_image),right_region),
                         left_region - left_region.min(),
                         m_search_region.size() + Vector2i(1,1),
                         m_kernel_size);

      // 4.0 ) Consistency check
      if ( m_consistency_threshold >= 0 ) {
        // Getting the crops correctly here is not important as we
        // will re-crop later. The important bit is aligning up the origins.
        ImageView<pixel_type> rl_result
          = calc_disparity(m_cost_type,
                           crop(m_prefilter.filter(m_right_image),right_region),
                           crop(m_prefilter.filter(m_left_image),
                                left_region - (m_search_region.size()+Vector2i(1,1))),
                           right_region - right_region.min(),
                           m_search_region.size() + Vector2i(1,1),
                           m_kernel_size) -
          pixel_type(m_search_region.size()+Vector2i(1,1));

        stereo::cross_corr_consistency_check( result, rl_result,
                                              m_consistency_threshold, false );
      }
      VW_ASSERT( bbox.size() == bounding_box(result).size(),
                 MathErr() << "CorrelationView::prerasterize got a bad return from best_of_search_convolution." );

      // 5.) Convert back to original coordinates
      result += pixel_type(m_search_region.min());

#if VW_DEBUG_LEVEL > 0
      watch.stop();
      vw_out(DebugMessage,"stereo") << "Tile " << bbox << " processed in " << watch.elapsed_seconds() << " s\n";
#endif

      return prerasterize_type( result, -bbox.min().x(), -bbox.min().y(), cols(), rows() );
    } // End function prerasterize

    template <class DestT>
    inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
      vw::rasterize(prerasterize(bbox), dest, bbox);
    }
  }; // End class CorrelationView

  template <class Image1T, class Image2T, class PreFilterT>
  CorrelationView<Image1T,Image2T,PreFilterT>
  correlate( ImageViewBase<Image1T> const& left,
             ImageViewBase<Image2T> const& right,
             PreFilterBase<PreFilterT> const& filter,
             BBox2i const& search_region, Vector2i const& kernel_size,
             CostFunctionType cost_type = ABSOLUTE_DIFFERENCE,
             float consistency_threshold = -1 ) {
    typedef CorrelationView<Image1T,Image2T,PreFilterT> result_type;
    return result_type( left.impl(), right.impl(), filter.impl(), search_region,
                        kernel_size, cost_type, consistency_threshold );
  }



  /// An image view for performing pyramid image correlation (Faster than CorrelationView).
  /// - TODO: What is prefilter?
  template <class Image1T, class Image2T, class Mask1T, class Mask2T, class PreFilterT>
  class PyramidCorrelationView : public ImageViewBase<PyramidCorrelationView<Image1T,Image2T, Mask1T, Mask2T,PreFilterT> > {

    Image1T          m_left_image;
    Image2T          m_right_image;
    Mask1T           m_left_mask;
    Mask2T           m_right_mask;
    PreFilterT       m_prefilter;
    BBox2i           m_search_region;
    Vector2i         m_kernel_size;
    CostFunctionType m_cost_type;
    int              m_corr_timeout;
    // How long it takes to do one corr op with given kernel and cost function
    double m_seconds_per_op;
    float  m_consistency_threshold; // < 0 = means don't do a consistency check
    int32  m_max_level_by_search;

    struct SubsampleMaskByTwoFunc : public ReturnFixedType<uint8> {
      BBox2i work_area() const { return BBox2i(0,0,2,2); }

      template <class PixelAccessorT>
      typename boost::remove_reference<typename PixelAccessorT::pixel_type>::type
      operator()( PixelAccessorT acc ) const {

        typedef typename PixelAccessorT::pixel_type PixelT;

        uint8 count = 0;
        if ( *acc ) count++;
        acc.next_col();
        if ( *acc ) count++;
        acc.advance(-1,1);
        if ( *acc ) count++;
        acc.next_col();
        if ( *acc ) count++;
        if ( count > 1 )
          return PixelT(ScalarTypeLimits<PixelT>::highest());
        return PixelT();
      }
    }; // End struct SubsampleMaskByTwoFunc

    template <class ViewT>
    SubsampleView<UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ZeroEdgeExtension>, SubsampleMaskByTwoFunc> >
    subsample_mask_by_two( ImageViewBase<ViewT> const& input ) const {
      return subsample(per_pixel_accessor_filter(input.impl(), SubsampleMaskByTwoFunc()),2);
    }

  public:
    typedef PixelMask<Vector2i> pixel_type;
    typedef PixelMask<Vector2i> result_type;
    typedef ProceduralPixelAccessor<PyramidCorrelationView> pixel_accessor;

    PyramidCorrelationView( ImageViewBase<Image1T> const& left,
                            ImageViewBase<Image2T> const& right,
                            ImageViewBase<Mask1T> const& left_mask,
                            ImageViewBase<Mask2T> const& right_mask,
                            PreFilterBase<PreFilterT> const& prefilter,
                            BBox2i const& search_region, Vector2i const& kernel_size,
                            CostFunctionType cost_type,
                            int corr_timeout, double seconds_per_op,
                            float consistency_threshold,
                            int32 max_pyramid_levels) :
      m_left_image(left.impl()), m_right_image(right.impl()),
      m_left_mask(left_mask.impl()), m_right_mask(right_mask.impl()),
      m_prefilter(prefilter.impl()), m_search_region(search_region), m_kernel_size(kernel_size),
      m_cost_type(cost_type),
      m_corr_timeout(corr_timeout), m_seconds_per_op(seconds_per_op),
      m_consistency_threshold(consistency_threshold){
      // Calculating max pyramid levels according to the supplied
      // search region.
      int32 largest_search = max( search_region.size() );
      m_max_level_by_search = std::floor(std::log(float(largest_search))/std::log(2.0f)) - 1;
      if ( m_max_level_by_search > max_pyramid_levels )
        m_max_level_by_search = max_pyramid_levels;
      if ( m_max_level_by_search < 0 )
        m_max_level_by_search = 0;
    }

    // Standard required ImageView interfaces
    inline int32 cols  () const { return m_left_image.cols(); }
    inline int32 rows  () const { return m_left_image.rows(); }
    inline int32 planes() const { return 1; }

    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }
    inline pixel_type operator()( int32 /*i*/, int32 /*j*/, int32 /*p*/ = 0) const {
      vw_throw( NoImplErr() << "PyramidCorrelationView::operator()(....) has not been implemented." );
      return pixel_type();
    }

    /// Block rasterization section that does actual work
    typedef CropView<ImageView<pixel_type> > prerasterize_type;
    inline prerasterize_type prerasterize(BBox2i const& bbox) const {

      time_t start, end;
      if (m_corr_timeout){
        std::time (&start);
      }

#if VW_DEBUG_LEVEL > 0
      Stopwatch watch;
      watch.start();
#endif

      // 1.0) Determining the number of levels to process
      //      There's a maximum base on kernel size. There's also
      //      maximum defined by the search range. Here we determine
      //      the maximum based on kernel size and current bbox.
      int32 smallest_bbox      = math::min(bbox.size());
      int32 largest_kernel     = math::max(m_kernel_size);
      int32 max_pyramid_levels = std::floor(log(smallest_bbox)/log(2.0f) - log(largest_kernel)/log(2.0f));
      if ( m_max_level_by_search < max_pyramid_levels )
        max_pyramid_levels = m_max_level_by_search;
      if ( max_pyramid_levels < 1 )
        max_pyramid_levels = 0;
      Vector2i half_kernel = m_kernel_size/2;

      // 2.0) Build the pyramids
      //      - Highest resolution image is stored at index zero.
      //      - There really ought to be a function call for this!
      std::vector<ImageView<typename Image1T::pixel_type> > left_pyramid      (max_pyramid_levels + 1 );
      std::vector<ImageView<typename Image2T::pixel_type> > right_pyramid     (max_pyramid_levels + 1 );
      std::vector<ImageView<typename Mask1T::pixel_type > > left_mask_pyramid (max_pyramid_levels + 1 );
      std::vector<ImageView<typename Mask2T::pixel_type > > right_mask_pyramid(max_pyramid_levels + 1 );
      int32 max_upscaling = 1 << max_pyramid_levels;
      BBox2i left_global_region, right_global_region;
      {
        left_global_region = bbox;
        left_global_region.min() -= half_kernel * max_upscaling;
        left_global_region.max() += half_kernel * max_upscaling;
        right_global_region = left_global_region + m_search_region.min();
        right_global_region.max() += m_search_region.size() + Vector2i(max_upscaling,max_upscaling);
        left_pyramid      [0] = crop(edge_extend(m_left_image), left_global_region);
        right_pyramid     [0] = crop(edge_extend(m_right_image),right_global_region);
        left_mask_pyramid [0] = crop(edge_extend(m_left_mask, ConstantEdgeExtension()), left_global_region);
        right_mask_pyramid[0] = crop(edge_extend(m_right_mask,ConstantEdgeExtension()), right_global_region);

#if VW_DEBUG_LEVEL > 0
        VW_OUT(DebugMessage,"stereo") << " > Left ROI: "    << left_global_region
                                      << "\n > Right ROI: " << right_global_region << "\n";
#endif

        // Fill in the nodata of the left and right images with a mean
        // pixel value. This helps with the edge quality of a DEM.
        typename Image1T::pixel_type left_mean;
        typename Image2T::pixel_type right_mean;
        try {
          left_mean  = mean_pixel_value(subsample(copy_mask(left_pyramid [0], create_mask(left_mask_pyramid [0],0)),2));
          right_mean = mean_pixel_value(subsample(copy_mask(right_pyramid[0], create_mask(right_mask_pyramid[0],0)),2));
        } catch ( const ArgumentErr& err ) {
          // Mean pixel value will throw an argument error if there
          // are no valid pixels. If that happens, it means either the
          // left or the right image is full masked.
#if VW_DEBUG_LEVEL > 0
          watch.stop();
          double elapsed = watch.elapsed_seconds();
          vw_out(DebugMessage,"stereo")
            << "Tile " << bbox << " has no data. Processed in "
            << elapsed << " s\n";
#endif
          return prerasterize_type(ImageView<pixel_type>(bbox.width(),
                                                         bbox.height()),
                                   -bbox.min().x(), -bbox.min().y(),
                                   cols(), rows() );
        } // End mean pixel insertion section
        left_pyramid [0] = apply_mask(copy_mask(left_pyramid [0],create_mask(left_mask_pyramid [0],0)), left_mean  );
        right_pyramid[0] = apply_mask(copy_mask(right_pyramid[0],create_mask(right_mask_pyramid[0],0)), right_mean );

        // Don't actually need the whole over cropped disparity
        // mask. We only need the active region. I over cropped before
        // just to calculate the mean color value options.
        BBox2i right_mask = bbox + m_search_region.min();
        right_mask.max() += m_search_region.size();
        left_mask_pyramid [0] = crop(left_mask_pyramid [0], bbox       - left_global_region.min());
        right_mask_pyramid[0] = crop(right_mask_pyramid[0], right_mask - right_global_region.min());

        // Szeliski's book recommended this simple kernel. This
        // operation is quickly becoming a time sink, we might
        // possibly want to write an integer optimized version.
        std::vector<typename DefaultKernelT<typename Image1T::pixel_type>::type > kernel(5);
        kernel[0] = kernel[4] = 1.0/16.0;
        kernel[1] = kernel[3] = 4.0/16.0;
        kernel[2] = 6.0/16.0;
        std::vector<uint8> mask_kern(max(m_kernel_size));
        std::fill(mask_kern.begin(), mask_kern.end(), 1 );

        // Build the pyramid first and then apply the filter to each level.
        for ( int32 i = 0; i < max_pyramid_levels; ++i ) {
          left_pyramid [i+1] = subsample(separable_convolution_filter(left_pyramid [i],kernel,kernel),2);
          right_pyramid[i+1] = subsample(separable_convolution_filter(right_pyramid[i],kernel,kernel),2);
          left_pyramid [i] = m_prefilter.filter(left_pyramid [i]);
          right_pyramid[i] = m_prefilter.filter(right_pyramid[i]);
          left_mask_pyramid [i+1] = subsample_mask_by_two(left_mask_pyramid [i]);
          right_mask_pyramid[i+1] = subsample_mask_by_two(right_mask_pyramid[i]);
        }
        left_pyramid [max_pyramid_levels] = m_prefilter.filter(left_pyramid [max_pyramid_levels]);
        right_pyramid[max_pyramid_levels] = m_prefilter.filter(right_pyramid[max_pyramid_levels]);
      } // Done building the pyramids!

      // 3.0) Actually perform correlation now
      ImageView<pixel_type > disparity;
      std::vector<SearchParam> zones; // Initial search region at lowest resolution level
      zones.push_back( SearchParam(bounding_box(left_mask_pyramid[max_pyramid_levels]),
                                   BBox2i(0,0,m_search_region.width ()/max_upscaling+1,
                                              m_search_region.height()/max_upscaling+1)) );

      // Perform correlation. Keep track of how much time elapsed
      // since we started and stop if we estimate that doing one more
      // image chunk will bring us over time.

      // To not slow us down with timing, we use some heuristics to
      // estimate how much time elapsed, as time to do an image chunk
      // is proportional with image area times search range area. This
      // is not completely accurate, so every now and then do actual
      // timing, no more often than once in measure_spacing seconds.
      double estim_elapsed = 0.0;
      int measure_spacing = 2; // seconds
      double prev_estim = estim_elapsed;

      // Loop down through all of the pyramid levels, low res to high res.
      for ( int32 level = max_pyramid_levels; level >= 0; --level) {

        int32 scaling = 1 << level;
        disparity.set_size( left_mask_pyramid[level] );
        Vector2i region_offset = max_upscaling*half_kernel/scaling;

        // 3.1) Process each zone with their refined search estimates
        // Do first the zones which take less time, as at some point
        // we may have to enforce the timeout.
        std::sort(zones.begin(), zones.end(), SearchParamLessThan());
        BOOST_FOREACH( SearchParam const& zone, zones ) {

          BBox2i left_region = zone.first + region_offset; // Kernel width offset
          left_region.min() -= half_kernel;
          left_region.max() += half_kernel;
          BBox2i right_region = left_region + zone.second.min();
          right_region.max() += zone.second.size();

          // Check timing estimate to see if we should go ahead with this zone or quit.
          double next_elapsed = m_seconds_per_op * search_volume(SearchParam(left_region, zone.second));
          if (m_corr_timeout > 0.0 && estim_elapsed + next_elapsed > m_corr_timeout){
            vw_out() << "Tile: " << bbox << " reached timeout: "
                     << m_corr_timeout << " s" << std::endl;
            break;
          }else
            estim_elapsed += next_elapsed;

          // See if it is time to actually accurately compute the time
          if (m_corr_timeout > 0.0 && estim_elapsed - prev_estim > measure_spacing){
            std::time (&end);
            double diff = std::difftime(end, start);
            estim_elapsed = diff;
            prev_estim = estim_elapsed;
          }

          // Compute left to right disparity vectors in this zone.
          crop(disparity, zone.first)
            = calc_disparity(m_cost_type,
                             crop(left_pyramid [level], left_region),
                             crop(right_pyramid[level], right_region),
                             left_region - left_region.min(),
                             zone.second.size(), m_kernel_size);

          // If at the last level and the user requested a left<->right consistency check,
          //   compute right to left disparity.
          if ( m_consistency_threshold >= 0 && level == 0 ) {

            // Check the time again before moving on with this
            double next_elapsed = m_seconds_per_op
              * search_volume(SearchParam(right_region, zone.second));
            if (m_corr_timeout > 0.0 && estim_elapsed + next_elapsed > m_corr_timeout){
              vw_out() << "Tile: " << bbox << " reached timeout: "
                       << m_corr_timeout << " s" << std::endl;
              break;
            }else{
              estim_elapsed += next_elapsed;
            }
            // Compute right to left disparity in this zone
            ImageView<pixel_type> rl_result
              = calc_disparity(m_cost_type,
                               crop(edge_extend(right_pyramid[level]), right_region),
                               crop(edge_extend(left_pyramid [level]),
                                    left_region - zone.second.size()),
                               right_region - right_region.min(),
                               zone.second.size(), m_kernel_size)
              - pixel_type(zone.second.size());

            // Find pixels where the disparity distance is greater than m_consistency_threshold
            stereo::cross_corr_consistency_check(crop(disparity,zone.first),
                                                  rl_result,
                                                 m_consistency_threshold, false);
          } // End of last level right to left disparity check

          // Fix the offsets to account for cropping.
          crop(disparity, zone.first) += pixel_type(zone.second.min());
        } // End of zone loop

        // 3.2a) Filter the disparity so we are not processing more than we need to.
        //       - Inner function filtering is only to catch "speckle" type noise of individual ouliers.
        //       - Outer function just merges the masks over the filtered disparity image.
        const int32 rm_half_kernel = 5;
        const float rm_min_matches_percent = 0.5;
        const float rm_threshold = 3.0;
        if ( level != 0 ) {
          disparity = disparity_mask(disparity_cleanup_using_thresh
                                       (disparity,
                                        rm_half_kernel, rm_half_kernel,
                                        rm_threshold,
                                        rm_min_matches_percent),
                                       left_mask_pyramid[level],
                                       right_mask_pyramid[level]);
        } else {
          // We don't do a single hot pixel check on the final level as it leaves a border.
          disparity = disparity_mask(rm_outliers_using_thresh
                                       (disparity,
                                        rm_half_kernel, rm_half_kernel,
                                        rm_threshold,
                                        rm_min_matches_percent),
                                       left_mask_pyramid[level],
                                       right_mask_pyramid[level]);
        }

        // 3.2b) Refine search estimates but never let them go beyond
        // the search region defined by the user
        if ( level != 0 ) {
          zones.clear();

          subdivide_regions( disparity, bounding_box(disparity),
                             zones, m_kernel_size );

          if (0) {
            BBox2i scaled = bbox/2;
            std::ostringstream ostr;
            ostr << "disparity_" << scaled.min()[0] << "_"
                 << scaled.min()[1] << "_" << scaled.max()[0] << "_"
                 << scaled.max()[1] << "_" << level;
            write_image( ostr.str() + ".tif", pixel_cast<PixelMask<Vector2f> >(disparity) );
            std::ofstream f( (ostr.str() + "_zone.txt").c_str() );
            BOOST_FOREACH( SearchParam& zone, zones ) {
              f << zone.first << " " << zone.second << "\n";
            }
            write_image( ostr.str() + "left.tif", normalize(left_pyramid[level]) );
            write_image( ostr.str() + "right.tif", normalize(right_pyramid[level]) );
            write_image( ostr.str() + "lmask.tif", left_mask_pyramid[level] );
            write_image( ostr.str() + "rmask.tif", right_mask_pyramid[level] );
            f.close();
          }
          scaling >>= 1;
          // Scale search range defines the maximum search range that
          // is possible in the next step. This (at lower levels) will
          // actually be larger than the search range that the user
          // specified. We are able to due this because we are taking
          // advantage of the half kernel padding needed at the hight
          // level of the pyramid.
          BBox2i scale_search_region(0,0,
                                     right_pyramid[level-1].cols() - left_pyramid[level-1].cols(),
                                     right_pyramid[level-1].rows() - left_pyramid[level-1].rows() );
          BBox2i next_zone_size = bounding_box( left_mask_pyramid[level-1] );
          BOOST_FOREACH( SearchParam& zone, zones ) {
            zone.first *= 2;
            zone.first.crop( next_zone_size );
            zone.second *= 2;
            zone.second.expand(2); // This is practically required. Our
            // correlation will fail if the search has only one
            // solution.
            zone.second.crop( scale_search_region );
          }
        }
      } // End of the level loop

      VW_ASSERT( bbox.size() == bounding_box(disparity).size(),
                 MathErr() << "PyramidCorrelation: Solved disparity doesn't match requested bbox size." );

#if VW_DEBUG_LEVEL > 0
      watch.stop();
      double elapsed = watch.elapsed_seconds();
      vw_out(DebugMessage,"stereo") << "Tile " << bbox << " processed in "
                                    << elapsed << " s\n";
      if (m_corr_timeout > 0.0){
        vw_out(DebugMessage,"stereo")
          << "Elapsed (actual/estimated/ratio): " << elapsed << ' '
          << estim_elapsed << ' ' << elapsed/estim_elapsed << std::endl;
      }
#endif

      // 5.0) Reposition our result back into the global
      // solution. Also we need to correct for the offset we applied
      // to the search region.
      return prerasterize_type(disparity + pixel_type(m_search_region.min()),
                               -bbox.min().x(), -bbox.min().y(),
                               cols(), rows() );
    } // End function prerasterize

    template <class DestT>
    inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
      vw::rasterize(prerasterize(bbox), dest, bbox);
    }
  }; // End class PyramidCorrelationView

  template <class Image1T, class Image2T, class Mask1T, class Mask2T, class PreFilterT>
  PyramidCorrelationView<Image1T,Image2T,Mask1T,Mask2T,PreFilterT>
  pyramid_correlate( ImageViewBase<Image1T> const& left,
                     ImageViewBase<Image2T> const& right,
                     ImageViewBase<Mask1T> const& left_mask,
                     ImageViewBase<Mask2T> const& right_mask,
                     PreFilterBase<PreFilterT> const& filter,
                     BBox2i const& search_region, Vector2i const& kernel_size,
                     CostFunctionType cost_type,
                     int corr_timeout, double seconds_per_op,
                     float consistency_threshold,
                     int32 max_pyramid_levels) {
    typedef PyramidCorrelationView<Image1T,Image2T,Mask1T,Mask2T,PreFilterT> result_type;
    return result_type( left.impl(), right.impl(), left_mask.impl(),
                        right_mask.impl(), filter.impl(), search_region,
                        kernel_size, cost_type,
                        corr_timeout, seconds_per_op,
                        consistency_threshold, max_pyramid_levels );
  }

}} // namespace vw::stereo

#endif//__VW_STEREO_CORRELATION_VIEW_H__
