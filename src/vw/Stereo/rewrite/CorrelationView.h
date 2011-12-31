// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __VW_STEREO_REWRITE_CORRELATION_VIEW_H__
#define __VW_STEREO_REWRITE_CORRELATION_VIEW_H__

#include <vw/Core/Exception.h>
#include <vw/Core/Stopwatch.h>
#include <vw/Image/Algorithms.h>
#include <vw/FileIO.h>
#include <vw/Stereo/rewrite/Correlation.h>
#include <boost/foreach.hpp>

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
                     PreFilterBase<PreFilterT> const& prefilter,
                     BBox2i const& search_region, Vector2i const& kernel_size,
                     CostFunctionType cost_type = ABSOLUTE_DIFFERENCE,
                     float consistency_threshold = -1 ) :
      m_left_image(left.impl()), m_right_image(right.impl()),
      m_prefilter(prefilter.impl()), m_search_region(search_region), m_kernel_size(kernel_size),
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

#if VW_DEBUG_LEVEL > 0
      Stopwatch watch;
      watch.start();
#endif

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
          best_of_search_convolution<NCCCost>( crop(m_prefilter.filter(m_left_image),left_region),
                                               crop(m_prefilter.filter(m_right_image),right_region),
                                               left_region - left_region.min(),
                                               m_search_region.size() + Vector2i(1,1),
                                               m_kernel_size );
        break;
      case SQUARED_DIFFERENCE:
        result =
          best_of_search_convolution<SquaredCost>( crop(m_prefilter.filter(m_left_image),left_region),
                                                   crop(m_prefilter.filter(m_right_image),right_region),
                                                   left_region - left_region.min(),
                                                   m_search_region.size() + Vector2i(1,1),
                                                   m_kernel_size );
        break;
      case ABSOLUTE_DIFFERENCE:
      default:
        result =
          best_of_search_convolution<AbsoluteCost>( crop(m_prefilter.filter(m_left_image),left_region),
                                                    crop(m_prefilter.filter(m_right_image),right_region),
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
            best_of_search_convolution<NCCCost>( crop(m_prefilter.filter(m_right_image),right_region),
                                                 crop(m_prefilter.filter(m_left_image),left_region-(m_search_region.size()+Vector2i(1,1))),
                                                 right_region - right_region.min(),
                                                 m_search_region.size() + Vector2i(1,1),
                                                 m_kernel_size ) -
            pixel_type(m_search_region.size()+Vector2i(1,1));
          break;
        case SQUARED_DIFFERENCE:
          // Getting the crops correctly here is not important as best
          // of search convolution will recrop. The important bit is
          // just aligning up the origins.
          rl_result =
            best_of_search_convolution<SquaredCost>( crop(m_prefilter.filter(m_right_image),right_region),
                                                     crop(m_prefilter.filter(m_left_image),left_region-(m_search_region.size()+Vector2i(1,1))),
                                                     right_region - right_region.min(),
                                                     m_search_region.size() + Vector2i(1,1),
                                                     m_kernel_size ) -
            pixel_type(m_search_region.size()+Vector2i(1,1));
          break;
        case ABSOLUTE_DIFFERENCE:
        default:
          // Getting the crops correctly here is not important as best
          // of search convolution will recrop. The important bit is
          // just aligning up the origins.
          rl_result =
            best_of_search_convolution<AbsoluteCost>( crop(m_prefilter.filter(m_right_image),right_region),
                                                      crop(m_prefilter.filter(m_left_image),left_region-(m_search_region.size()+Vector2i(1,1))),
                                                      right_region - right_region.min(),
                                                      m_search_region.size() + Vector2i(1,1),
                                                      m_kernel_size ) -
            pixel_type(m_search_region.size()+Vector2i(1,1));
        }
        stereo::cross_corr_consistency_check( result, rl_result,
                                              m_consistency_threshold, false );
      }
      VW_DEBUG_ASSERT( bbox.size() == bounding_box(result).size(),
                       MathErr() << "CorrelationView::prerasterize got a bad return from best_of_search_convolution." );

      // 5.) Convert back to original coordinates
      result += pixel_type(m_search_region.min());

#if VW_DEBUG_LEVEL > 0
      watch.stop();
      vw_out(DebugMessage,"stereo") << "Tile " << bbox << " processed in " << watch.elapsed_seconds() << "\n";
#endif

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
             PreFilterBase<PreFilterT> const& filter,
             BBox2i const& search_region, Vector2i const& kernel_size,
             CostFunctionType cost_type = ABSOLUTE_DIFFERENCE,
             float consistency_threshold = -1 ) {
    typedef CorrelationView<Image1T,Image2T,PreFilterT> result_type;
    return result_type( left.impl(), right.impl(), filter.impl(), search_region,
                        kernel_size, cost_type, consistency_threshold );
  }

  /// An image view for performing pyramid image correlation (Faster
  /// than CorrelationView)
  template <class Image1T, class Image2T, class PreFilterT>
  class PyramidCorrelationView : public ImageViewBase<PyramidCorrelationView<Image1T,Image2T,PreFilterT> > {

    Image1T m_left_image;
    Image2T m_right_image;
    PreFilterT m_prefilter;
    BBox2i m_search_region;
    Vector2i m_kernel_size;
    rewrite::CostFunctionType m_cost_type;
    float m_consistency_threshold; // < 0 = means don't do a consistency check
    int32 m_max_level_by_search;

  public:
    typedef PixelMask<Vector2i> pixel_type;
    typedef PixelMask<Vector2i> result_type;
    typedef ProceduralPixelAccessor<PyramidCorrelationView> pixel_accessor;

    PyramidCorrelationView( ImageViewBase<Image1T> const& left,
                            ImageViewBase<Image2T> const& right,
                            PreFilterBase<PreFilterT> const& prefilter,
                            BBox2i const& search_region, Vector2i const& kernel_size,
                            CostFunctionType cost_type = ABSOLUTE_DIFFERENCE,
                            float consistency_threshold = -1 ) :
      m_left_image(left.impl()), m_right_image(right.impl()),
      m_prefilter(prefilter.impl()), m_search_region(search_region), m_kernel_size(kernel_size),
      m_cost_type(cost_type), m_consistency_threshold(consistency_threshold) {
      // Calculating max pyramid levels according to the supplied
      // search region.
      int32 largest_search = max( search_region.size() );
      m_max_level_by_search = std::floor(std::log(float(largest_search))/std::log(2.0f)) - 1;
      if ( m_max_level_by_search < 0 )
        m_max_level_by_search = 0;
    }

    // Standard required ImageView interfaces
    inline int32 cols() const { return m_left_image.cols(); }
    inline int32 rows() const { return m_left_image.rows(); }
    inline int32 planes() const { return 1; }

    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }
    inline pixel_type operator()( int32 /*i*/, int32 /*j*/, int32 /*p*/ = 0) const {
      vw_throw( NoImplErr() << "PyramidCorrelationView::operator()(....) has not been implemented." );
      return pixel_type();
    }

    // Block rasterization section that does actual work
    typedef CropView<ImageView<pixel_type> > prerasterize_type;
    inline prerasterize_type prerasterize(BBox2i const& bbox) const {

#if VW_DEBUG_LEVEL > 0
      Stopwatch watch;
      watch.start();
#endif

      // 1.0) Determining the number of levels to process
      //      There's a maximum base on kernel size. There's also
      //      maximum defined by the search range. Here we determine
      //      the maximum based on kernel size and current bbox.
      int32 smallest_bbox = math::min(bbox.size());
      int32 largest_kernel = math::max(m_kernel_size);
      int32 max_pyramid_levels = std::floor(log(smallest_bbox)/log(2.0f) - log(largest_kernel)/log(2.0f));
      if ( m_max_level_by_search < max_pyramid_levels )
        max_pyramid_levels = m_max_level_by_search;
      if ( max_pyramid_levels < 1 )
        max_pyramid_levels = 0;
      Vector2i half_kernel = m_kernel_size/2;

      // 2.0) Build the pyramid
      std::vector<ImageView<typename Image1T::pixel_type> > left_pyramid(max_pyramid_levels + 1 );
      std::vector<ImageView<typename Image2T::pixel_type> > righ_pyramid(max_pyramid_levels + 1 );
      int32 max_upscaling = 1 << max_pyramid_levels;
      BBox2i left_global_region, right_global_region;
      {
        // Here we apply the prefilter first. However it might be that
        // we should build the gaussian pyramid first and then apply
        // the filter to each region separately.
        left_global_region = bbox;
        left_global_region.min() -= half_kernel * max_upscaling;
        left_global_region.max() += half_kernel * max_upscaling;
        right_global_region = left_global_region + m_search_region.min();
        right_global_region.max() += m_search_region.size() + Vector2i(max_upscaling,max_upscaling);
        left_pyramid[0] = crop(edge_extend(m_left_image),left_global_region);
        righ_pyramid[0] = crop(edge_extend(m_right_image),right_global_region);

#if VW_DEBUG_LEVEL > 0
        vw_out(DebugMessage,"stereo") << " > Left ROI: " << left_global_region
                                      << "\n > Right ROI: " << right_global_region << "\n";
#endif

        // Szeliski's book recommended this simple kernel. This
        // operation is quickly becoming a time sink, we might
        // possibly want to write an integer optimized version.
        std::vector<typename DefaultKernelT<typename Image1T::pixel_type>::type > kernel(5);
        kernel[0] = kernel[4] = 1.0/16.0;
        kernel[1] = kernel[3] = 4.0/16.0;
        kernel[2] = 6.0/16.0;

        for ( int32 i = 0; i < max_pyramid_levels; ++i ) {
          left_pyramid[i+1] = subsample(separable_convolution_filter(left_pyramid[i],kernel,kernel),2);
          righ_pyramid[i+1] = subsample(separable_convolution_filter(righ_pyramid[i],kernel,kernel),2);
          left_pyramid[i] = m_prefilter.filter(left_pyramid[i]);
          righ_pyramid[i] = m_prefilter.filter(righ_pyramid[i]);
        }
        left_pyramid[max_pyramid_levels] = m_prefilter.filter(left_pyramid[max_pyramid_levels]);
        righ_pyramid[max_pyramid_levels] = m_prefilter.filter(righ_pyramid[max_pyramid_levels]);
      }

      // 3.0) Actually perform correlation now
      ImageView<pixel_type > disparity;
      std::list<SearchParam> zones;
      zones.push_back( SearchParam( BBox2i(0,0,bbox.width()/max_upscaling,
                                             bbox.height()/max_upscaling),
                                    BBox2i(0,0,m_search_region.width()/max_upscaling+1,
                                           m_search_region.height()/max_upscaling+1)) );
      for ( int32 level = max_pyramid_levels; level >= 0; --level) {
        int32 scaling = 1 << level;
        disparity.set_size( bbox.width()/scaling, bbox.height()/scaling );

        // 3.1) Process each zone with their refined search estimates
        BOOST_FOREACH( SearchParam const& zone, zones ) {
          BBox2i left_region = zone.first + (max_upscaling/scaling)*half_kernel;
          left_region.min() -= half_kernel;
          left_region.max() += half_kernel;
          BBox2i righ_region = left_region + zone.second.min();
          righ_region.max() += zone.second.size();

          switch ( m_cost_type ) {
          case CROSS_CORRELATION:
            crop(disparity,zone.first) =
              best_of_search_convolution<NCCCost>(
                crop(left_pyramid[level],left_region),
                crop(righ_pyramid[level],righ_region),
                left_region - left_region.min(),
                zone.second.size(), m_kernel_size );
            break;
          case SQUARED_DIFFERENCE:
            crop(disparity,zone.first) =
              best_of_search_convolution<SquaredCost>(
                crop(left_pyramid[level],left_region),
                crop(righ_pyramid[level],righ_region),
                left_region - left_region.min(),
                zone.second.size(), m_kernel_size );
            break;
          case ABSOLUTE_DIFFERENCE:
          default:
            crop(disparity,zone.first) =
              best_of_search_convolution<AbsoluteCost>(
                crop(left_pyramid[level],left_region),
                crop(righ_pyramid[level],righ_region),
                left_region - left_region.min(),
                zone.second.size(), m_kernel_size );
          }

          if ( m_consistency_threshold >= 0 && level == 0 ) {
            ImageView<pixel_type> rl_result;

            switch ( m_cost_type ) {
            case CROSS_CORRELATION:
              rl_result =
                best_of_search_convolution<NCCCost>(
                  crop(edge_extend(righ_pyramid[level]),righ_region),
                  crop(edge_extend(left_pyramid[level]),left_region - zone.second.size()),
                  righ_region - righ_region.min(),
                  zone.second.size(), m_kernel_size ) - pixel_type(zone.second.size());
              break;
            case SQUARED_DIFFERENCE:
              rl_result =
                best_of_search_convolution<SquaredCost>(
                  crop(edge_extend(righ_pyramid[level]),righ_region),
                  crop(edge_extend(left_pyramid[level]),left_region - zone.second.size()),
                  righ_region - righ_region.min(),
                  zone.second.size(), m_kernel_size ) - pixel_type(zone.second.size());
              break;
            case ABSOLUTE_DIFFERENCE:
            default:
              rl_result =
               best_of_search_convolution<AbsoluteCost>(
                  crop(edge_extend(righ_pyramid[level]),righ_region),
                  crop(edge_extend(left_pyramid[level]),left_region - zone.second.size()),
                  righ_region - righ_region.min(),
                  zone.second.size(), m_kernel_size ) - pixel_type(zone.second.size());
            }

            stereo::cross_corr_consistency_check( crop(disparity,zone.first),
                                                  rl_result, m_consistency_threshold, false );
          }

          // Fix the offset
          crop(disparity,zone.first) += pixel_type(zone.second.min());

        } // end of zone loop

        // 3.2) Refine search estimates but never let them go beyond
        // the search region defined by the user
        if ( level != 0 ) {
          zones.clear();

          // 3.2a) Filter the disparity so we are not processing more
          // that we need to.
          const int32 rm_half_kernel = 5;
          const float rm_min_matches_percent = 0.5;
          const float rm_threshold = 3.0;
          disparity = disparity_clean_up(disparity,
                                         rm_half_kernel, rm_half_kernel,
                                         rm_threshold,
                                         rm_min_matches_percent);

          subdivide_regions( disparity, bounding_box(disparity), zones, m_kernel_size );

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
            write_image( ostr.str() + "right.tif", normalize(righ_pyramid[level]) );
            f.close();
          }
          scaling >>= 1;
          BBox2i scale_search_region = (m_search_region - m_search_region.min())/scaling;
          scale_search_region.max() += Vector2i(1,1);
          BOOST_FOREACH( SearchParam& zone, zones ) {
            zone.first *= 2;
            zone.second *= 2;
            zone.second.expand(1); // This is practically required. Our
            // correlation will fail if the search has only one
            // solution. There's also a fudge factor here.
            zone.second.crop( scale_search_region );
          }
        }
      }

      VW_DEBUG_ASSERT( bbox.size() == bounding_box(disparity).size(),
                       MathErr() << "PyramidCorrelation: Solved disparity doesn't match requested bbox size." );

#if VW_DEBUG_LEVEL > 0
      watch.stop();
      vw_out(DebugMessage,"stereo") << "Tile " << bbox << " processed in " << watch.elapsed_seconds() << "\n";
#endif

      // 4.0) Reposition our result back into the global
      // solution. Also we need to correct for the offset we applied
      // to the search region.
      return prerasterize_type(disparity + pixel_type(m_search_region.min()),
                               -bbox.min().x(), -bbox.min().y(),
                               cols(), rows() );
    }

    template <class DestT>
    inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
      vw::rasterize(prerasterize(bbox), dest, bbox);
    }
  };

  template <class Image1T, class Image2T, class PreFilterT>
  PyramidCorrelationView<Image1T,Image2T,PreFilterT>
  pyramid_correlate( ImageViewBase<Image1T> const& left,
                     ImageViewBase<Image2T> const& right,
                     PreFilterBase<PreFilterT> const& filter,
                     BBox2i const& search_region, Vector2i const& kernel_size,
                     CostFunctionType cost_type = ABSOLUTE_DIFFERENCE,
                     float consistency_threshold = -1 ) {
    typedef PyramidCorrelationView<Image1T,Image2T,PreFilterT> result_type;
    return result_type( left.impl(), right.impl(), filter.impl(), search_region,
                        kernel_size, cost_type, consistency_threshold );
  }

  /// An image view for performing pyramid image correlation (Faster
  /// than CorrelationView). This version accepts a mask and is good
  /// for images that have been map projected and have odd sections
  /// that are invalid.
  template <class Image1T, class Image2T, class Mask1T, class Mask2T, class PreFilterT>
  class PyramidMaskCorrelationView : public ImageViewBase<PyramidMaskCorrelationView<Image1T,Image2T, Mask1T, Mask2T,PreFilterT> > {

    Image1T m_left_image;
    Image2T m_right_image;
    Mask1T  m_left_mask;
    Mask2T  m_right_mask;
    PreFilterT m_prefilter;
    BBox2i m_search_region;
    Vector2i m_kernel_size;
    rewrite::CostFunctionType m_cost_type;
    float m_consistency_threshold; // < 0 = means don't do a consistency check
    int32 m_max_level_by_search;

  public:
    typedef PixelMask<Vector2i> pixel_type;
    typedef PixelMask<Vector2i> result_type;
    typedef ProceduralPixelAccessor<PyramidMaskCorrelationView> pixel_accessor;

    PyramidMaskCorrelationView( ImageViewBase<Image1T> const& left,
                                ImageViewBase<Image2T> const& right,
                                ImageViewBase<Mask1T> const& left_mask,
                                ImageViewBase<Mask2T> const& right_mask,
                                PreFilterBase<PreFilterT> const& prefilter,
                                BBox2i const& search_region, Vector2i const& kernel_size,
                                CostFunctionType cost_type = ABSOLUTE_DIFFERENCE,
                                float consistency_threshold = -1 ) :
      m_left_image(left.impl()), m_right_image(right.impl()),
      m_left_mask(left_mask.impl()), m_right_mask(right_mask.impl()),
      m_prefilter(prefilter.impl()), m_search_region(search_region), m_kernel_size(kernel_size),
      m_cost_type(cost_type), m_consistency_threshold(consistency_threshold) {
      // Calculating max pyramid levels according to the supplied
      // search region.
      int32 largest_search = max( search_region.size() );
      m_max_level_by_search = std::floor(std::log(float(largest_search))/std::log(2.0f)) - 1;
      if ( m_max_level_by_search < 0 )
        m_max_level_by_search = 0;
    }

    // Standard required ImageView interfaces
    inline int32 cols() const { return m_left_image.cols(); }
    inline int32 rows() const { return m_left_image.rows(); }
    inline int32 planes() const { return 1; }

    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }
    inline pixel_type operator()( int32 /*i*/, int32 /*j*/, int32 /*p*/ = 0) const {
      vw_throw( NoImplErr() << "PyramidMaskCorrelationView::operator()(....) has not been implemented." );
      return pixel_type();
    }

    // Block rasterization section that does actual work
    typedef CropView<ImageView<pixel_type> > prerasterize_type;
    inline prerasterize_type prerasterize(BBox2i const& bbox) const {

#if VW_DEBUG_LEVEL > 0
      Stopwatch watch;
      watch.start();
#endif

      // 1.0) Determining the number of levels to process
      //      There's a maximum base on kernel size. There's also
      //      maximum defined by the search range. Here we determine
      //      the maximum based on kernel size and current bbox.
      int32 smallest_bbox = math::min(bbox.size());
      int32 largest_kernel = math::max(m_kernel_size);
      int32 max_pyramid_levels = std::floor(log(smallest_bbox)/log(2.0f) - log(largest_kernel)/log(2.0f));
      if ( m_max_level_by_search < max_pyramid_levels )
        max_pyramid_levels = m_max_level_by_search;
      if ( max_pyramid_levels < 1 )
        max_pyramid_levels = 0;
      Vector2i half_kernel = m_kernel_size/2;

      // 2.0) Build the pyramid
      std::vector<ImageView<typename Image1T::pixel_type> > left_pyramid(max_pyramid_levels + 1 );
      std::vector<ImageView<typename Image2T::pixel_type> > righ_pyramid(max_pyramid_levels + 1 );
      std::vector<ImageView<typename Mask1T::pixel_type> > left_mask_pyramid(max_pyramid_levels + 1 );
      std::vector<ImageView<typename Mask2T::pixel_type> > righ_mask_pyramid(max_pyramid_levels + 1 );
      int32 max_upscaling = 1 << max_pyramid_levels;
      BBox2i left_global_region, right_global_region;
      {
        // Here we apply the prefilter first. However it might be that
        // we should build the gaussian pyramid first and then apply
        // the filter to each region separately.
        left_global_region = bbox;
        left_global_region.min() -= half_kernel * max_upscaling;
        left_global_region.max() += half_kernel * max_upscaling;
        right_global_region = left_global_region + m_search_region.min();
        right_global_region.max() += m_search_region.size() + Vector2i(max_upscaling,max_upscaling);
        left_pyramid[0] = crop(edge_extend(m_left_image),left_global_region);
        righ_pyramid[0] = crop(edge_extend(m_right_image),right_global_region);
        BBox2i right_mask = bbox + m_search_region.min();
        right_mask.max() += m_search_region.size();
        left_mask_pyramid[0] = crop(edge_extend(m_left_mask),bbox);
        righ_mask_pyramid[0] = crop(edge_extend(m_right_mask),right_mask);

#if VW_DEBUG_LEVEL > 0
        vw_out(DebugMessage,"stereo") << " > Left ROI: " << left_global_region
                                      << "\n > Right ROI: " << right_global_region << "\n";
#endif

        // Szeliski's book recommended this simple kernel. This
        // operation is quickly becoming a time sink, we might
        // possibly want to write an integer optimized version.
        std::vector<typename DefaultKernelT<typename Image1T::pixel_type>::type > kernel(5);
        kernel[0] = kernel[4] = 1.0/16.0;
        kernel[1] = kernel[3] = 4.0/16.0;
        kernel[2] = 6.0/16.0;
        std::vector<uint8> mask_kern(max(m_kernel_size));
        std::fill(mask_kern.begin(), mask_kern.end(), 1 );

        for ( int32 i = 0; i < max_pyramid_levels; ++i ) {
          left_pyramid[i+1] = subsample(separable_convolution_filter(left_pyramid[i],kernel,kernel),2);
          righ_pyramid[i+1] = subsample(separable_convolution_filter(righ_pyramid[i],kernel,kernel),2);
          left_pyramid[i] = m_prefilter.filter(left_pyramid[i]);
          righ_pyramid[i] = m_prefilter.filter(righ_pyramid[i]);
          left_mask_pyramid[i+1] = threshold(separable_convolution_filter(subsample(threshold(left_mask_pyramid[i],0,1),2),mask_kern,mask_kern), 0, 255);
          righ_mask_pyramid[i+1] = threshold(separable_convolution_filter(subsample(threshold(righ_mask_pyramid[i],0,1),2),mask_kern,mask_kern), 0, 255);
        }
        left_pyramid[max_pyramid_levels] = m_prefilter.filter(left_pyramid[max_pyramid_levels]);
        righ_pyramid[max_pyramid_levels] = m_prefilter.filter(righ_pyramid[max_pyramid_levels]);
      }

      // 3.0) Actually perform correlation now
      ImageView<pixel_type > disparity;
      std::list<SearchParam> zones;
      zones.push_back( SearchParam( bounding_box(left_mask_pyramid[max_pyramid_levels]),
                                    BBox2i(0,0,m_search_region.width()/max_upscaling+1,
                                           m_search_region.height()/max_upscaling+1)) );
      for ( int32 level = max_pyramid_levels; level >= 0; --level) {
        int32 scaling = 1 << level;
        disparity.set_size( left_mask_pyramid[level] );
        Vector2i region_offset = max_upscaling*half_kernel/scaling;

        // 3.1) Process each zone with their refined search estimates
        BOOST_FOREACH( SearchParam const& zone, zones ) {
          BBox2i left_region = zone.first + region_offset;
          left_region.min() -= half_kernel;
          left_region.max() += half_kernel;
          BBox2i righ_region = left_region + zone.second.min();
          righ_region.max() += zone.second.size();

          switch ( m_cost_type ) {
          case CROSS_CORRELATION:
            crop(disparity,zone.first) =
              best_of_search_convolution<NCCCost>(
                crop(left_pyramid[level],left_region),
                crop(righ_pyramid[level],righ_region),
                left_region - left_region.min(),
                zone.second.size(), m_kernel_size );
            break;
          case SQUARED_DIFFERENCE:
            crop(disparity,zone.first) =
              best_of_search_convolution<SquaredCost>(
                crop(left_pyramid[level],left_region),
                crop(righ_pyramid[level],righ_region),
                left_region - left_region.min(),
                zone.second.size(), m_kernel_size );
            break;
          case ABSOLUTE_DIFFERENCE:
          default:
            crop(disparity,zone.first) =
              best_of_search_convolution<AbsoluteCost>(
                crop(left_pyramid[level],left_region),
                crop(righ_pyramid[level],righ_region),
                left_region - left_region.min(),
                zone.second.size(), m_kernel_size );
          }

          if ( m_consistency_threshold >= 0 && level == 0 ) {
            ImageView<pixel_type> rl_result;

            switch ( m_cost_type ) {
            case CROSS_CORRELATION:
              rl_result =
                best_of_search_convolution<NCCCost>(
                  crop(edge_extend(righ_pyramid[level]),righ_region),
                  crop(edge_extend(left_pyramid[level]),left_region - zone.second.size()),
                  righ_region - righ_region.min(),
                  zone.second.size(), m_kernel_size ) - pixel_type(zone.second.size());
              break;
            case SQUARED_DIFFERENCE:
              rl_result =
                best_of_search_convolution<SquaredCost>(
                  crop(edge_extend(righ_pyramid[level]),righ_region),
                  crop(edge_extend(left_pyramid[level]),left_region - zone.second.size()),
                  righ_region - righ_region.min(),
                  zone.second.size(), m_kernel_size ) - pixel_type(zone.second.size());
              break;
            case ABSOLUTE_DIFFERENCE:
            default:
              rl_result =
               best_of_search_convolution<AbsoluteCost>(
                  crop(edge_extend(righ_pyramid[level]),righ_region),
                  crop(edge_extend(left_pyramid[level]),left_region - zone.second.size()),
                  righ_region - righ_region.min(),
                  zone.second.size(), m_kernel_size ) - pixel_type(zone.second.size());
            }

            stereo::cross_corr_consistency_check( crop(disparity,zone.first),
                                                  rl_result, m_consistency_threshold, false );
          }

          // Fix the offset
          crop(disparity,zone.first) += pixel_type(zone.second.min());
        } // end of zone loop

        // 3.2) Refine search estimates but never let them go beyond
        // the search region defined by the user
        if ( level != 0 ) {
          zones.clear();

          // 3.2a) Filter the disparity so we are not processing more
          // that we need to.
          const int32 rm_half_kernel = 5;
          const float rm_min_matches_percent = 0.5;
          const float rm_threshold = 3.0;
          disparity = disparity_mask(disparity_clean_up(disparity,
                                                        rm_half_kernel, rm_half_kernel,
                                                        rm_threshold,
                                                        rm_min_matches_percent),
                                     left_mask_pyramid[level],
                                     righ_mask_pyramid[level]);

          subdivide_regions( disparity, bounding_box(disparity), zones, m_kernel_size );

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
            write_image( ostr.str() + "right.tif", normalize(righ_pyramid[level]) );
            write_image( ostr.str() + "lmask.tif", left_mask_pyramid[level] );
            write_image( ostr.str() + "rmask.tif", righ_mask_pyramid[level] );
            f.close();
          }
          scaling >>= 1;
          BBox2i scale_search_region = (m_search_region - m_search_region.min())/scaling;
          scale_search_region.max() += Vector2i(1,1);
          BBox2i next_zone_size = bounding_box( left_mask_pyramid[level-1] );
          BOOST_FOREACH( SearchParam& zone, zones ) {
            zone.first *= 2;
            zone.first.crop( next_zone_size );
            zone.second *= 2;
            zone.second.expand(1); // This is practically required. Our
            // correlation will fail if the search has only one
            // solution. There's also a fudge factor here.
            zone.second.crop( scale_search_region );
          }
        }
      }

      VW_DEBUG_ASSERT( bbox.size() == bounding_box(disparity).size(),
                       MathErr() << "PyramidCorrelation: Solved disparity doesn't match requested bbox size." );

#if VW_DEBUG_LEVEL > 0
      watch.stop();
      vw_out(DebugMessage,"stereo") << "Tile " << bbox << " processed in " << watch.elapsed_seconds() << "\n";
#endif

      // 4.0) Reposition our result back into the global
      // solution. Also we need to correct for the offset we applied
      // to the search region.
      return prerasterize_type(disparity + pixel_type(m_search_region.min()),
                               -bbox.min().x(), -bbox.min().y(),
                               cols(), rows() );
    }

    template <class DestT>
    inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
      vw::rasterize(prerasterize(bbox), dest, bbox);
    }
  };

  template <class Image1T, class Image2T, class Mask1T, class Mask2T, class PreFilterT>
  PyramidMaskCorrelationView<Image1T,Image2T,Mask1T,Mask2T,PreFilterT>
  pyramid_correlate( ImageViewBase<Image1T> const& left,
                     ImageViewBase<Image2T> const& right,
                     ImageViewBase<Mask1T> const& left_mask,
                     ImageViewBase<Mask2T> const& right_mask,
                     PreFilterBase<PreFilterT> const& filter,
                     BBox2i const& search_region, Vector2i const& kernel_size,
                     CostFunctionType cost_type = ABSOLUTE_DIFFERENCE,
                     float consistency_threshold = -1 ) {
    typedef PyramidMaskCorrelationView<Image1T,Image2T,Mask1T,Mask2T,PreFilterT> result_type;
    return result_type( left.impl(), right.impl(), left_mask.impl(),
                        right_mask.impl(), filter.impl(), search_region,
                        kernel_size, cost_type, consistency_threshold );
  }

}}} // namespace vw::stereo::rewrite

#endif//__VW_STEREO_REWRITE_CORRELATION_VIEW_H__
