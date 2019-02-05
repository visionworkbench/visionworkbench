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
#include <vw/Image/ErodeView.h>
#include <vw/Image/PerPixelAccessorViews.h>
#include <vw/Image/ImageIO.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Stereo/Correlation.h>
#include <vw/Stereo/Correlate.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Stereo/PreFilter.h>
#include <boost/foreach.hpp>
#include <ctime>

#include <vw/Stereo/SGM.h>

namespace vw {
namespace stereo {

  // TODO: Maybe delete this class, it is not used or maintained.
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
    inline prerasterize_type prerasterize(BBox2i const& bbox) const;

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

//=================================================================================


  enum CorrelationAlgorithm {
    CORRELATION_WINDOW    = 0, ///< Use a local sliding search window.
    CORRELATION_SGM       = 1, ///< Use a Semi-Global Matching algorithm.
    CORRELATION_MGM       = 2, ///< Use the More-Global Matching algorithm.
    CORRELATION_FINAL_MGM = 3  ///< SGM until resolution level 0, then MGM.
  };

  /// An image view for performing pyramid image correlation (Faster than CorrelationView).
  template <class Image1T, class Image2T, class Mask1T, class Mask2T>
  class PyramidCorrelationView : 
      public ImageViewBase<PyramidCorrelationView<Image1T,Image2T, Mask1T, Mask2T> > {

  public: // Definitions
  typedef PixelMask<Vector2i> pixel_typeI;
  // TODO: Can we delete result_type from all image views???
    typedef PixelMask<Vector2f> pixel_type; // TODO: Safe to ever have differ?
    typedef PixelMask<Vector2f> result_type;
    typedef ProceduralPixelAccessor<PyramidCorrelationView> pixel_accessor;

  public: // Functions

    /// Initialize the view
    /// - Set blob_filter_area > 0 to filter out disparity blobs.
    PyramidCorrelationView( ImageViewBase<Image1T> const& left,
                            ImageViewBase<Image2T> const& right,
                            ImageViewBase<Mask1T > const& left_mask,
                            ImageViewBase<Mask2T > const& right_mask,
                            PrefilterModeType prefilter_mode, float prefilter_width,
                            BBox2i const& search_region, Vector2i const& kernel_size,
                            stereo::CostFunctionType cost_type,
                            int   corr_timeout, double seconds_per_op,
                            float consistency_threshold,
                            int min_consistency_level,
                            int filter_half_kernel,
                            int32 max_pyramid_levels,
                            CorrelationAlgorithm  algorithm = CORRELATION_WINDOW,
                            int   collar_size        = 0,
                            SemiGlobalMatcher::SgmSubpixelMode sgm_subpixel_mode = SemiGlobalMatcher::SUBPIXEL_LC_BLEND,
                            Vector2i  sgm_search_buffer = Vector2i(2,2),
                            size_t memory_limit_mb=6000,
                            int   blob_filter_area   = 0,
                            bool  write_debug_images = false) :
      m_left_image(left.impl()),     m_right_image(right.impl()),
      m_left_mask(left_mask.impl()), m_right_mask(right_mask.impl()),
      m_prefilter_mode(prefilter_mode), m_prefilter_width(prefilter_width),
      m_search_region(search_region), m_kernel_size(kernel_size),
      m_cost_type(cost_type),
      m_corr_timeout(corr_timeout), m_seconds_per_op(seconds_per_op),
      m_consistency_threshold(consistency_threshold),
      m_min_consistency_level(min_consistency_level),
      m_filter_half_kernel(filter_half_kernel), 
      m_blob_filter_area(blob_filter_area),
      m_algorithm(algorithm),
      m_collar_size(collar_size),
      m_sgm_subpixel_mode(sgm_subpixel_mode),
      m_sgm_search_buffer(sgm_search_buffer),
      m_memory_limit_mb(memory_limit_mb),
      m_write_debug_images(write_debug_images){

      // Quit if an invalid area was passed in
      double area = search_region.area();
      if (area != area) {
        vw_throw( ArgumentErr() << "PyramidCorrelationView: Invalid search region: "
                                << search_region.min() << ", " << search_region.max());
      }

      if (algorithm != CORRELATION_WINDOW)
        m_prefilter_mode = PREFILTER_NONE; // SGM/MGM works best with no prefilter
      
      // Calculating max pyramid levels according to the supplied search region.
      int32 largest_search = max( search_region.size() );
      m_max_level_by_search = std::floor(std::log(float(largest_search))/std::log(2.0f)) - 1;
      if ( m_max_level_by_search > max_pyramid_levels )
        m_max_level_by_search = max_pyramid_levels;
      if ( m_max_level_by_search < 0 )
        m_max_level_by_search = 0;
    } // End constructor

    // Standard required ImageView interfaces
    inline int32 cols  () const { return m_left_image.cols(); }
    inline int32 rows  () const { return m_left_image.rows(); }
    inline int32 planes() const { return 1; }

    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }
    inline result_type operator()( int32 /*i*/, int32 /*j*/, int32 /*p*/ = 0) const {
      vw_throw( NoImplErr() << "NewCorrelationView::operator()(....) has not been implemented." );
      return result_type();
    }

    /// Block rasterization section that does actual work
    typedef CropView<ImageView<result_type> > prerasterize_type;
    inline prerasterize_type prerasterize(BBox2i const& bbox) const;

    template <class DestT>
    inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
    
      vw_out(VerboseDebugMessage, "stereo") << "Rasterize called with box: " << bbox << std::endl;
    
      BBox2i proc_bbox = bbox;
      if (m_collar_size > 0)
        proc_bbox.expand(m_collar_size);
      vw_out(VerboseDebugMessage, "stereo") << "Collared raster box: " << proc_bbox << std::endl;
      vw::rasterize(prerasterize(proc_bbox), dest, bbox);
    }



  private: // Variables

    Image1T          m_left_image;
    Image2T          m_right_image;
    Mask1T           m_left_mask;
    Mask2T           m_right_mask;
    
    // These two variables pick a prefilter which is applied to each pyramid level
    PrefilterModeType m_prefilter_mode; ///< See Prefilter.h for the types
    float m_prefilter_width;     ///< Preprocessing filter width
    
    BBox2i           m_search_region;
    Vector2i         m_kernel_size;
    stereo::CostFunctionType m_cost_type;
    int              m_corr_timeout;
    // How long it takes to do one corr op with given kernel and cost function
    double m_seconds_per_op;
    float  m_consistency_threshold; // < 0 = means don't do a consistency check
    int    m_min_consistency_level; ///< Min level to do consistency check at
    int    m_filter_half_kernel; // Require half of pixels in this radius to be within 3 dist
    int32  m_max_level_by_search;
    
    /// Remove disparity blobs this size or less at the top level.
    /// - Set to zero to disable blob filtering.
    int    m_blob_filter_area;
    
    CorrelationAlgorithm m_algorithm; ///< Store the algorithm choice
    int m_collar_size;     ///< Expand the size of the image for each tile before correlating
    SemiGlobalMatcher::SgmSubpixelMode m_sgm_subpixel_mode; ///< Subpixel mode used by SGM algorithms
    Vector2i m_sgm_search_buffer;
    size_t m_memory_limit_mb;

    bool m_write_debug_images; ///< If true, write out a bunch of intermediate images.

  private: // Functions

    /// Downsample a mask by two.
    /// - If at least two mask pixels in a 2x2 region are on, the output pixel is on.
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

    /// Create the image pyramids needed by the prerasterize function.
    /// - Most of this function is spent figuring out the correct ROIs to use.
    bool build_image_pyramids(BBox2i const& bbox, int32 const max_pyramid_levels,
                              std::vector<ImageView<typename Image1T::pixel_type> > & left_pyramid,
                              std::vector<ImageView<typename Image2T::pixel_type> > & right_pyramid,
                              std::vector<ImageView<typename Mask1T::pixel_type > > & left_mask_pyramid,
                              std::vector<ImageView<typename Mask2T::pixel_type > > & right_mask_pyramid) const;

    /// Filter out isolated blobs of valid disparity regions which are usually wrong.
    /// - Using this can decrease run time in images with lots of little disparity islands.
    void disparity_blob_filter(ImageView<pixel_typeI > &disparity, int level,
                               int max_blob_area) const;



  }; // End class PyramidCorrelationView

  template <class Image1T, class Image2T, class Mask1T, class Mask2T>
  PyramidCorrelationView<Image1T,Image2T,Mask1T,Mask2T>
  pyramid_correlate( ImageViewBase<Image1T> const& left,
                     ImageViewBase<Image2T> const& right,
                     ImageViewBase<Mask1T > const& left_mask,
                     ImageViewBase<Mask2T > const& right_mask,
                     PrefilterModeType prefilter_mode, float prefilter_width,
                     BBox2i const& search_region, Vector2i const& kernel_size,
                     stereo::CostFunctionType cost_type,
                     int corr_timeout, double seconds_per_op,
                     float consistency_threshold,
                     int min_consistency_level,
                     int filter_half_kernel,
                     int32 max_pyramid_levels,
                     CorrelationAlgorithm  algorithm = CORRELATION_WINDOW,
                     int   collar_size        = 0,
                     SemiGlobalMatcher::SgmSubpixelMode sgm_subpixel_mode = SemiGlobalMatcher::SUBPIXEL_LC_BLEND,
                     Vector2i sgm_search_buffer=Vector2i(2,2),
                     size_t memory_limit_mb=6000,
                     int   blob_filter_area   = 0,
                     bool  write_debug_images =false) {
    typedef PyramidCorrelationView<Image1T,Image2T,Mask1T,Mask2T> result_type;
    return result_type( left.impl(),      right.impl(), 
                        left_mask.impl(), right_mask.impl(),
                        prefilter_mode, prefilter_width,
                        search_region,
                        kernel_size, cost_type,
                        corr_timeout, seconds_per_op,
                        consistency_threshold, min_consistency_level,
                        filter_half_kernel, max_pyramid_levels,
                        algorithm, collar_size, sgm_subpixel_mode, sgm_search_buffer, memory_limit_mb, blob_filter_area,
                        write_debug_images);
  }

}} // namespace vw::stereo

#include <vw/Stereo/CorrelationView.tcc>

#endif//__VW_STEREO_CORRELATION_VIEW_H__
