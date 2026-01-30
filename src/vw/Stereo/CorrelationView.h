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
#include <vw/Image/PerPixelAccessorViews.h>
#include <vw/Image/ImageIO.h>
#include <vw/Stereo/PreFilter.h>
#include <vw/Stereo/SGM.h>
#include <vw/Stereo/CorrelationAlgorithms.h>
#include <vw/Image/Manipulation.h>

namespace vw { namespace stereo {

  typedef ImageViewRef<PixelGray<float>> PixelGrayImageRef;
  typedef ImageViewRef<uint8> Int8ImageRef;
  
  /// An image view for performing pyramid image correlation
  class PyramidCorrelationView: public ImageViewBase<PyramidCorrelationView> {

  public: // Definitions
  typedef PixelMask<Vector2i> pixel_typeI;
    // TODO: Can we delete result_type from all image views???
    typedef PixelMask<Vector2f> pixel_type; // TODO: Safe to ever have differ?
    typedef PixelMask<Vector2f> result_type;
    typedef ProceduralPixelAccessor<PyramidCorrelationView> pixel_accessor;

  public: // Functions

    /// Initialize the view
    /// - Set blob_filter_area > 0 to filter out disparity blobs.
    PyramidCorrelationView(PixelGrayImageRef const& left,
                           PixelGrayImageRef const& right,
                           Int8ImageRef const& left_mask,
                           Int8ImageRef const& right_mask,
                           PrefilterModeType prefilter_mode, float prefilter_width,
                           BBox2i const& search_region, Vector2i const& kernel_size,
                           stereo::CostFunctionType cost_type,
                           int corr_timeout, double seconds_per_op,
                           float consistency_threshold,
                           int min_consistency_level,
                           int filter_half_kernel,
                           int32 max_pyramid_levels,
                           CorrelationAlgorithm  algorithm = VW_CORRELATION_BM,
                           int collar_size = 0,
                           SemiGlobalMatcher::SgmSubpixelMode sgm_subpixel_mode
                           = SemiGlobalMatcher::SUBPIXEL_LC_BLEND,
                           Vector2i sgm_search_buffer = Vector2i(2,2),
                           size_t memory_limit_mb=6000,
                           int blob_filter_area = 0,
                           ImageView<PixelMask<float>> * lr_disp_diff = NULL,
                           Vector2i const& region_ul = Vector2i(0, 0),
                           bool write_debug_images = false):
      m_left_image(left.impl()), m_right_image(right.impl()),
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
      m_lr_disp_diff(lr_disp_diff),
      m_region_ul(region_ul),
      m_write_debug_images(write_debug_images) {

      // Quit if an invalid area was passed in
      double area = search_region.area();
      if (area != area) {
        vw_throw(ArgumentErr() << "PyramidCorrelationView: Invalid search region: "
                               << search_region.min() << ", " << search_region.max());
      }

      if (algorithm != VW_CORRELATION_BM)
        m_prefilter_mode = PREFILTER_NONE; // SGM/MGM works best with no prefilter
      
      // Calculating max pyramid levels according to the supplied search region.
      int32 largest_search = max(search_region.size());
      m_max_level_by_search = std::floor(std::log(float(largest_search))/std::log(2.0f)) - 1;
      if (m_max_level_by_search > max_pyramid_levels)
        m_max_level_by_search = max_pyramid_levels;
      if (m_max_level_by_search < 0)
        m_max_level_by_search = 0;
    } // End constructor

    // Standard required ImageView interfaces
    inline int32 cols  () const { return m_left_image.cols(); }
    inline int32 rows  () const { return m_left_image.rows(); }
    inline int32 planes() const { return 1; }

    inline pixel_accessor origin() const { return pixel_accessor(*this, 0, 0); }
    inline result_type operator()(int32 /*i*/, int32 /*j*/, int32 /*p*/ = 0) const {
      vw_throw(NoImplErr() << "NewCorrelationView::operator() is not implemented.");
      return result_type();
    }

    /// Block rasterization section that does actual work
    typedef CropView<ImageView<result_type>> prerasterize_type;
    prerasterize_type prerasterize(BBox2i const& bbox) const;

    template <class DestT>
    inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
    
      vw_out(VerboseDebugMessage, "stereo") << "Rasterize called with box: " << bbox << "\n";
    
      BBox2i proc_bbox = bbox;
      if (m_collar_size > 0)
        proc_bbox.expand(m_collar_size);
      vw_out(VerboseDebugMessage, "stereo") << "Collared raster box: " << proc_bbox << "\n";
      vw::rasterize(prerasterize(proc_bbox), dest, bbox);
    }

  private: // Variables

    PixelGrayImageRef m_left_image;
    PixelGrayImageRef m_right_image;
    Int8ImageRef      m_left_mask;
    Int8ImageRef      m_right_mask;
    
    // These two variables pick a prefilter which is applied to each pyramid level
    PrefilterModeType m_prefilter_mode; ///< See Prefilter.h for the types
    float m_prefilter_width;            ///< Preprocessing filter width
    
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
    SemiGlobalMatcher::SgmSubpixelMode m_sgm_subpixel_mode; ///< Subpixel mode used by SGM
    Vector2i m_sgm_search_buffer;
    size_t m_memory_limit_mb;

    // Return here the L-R and R-L disparity discrepancy, if provided as a non-null pointer.
    ImageView<PixelMask<float>> * m_lr_disp_diff;

    // The upper-left corner of the region that will contain all
    // pixels to be processed in this class. This will help locating
    // where to populate the pixels of m_lr_disp_diff.
    Vector2i m_region_ul;
    
    bool m_write_debug_images; ///< If true, write out a bunch of intermediate images.

  private: // Functions

    /// Create the image pyramids needed by the prerasterize function.
    /// - Most of this function is spent figuring out the correct ROIs to use.
    bool build_image_pyramids
    (BBox2i const& bbox, int32 const max_pyramid_levels,
     std::vector<ImageView<PixelGray<float>>> & left_pyramid,
     std::vector<ImageView<PixelGray<float>>> & right_pyramid,
     std::vector<ImageView<uint8>> & left_mask_pyramid,
     std::vector<ImageView<uint8>> & right_mask_pyramid) const;
    
    /// Filter out isolated blobs of valid disparity regions which are usually wrong.
    /// - Using this can decrease run time in images with lots of little disparity islands.
    void disparity_blob_filter(ImageView<pixel_typeI > &disparity, int level,
                               int max_blob_area) const;

  }; // End class PyramidCorrelationView

  inline PyramidCorrelationView
  pyramid_correlate(PixelGrayImageRef const& left,
                    PixelGrayImageRef const& right,
                    Int8ImageRef const& left_mask,
                    Int8ImageRef const& right_mask,
                    PrefilterModeType prefilter_mode, float prefilter_width,
                    BBox2i const& search_region, Vector2i const& kernel_size,
                    stereo::CostFunctionType cost_type,
                    int corr_timeout, double seconds_per_op,
                    float consistency_threshold,
                    int min_consistency_level,
                    int filter_half_kernel,
                    int32 max_pyramid_levels,
                    CorrelationAlgorithm  algorithm = VW_CORRELATION_BM,
                    int collar_size = 0,
                    SemiGlobalMatcher::SgmSubpixelMode sgm_subpixel_mode
                    = SemiGlobalMatcher::SUBPIXEL_LC_BLEND,
                    Vector2i sgm_search_buffer = Vector2i(2,2),
                    size_t memory_limit_mb = 6000,
                    int blob_filter_area = 0,
                    ImageView<PixelMask<float>> * lr_disp_diff = NULL,
                    Vector2i const& region_ul = Vector2i(0, 0),
                    bool  write_debug_images = false) {
    typedef PyramidCorrelationView result_type;
    return result_type(left.impl(), right.impl(), 
                       left_mask.impl(), right_mask.impl(),
                       prefilter_mode, prefilter_width,
                       search_region,
                       kernel_size, cost_type,
                       corr_timeout, seconds_per_op,
                       consistency_threshold, min_consistency_level,
                       filter_half_kernel, max_pyramid_levels,
                       algorithm, collar_size, sgm_subpixel_mode,
                       sgm_search_buffer, memory_limit_mb, blob_filter_area,
                       lr_disp_diff, region_ul, write_debug_images);
  }
  
}} // namespace vw::stereo

#endif//__VW_STEREO_CORRELATION_VIEW_H__
