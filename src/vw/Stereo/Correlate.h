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


#ifndef __VW_STEREO_CORRELATE_H__
#define __VW_STEREO_CORRELATE_H__

#include <vw/Core/Exception.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Log.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/PixelMask.h>
#include <vw/Math/Vector.h>

namespace vw { namespace stereo {

  /// ???
  int
  adjust_weight_image(ImageView<float> &weight,
                      ImageView<PixelMask<Vector2f>> const& disparity_map_patch,
                      ImageView<float> const& weight_template);


  namespace detail {
    ImageView<float>
    compute_spatial_weight_image(int32 kern_width, int32 kern_height,
                                 float two_sigma_sqr);
  }

  /// Check for consistency between a left-to-right and right-to-left
  /// pair of stereo results.  The results are expected to be opposite in
  /// sign but equal in magnitude. Save the absolute discrepancy
  /// among the two disparities in lr_disp_diff, if provided.
  template <class ImageT1, class ImageT2>
  void cross_corr_consistency_check(ImageViewBase<ImageT1> const& l2r,
                                    ImageViewBase<ImageT2> const& r2l,
                                    float cross_corr_threshold,
                                    ImageView<PixelMask<float>> * lr_disp_diff = NULL,
                                    Vector2i ul_corner_offset = Vector2i(),
                                    bool verbose = false) {
    int32 l2r_rows = l2r.impl().rows(), l2r_cols = l2r.impl().cols(),
      r2l_rows = r2l.impl().rows(), r2l_cols = r2l.impl().cols();
    size_t count = 0, match_count = 0;

    if (verbose)
      vw_out(VerboseDebugMessage, "stereo") << "\tCrosscorr threshold: "
                                            << cross_corr_threshold << "\n";
    VW_DEBUG_ASSERT(cross_corr_threshold >= 0.0,
                    ArgumentErr() << "cross_corr_consistency_check: the threshold is less than 0.");
  
    typename ImageT1::pixel_accessor l2r_row = l2r.impl().origin();
    for ( int32 r = 0; r < l2r_rows; ++r ) {
      typename ImageT1::pixel_accessor l2r_col = l2r_row;
      for ( int32 c = 0; c < l2r_cols; ++c ) {

        // The corresponding disparity scores will be at different
        // pixels according to the disparity values.
        int32 r2l_x = c + (*l2r_col)[0];
        int32 r2l_y = r + (*l2r_col)[1];

        if ( r2l_x < 0 || r2l_x >= r2l_cols ||
             r2l_y < 0 || r2l_y >= r2l_rows ) {
          // Verify that we are in image bounds
          invalidate( *l2r_col );
        } else if ( !is_valid( *l2r_col ) ||
                    !is_valid( r2l.impl()(r2l_x,r2l_y ) ) ) {
          // Verify that both are not missing
          invalidate(*l2r_col);
        } else {

          float disp_diff = std::max(fabs((*l2r_col)[0] + r2l.impl()(r2l_x,r2l_y)[0]),
                                     fabs((*l2r_col)[1] + r2l.impl()(r2l_x,r2l_y)[1]));
          if (cross_corr_threshold >= disp_diff) {
          
            // Actually check the correlation consistency
            //
            // Since the hdisp for the R2L and L2R buffers will be opposite
            // in sign, we determine their similarity by *summing* them, rather
            // than differencing them as you might expect.
            count++;
            match_count++;

            // Save the difference if a buffer is provided. At pixels where nothing is saved
            // the original nodata value will be kept.
            if (lr_disp_diff != NULL) 
              (*lr_disp_diff)(c + ul_corner_offset[0], r + ul_corner_offset[1])
                = PixelMask<float>(disp_diff); // it becomes a valid pixel with this value
          
          } else {
            match_count++;
            invalidate(*l2r_col);
          }
        }
      
        l2r_col.next_col();
      }
      l2r_row.next_row();
    }

    if (verbose)
      vw_out(VerboseDebugMessage, "stereo") << "\tCross-correlation retained " << count
                                            << " / " << match_count << " matches ("
                                            << ((float)count/match_count*100) << " percent).\n";
  }

  /// Fast affine-EM implementation
  /// In this version we don't keep around future research ideas
  /// since they are slow.
  void subpixel_optimized_affine_2d_EM(ImageView<PixelMask<Vector2f>> & disparity_map,
                                       ImageView<float> const& left_image,
                                       ImageView<float> const& right_image,
                                       int32  kern_width, int32 kern_height,
                                       BBox2i region_of_interest,
                                       bool   do_horizontal_subpixel,
                                       bool   do_vertical_subpixel,
                                       bool /*verbose*/ );

  /// Research-y implementation of affine EM correlation.
  void
  subpixel_correlation_affine_2d_EM(ImageView<PixelMask<Vector2f> > &disparity_map,
                                    ImageView<float> const& left_image,
                                    ImageView<float> const& right_image,
                                    int32 kern_width, int32 kern_height,
                                    BBox2i region_of_interest,
                                    bool do_horizontal_subpixel,
                                    bool do_vertical_subpixel,
                                    bool /*verbose*/ );

  /// Affine subpixel correlation function
  /// - Similar to the _EM function but about 5X faster and sometimes less accurate.
  void
  subpixel_optimized_affine_2d(ImageView<PixelMask<Vector2f> > &disparity_map,
                              ImageView<float> const& left_image,
                              ImageView<float> const& right_image,
                              int32  kern_width, int32 kern_height,
                              BBox2i region_of_interest,
                              bool   do_horizontal_subpixel,
                              bool   do_vertical_subpixel,
                              bool /*verbose*/ );

  /// Lucas-Kanade subpixel correlation function
  /// - So far this one does not work that well.
  void
  subpixel_optimized_LK_2d(ImageView<PixelMask<Vector2f>> & disparity_map,
                          ImageView<float> const& left_image,
                          ImageView<float> const& right_image,
                          int32  kern_width, int32 kern_height,
                          BBox2i region_of_interest,
                          bool   do_horizontal_subpixel,
                          bool   do_vertical_subpixel);

}} // namespace vw::stereo

#endif // __VW_STEREO_CORRELATOR_H__
