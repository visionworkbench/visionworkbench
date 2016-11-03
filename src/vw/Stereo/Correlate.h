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
#include <vw/Image/Interpolation.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/PixelMask.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/LinearAlgebra.h>

namespace vw {
namespace stereo {

  int
  adjust_weight_image(ImageView<float> &weight,
                      ImageView<PixelMask<Vector2f> > const& disparity_map_patch,
                      ImageView<float> const& weight_template);

  namespace detail {
    ImageView<float>
    compute_spatial_weight_image(int32 kern_width, int32 kern_height,
                                 float two_sigma_sqr);
  }

  /// Check for consistency between a left-to-right and right-to-left
  ///  pair of stereo results.  The results are expected to be opposite in
  ///  sign but equal in magnitude.
  /// - Set aligned_images if the RL image has been shifted so that
  ///   pixels in RL should be in the same location as in LR.
  template <class ImageT1, class ImageT2>
  void cross_corr_consistency_check( ImageViewBase<ImageT1> const& l2r,
                                     ImageViewBase<ImageT2> const& r2l,
                                     float cross_corr_threshold, 
                                     bool aligned_images = false,
                                     bool verbose = false ) {
    int32 l2r_rows = l2r.impl().rows(), l2r_cols = l2r.impl().cols(),
          r2l_rows = r2l.impl().rows(), r2l_cols = r2l.impl().cols();
    size_t count = 0, match_count = 0;

    if (verbose)
      vw_out(DebugMessage, "stereo") << "\tCrosscorr threshold: "
                                    << cross_corr_threshold << "\n";
    VW_DEBUG_ASSERT( cross_corr_threshold >= 0,
                     ArgumentErr() << "cross_corr_consistency_check: the threshold is less than 0." );

    typename ImageT1::pixel_accessor l2r_row = l2r.impl().origin();
    for ( int32 r = 0; r < l2r_rows; ++r ) {
      typename ImageT1::pixel_accessor l2r_col = l2r_row;
      for ( int32 c = 0; c < l2r_cols; ++c ) {

        // Under our regular workflow the corresponding disparity scores will be
        //  at different pixels, but that is not always the case.
        int32 r2l_x = c;
        int32 r2l_y = r;
        if (!aligned_images) {
          r2l_x += (*l2r_col)[0];
          r2l_y += (*l2r_col)[1];
        }

        if ( r2l_x < 0 || r2l_x >= r2l_cols ||
             r2l_y < 0 || r2l_y >= r2l_rows ) {
          // Verify that we are in image bounds
          invalidate( *l2r_col );
        } else if ( !is_valid( *l2r_col ) ||
                    !is_valid( r2l.impl()(r2l_x,r2l_y ) ) ) {
          // Verify that both are not missing
          invalidate( *l2r_col );
        } else if ( cross_corr_threshold >= fabs((*l2r_col)[0] + r2l.impl()(r2l_x,r2l_y)[0]) &&
                    cross_corr_threshold >= fabs((*l2r_col)[1] + r2l.impl()(r2l_x,r2l_y)[1]) ) {
          // Actually check the correlation consistency
          //
          // Since the hdisp for the R2L and L2R buffers will be opposite
          // in sign, we determine their similarity by *summing* them, rather
          // than differencing them as you might expect.
          count++;
          match_count++;
        } else {
          match_count++;
          invalidate( *l2r_col );
        }

        l2r_col.next_col();
      }
      l2r_row.next_row();
    }

    if (verbose)
      vw_out(DebugMessage, "stereo") << "\tCross-correlation retained " << count
                                     << " / " << match_count << " matches ("
                                     << ((float)count/match_count*100) << " percent).\n";
  }

#include <vw/Stereo/Correlate.tcc>


}} // namespace vw::stereo

#endif // __VW_STEREO_CORRELATOR_H__
