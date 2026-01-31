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

#ifndef __VW_STEREO_PHASESUBPIXEL_VIEW__
#define __VW_STEREO_PHASESUBPIXEL_VIEW__

#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/PixelMask.h>
#include <opencv2/core.hpp>

// Implement the subpixel phase correlation algorithm from the following paper:
// Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup,
// "Efficient subpixel image registration algorithms," Opt. Lett. 33, 156-158 (2008).
// Adapted from the MATLAB implementation available here:
// https://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation

namespace vw {
namespace stereo {

/// Use matrix multiplication to upsample a DFT in only a small region.
///  It is usually faster than doing the equivalent series of steps:
///   1) pad_fourier_transform(input, input.rows()*upscale, input.cols()*upscale)
///      dimension. fftshift(inverse) to bring the center of the image to (0,0).
///   2) dft(upsampled_image)
///   3) crop(dft_result, col_offset, row_offset, upsampled_width, upsampled_height)
cv::Mat partial_upsample_dft(cv::Mat const& input, int upsampled_height, int upsampled_width, 
                             int upscale, int row_offset=0, int col_offset=0);

/// Compute the subpixel translation between two images using a two-pass
/// frequency based method. The images must be the same size. Maximum accuracy
/// is 1/subpixel_accuracy. Images are converted to float internally, so use
/// ImageViewRef<float> for efficiency.
void phase_correlation_subpixel(ImageViewRef<float> const& left_image,
                                ImageViewRef<float> const& right_image,
                                Vector2f &offset,
                                int subpixel_accuracy=10,
                                bool debug = false);

/// Update the values in disparity_map according to region_of_interest.
/// - This function is set up to work with the PyramidSubpixelView class.
/// - use_second_refinement improves results by repeating the computation.
/// - Images must be float for phase correlation.
void subpixel_phase_2d(ImageView<PixelMask<Vector2f> > &disparity_map,
                       ImageView<float> const& left_image,
                       ImageView<float> const& right_image,
                       int32  kern_width, int32 kern_height,
                       BBox2i region_of_interest,
                       int  subpixel_accuracy = 20,
                       bool use_second_refinement = true);


}} // namespace vw::stereo

#endif // __VW_STEREO_PHASESUBPIXEL_VIEW__
