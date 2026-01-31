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

// \file Fourier.h
//
// Functions dealing with Fourier transforms. To make things easier, these all
// use OpenCV code.

#ifndef __VW_IMAGE_FOURIER_H__
#define __VW_IMAGE_FOURIER_H__

#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/ImageViewBase.h>
#include <opencv2/core.hpp>

namespace vw {

/// Take the Discrete Fourier Transform of a VW image and return it in OpenCV format.
/// Input is converted to float internally, so use ImageViewRef<float> for efficiency.
void get_dft(ImageViewRef<float> const& input_view, cv::Mat &output_image);

/// Extract the magnitude from a complex image.
void get_magnitude(cv::Mat const& complex, cv::Mat & mag);

/// Rearrange the quadrants of Fourier image so that the origin is at the image center
/// - Equivalent to FFTSHIFT in Matlab and numpy.
/// - Set reverse=true to undo an earlier fftshift call.
cv::Mat fftshift(cv::Mat const& in, bool reverse=false);

/// Convert the magnitude image into an easy to view format.
void get_pretty_magnitude(cv::Mat &magI, bool do_fftshift=true);

/// Increase the size of a frequency domain image with zero padding,
///  keeping everything in the correct place.
cv::Mat pad_fourier_transform(cv::Mat const& input, int new_width, int new_height);

/// Applies a raised cosine filter across an entire image.
/// - The filter width is the same as the image width so no convolution is performed.
/// - Inputs must be type CV_32FC2
/// - Set center_zero_freq if you fftshifted the low frequencies into the center of the image.
void apply_raised_cosine_filter(cv::Mat &input, const float beta, bool center_zero_freq=false);

} // namespace vw

#endif // __VW_IMAGE_FOURIER_H__
