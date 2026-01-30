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
// Functions dealing with Fourier transforms.
// - To make things easier, these all use OpenCV code.
//
//
#ifndef __VW_IMAGE_FOURIER_H__
#define __VW_IMAGE_FOURIER_H__

#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/ImageIO.h>

#include <complex>
#include <vw/Image/ImageResourceImpl.h>
#include <vw/Image/ImageResourceView.h>
#include <vw/FileIO/DiskImageResource.h>
#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp" // DEBUG
#include <boost/core/null_deleter.hpp>
  
namespace vw {

// Helper function to convert VW image to OpenCV format
// (Moved here from InterestPoint/Detector.h to avoid that dependency)
template <class ViewT>
inline void fourier_get_opencv_wrapper(ImageViewBase<ViewT> const& input_image,
                                       cv::Mat & cv_image,
                                       ImageView<vw::uint8> &image_buffer,
                                       cv::Mat & cv_mask,
                                       bool normalize = true) {
  // Rasterize the input image
  ImageView<typename ViewT::pixel_type> input_buffer = input_image.impl();

  if (normalize) // Convert to uint8 with 2%-98% intensity scaling
    percentile_scale_convert(input_buffer, image_buffer, 0.02, 0.98);
  else {
    double standard_min = ChannelRange<typename ViewT::pixel_type>::min();
    double standard_max = ChannelRange<typename ViewT::pixel_type>::max();
    image_buffer = pixel_cast_rescale<vw::uint8>(clamp(input_buffer, standard_min, standard_max));
  }

  // Create OpenCV wrapper for the buffer
  int     cv_data_type = CV_8UC1;
  void*   raw_data_ptr = reinterpret_cast<void*>(image_buffer.data());
  size_t  pixel_size   = sizeof(vw::uint8);
  size_t  step_size    = image_buffer.cols() * pixel_size;

  cv_image = cv::Mat(image_buffer.rows(), image_buffer.cols(),
                     cv_data_type, raw_data_ptr, step_size);

  // No mask handling needed for Fourier transforms
  cv_mask = cv::Mat();
}

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

/// Get the magnitude of a frequency domain image and write it to disk.
void save_mag_from_ft(cv::Mat const& input, std::string const& path, bool do_fftshift=true);

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
