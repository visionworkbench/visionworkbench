// __BEGIN_LICENSE__
//  Copyright (c) 2006-2026, United States Government as represented by the
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

#include <vw/Image/Fourier.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/ImageChannels.h>
#include <vw/Image/ImageIO.h>
#include <vw/Image/ImageResourceView.h>
#include <vw/Image/ImageResourceImpl.h>
#include <vw/FileIO/DiskImageResource.h>
#include <opencv2/imgproc.hpp>
#include <boost/core/null_deleter.hpp>

namespace vw {

// Helper function to convert VW image to OpenCV format
void fourier_get_opencv_wrapper(ImageViewRef<float> const& input_image,
                                cv::Mat & cv_image,
                                ImageView<vw::uint8> &image_buffer,
                                cv::Mat & cv_mask,
                                bool normalize = true) {

  // Rasterize the input image
  ImageView<float> input_buffer = input_image.impl();

  if (normalize) // Convert to uint8 with 2%-98% intensity scaling
    percentile_scale_convert(input_buffer, image_buffer, 0.02, 0.98);
  else {
    double standard_min = ChannelRange<float>::min();
    double standard_max = ChannelRange<float>::max();
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

void get_dft(ImageViewRef<float> const& input_view, cv::Mat &output_image) {

  cv::Mat I, cv_mask;
  ImageView<vw::uint8> buffer_view;
  fourier_get_opencv_wrapper(input_view, I, buffer_view, cv_mask, true);

  cv::Mat padded;
  int m = I.rows;
  int n = I.cols;
  int pad_bottom = m - I.rows;
  int pad_right  = n - I.cols;
  cv::copyMakeBorder(I, padded, 0, pad_bottom, 0, pad_right, cv::BORDER_CONSTANT, cv::Scalar::all(0));

  cv::Mat planes[] = {cv::Mat_<float>(padded), cv::Mat::zeros(padded.size(), CV_32F)};
  cv::merge(planes, 2, output_image);

  cv::dft(output_image, output_image, 0, I.rows);
}

void get_magnitude(cv::Mat const& complex, cv::Mat & mag) {

  cv::Mat planes[2];
  cv::split(complex, planes);
  cv::magnitude(planes[0], planes[1], planes[0]);
  mag = planes[0];
}

cv::Mat fftshift(cv::Mat const& in, bool reverse) {

  // Allocate the output image, same size as the input image.
  const int h = in.rows;
  const int w = in.cols;
  cv::Mat out(h, w, in.type());

  // Set up the dividing indices for the input and output images, handling odd sizes.
  int cxo = in.cols/2;
  int cyo = in.rows/2;

  int cxi = cxo;
  int cyi = cyo;
  if (in.cols % 2 != 0) ++cxi;
  if (in.rows % 2 != 0) ++cyi;
  
  if (reverse) {
    std::swap(cxo, cxi);
    std::swap(cyo, cyi);
  }

  // Create quadrants for the input and output images
  cv::Mat q0_in(in, cv::Rect(  0,   0, cxi,   cyi  ));
  cv::Mat q1_in(in, cv::Rect(cxi,   0, w-cxi, cyi  ));
  cv::Mat q2_in(in, cv::Rect(  0, cyi, cxi,   h-cyi));
  cv::Mat q3_in(in, cv::Rect(cxi, cyi, w-cxi, h-cyi));

  cv::Mat q0_out(out, cv::Rect(  0,   0, cxo,   cyo  ));
  cv::Mat q1_out(out, cv::Rect(cxo,   0, w-cxo, cyo  ));
  cv::Mat q2_out(out, cv::Rect(  0, cyo, cxo,   h-cyo));
  cv::Mat q3_out(out, cv::Rect(cxo, cyo, w-cxo, h-cyo));

  // Copy all of the opposite quadrants.
  if (!q0_in.empty() && !q3_out.empty()) q0_in.copyTo(q3_out);
  if (!q1_in.empty() && !q2_out.empty()) q1_in.copyTo(q2_out);
  if (!q2_in.empty() && !q1_out.empty()) q2_in.copyTo(q1_out);
  if (!q3_in.empty() && !q0_out.empty()) q3_in.copyTo(q0_out);

  return out;
}

void get_pretty_magnitude(cv::Mat &magI, bool do_fftshift) {

  magI += cv::Scalar::all(1);
  cv::log(magI, magI);
  if (do_fftshift)
    magI = fftshift(magI);
}

cv::Mat pad_fourier_transform(cv::Mat const& input, int new_width, int new_height) {

  if ((new_width < input.cols) || (new_height < input.rows))
    vw_throw(ArgumentErr() << "pad_fourier_transform cannot shrink the image!\n");

  if ((new_width == input.cols) && (new_height == input.rows))
    return input;

  cv::Mat temp(new_height, new_width, input.type());
  temp = 0.0;

  cv::Mat in_copy = fftshift(input);

  int in_center_x = floor(input.cols/2)+1;
  int in_center_y = floor(input.rows/2)+1;

  int out_center_x = floor(new_width/2)+1;
  int out_center_y = floor(new_height/2)+1;

  int cdx = out_center_x - in_center_x;
  int cdy = out_center_y - in_center_y;

  float scale = static_cast<float>(new_width*new_height)/
                static_cast<float>(input.cols*input.rows);

  cv::Mat out_ref(temp, cv::Rect(cdx, cdy, input.cols, input.rows));

  in_copy.copyTo(out_ref);

  return fftshift(temp, true)*scale;
}

void apply_raised_cosine_filter(cv::Mat &input, const float beta, bool center_zero_freq) {

  const int COMPLEX_TYPE_CV = CV_32FC2;
  typedef std::complex<float> c_type;

  const float Y = input.rows;
  const float X = input.cols;

  const float y_center = Y / 2.0;
  const float x_center = X / 2.0;
  const float y_cutoff = Y * (0.5 - beta);
  const float x_cutoff = X * (0.5 - beta);
  const float y_coeff  = M_PI / (2.0*beta*Y);
  const float x_coeff  = M_PI / (2.0*beta*X);

  float y, x;
  for (int r=0; r<input.rows; ++r) {

    // Compute the Y weight for this row.
    if (center_zero_freq)
      y = fabs(r - y_center);
    else
      y = std::min(r, input.rows-1-r);
    float vert_weight = 0;
    if (y < y_cutoff)
      vert_weight = 1.0;
    else {
      if (y < y_center) {
        float temp = cos(y_coeff * (y - y_cutoff));
        vert_weight = temp*temp;
      }
    }

    for (int c=0; c<input.cols; ++c) {

      // Compute the X weight for this column.
      if (center_zero_freq)
        x = fabs(c - x_center);
      else
        x = std::min(c, input.cols-1-c);
      float horiz_weight = 0;
      if (x < x_cutoff)
        horiz_weight = 1.0;
      else {
        if (x < x_center) {
          float temp = cos(x_coeff * (x - x_cutoff));
          horiz_weight = temp*temp;
        }
      }
      input.at<c_type>(r,c) *= (vert_weight*horiz_weight);
    }
  }
  return;
}

} // namespace vw
