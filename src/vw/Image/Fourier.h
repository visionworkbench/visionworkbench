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

#include <complex>
#include <vw/InterestPoint/Detector.h>  // TODO: Move get_opencv_wrapper out of here!
#include <vw/Image/ImageResourceImpl.h>
#include <vw/Image/ImageResourceView.h>
#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp" // DEBUG
#include <boost/core/null_deleter.hpp>
  
namespace vw {


//TODO: Split off into a .cc file!
  
/// Take the Discrete Fourier Transform of a VW image and return it in OpenCV format.
template <class T>
void get_dft(ImageViewBase<T> const& input_view, cv::Mat &output_image) {

  cv::Mat I, cv_mask;
  ImageView<vw::uint8> buffer_view;
  vw::ip::get_opencv_wrapper(input_view, I, buffer_view, cv_mask, true);

  //std::cout << "Input view size: rows = " << input_view.rows() <<
  //                            ", cols = " << input_view.cols() << std::endl;
  
  cv::Mat padded;                            //expand input image to optimal size
  //int m = cv::getOptimalDFTSize( I.rows );
  //int n = cv::getOptimalDFTSize( I.cols ); // on the border add zero values
  int m = I.rows; // TODO: Fix padding so it does not affect the results!
  int n = I.cols;
  int pad_bottom = m - I.rows;
  int pad_right  = n - I.cols;
  cv::copyMakeBorder(I, padded, 0, pad_bottom, 0, pad_right, cv::BORDER_CONSTANT, cv::Scalar::all(0));

  //cv::imwrite("/home/smcmich1/data/subpixel/padded.tif", padded);

  //std::cout << "Padded image size: rows = " << padded.rows <<
  //                              ", cols = " << padded.cols << std::endl;

  cv::Mat planes[] = {cv::Mat_<float>(padded), cv::Mat::zeros(padded.size(), CV_32F)};
  cv::merge(planes, 2, output_image);         // Add to the expanded another plane with zeros
  
  cv::dft(output_image, output_image, 0, I.rows);            // this way the result may fit in the source matrix
}

/// Extract the magnitude from a complex image.
inline
void get_magnitude(cv::Mat const& complex, cv::Mat & mag) {

  cv::Mat planes[2];
  cv::split(complex, planes);                   // planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))
  cv::magnitude(planes[0], planes[1], planes[0]);// planes[0] = magnitude
  mag = planes[0];
}

/// Rearrange the quadrants of Fourier image  so that the origin is at the image center
/// - Equivalent to FFTSHIFT in Matlab and numpy.
/// - Set reverse=true to undo an earlier fftshift call.
inline
cv::Mat fftshift(cv::Mat const& in, bool reverse=false) {

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
  cv::Mat q0_in(in, cv::Rect(  0,   0, cxi,   cyi  )); // Top-Left - Create a ROI per quadrant
  cv::Mat q1_in(in, cv::Rect(cxi,   0, w-cxi, cyi  )); // Top-Right
  cv::Mat q2_in(in, cv::Rect(  0, cyi, cxi,   h-cyi)); // Bottom-Left
  cv::Mat q3_in(in, cv::Rect(cxi, cyi, w-cxi, h-cyi)); // Bottom-Right

  cv::Mat q0_out(out, cv::Rect(  0,   0, cxo,   cyo  ));   // Top-Left - Create a ROI per quadrant
  cv::Mat q1_out(out, cv::Rect(cxo,   0, w-cxo, cyo  ));  // Top-Right
  cv::Mat q2_out(out, cv::Rect(  0, cyo, cxo,   h-cyo));  // Bottom-Left
  cv::Mat q3_out(out, cv::Rect(cxo, cyo, w-cxo, h-cyo)); // Bottom-Right

  // Copy all of the opposite quadrants.
  if (!q0_in.empty() && !q3_out.empty()) q0_in.copyTo(q3_out);
  if (!q1_in.empty() && !q2_out.empty()) q1_in.copyTo(q2_out);
  if (!q2_in.empty() && !q1_out.empty()) q2_in.copyTo(q1_out);
  if (!q3_in.empty() && !q0_out.empty()) q3_in.copyTo(q0_out);

  return out;
}

/// Convert the magnitude image into an easy to view format.
inline
void get_pretty_magnitude(cv::Mat &magI, bool do_fftshift=true) {

  magI += cv::Scalar::all(1); // Switch to logarithmic scale
  cv::log(magI, magI);
  if (do_fftshift)
    magI = fftshift(magI);
}

/// Get the magnitude of a frequency domain image and write it to disk.
inline
void save_mag_from_ft(cv::Mat const& input, std::string const& path, bool do_fftshift=true) {
  cv::Mat mag;
  get_magnitude(input, mag);
  get_pretty_magnitude(mag, do_fftshift);
  boost::shared_ptr<cv::Mat> ocv_ptr(&mag, boost::null_deleter());
  ImageResourceView<float> ocv_view(new ImageResourceOpenCV(ocv_ptr));
  write_image(path, ocv_view);
}

/// Increase the size of a frequency domain image with zero padding,
///  keeping everything in the correct place.
inline
cv::Mat pad_fourier_transform(cv::Mat const& input, int new_width, int new_height) {

  if ( (new_width < input.cols) || (new_height < input.rows) ) {
    vw_throw(ArgumentErr() << "pad_fourier_transform cannot shrink the image!\n");
  }

  if ((new_width == input.cols) && (new_height == input.rows))
    return input; // No change case

  cv::Mat temp(new_height, new_width, input.type());
  temp = 0.0;

  cv::Mat in_copy = fftshift(input); // Move zero frequency to the center to make things easier.

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

  // Reverse fftshift to get things back how they were.
  return fftshift(temp, true)*scale;
}

/// Applies a raised cosine filter across an entire image.
/// - The filter width is the same as the image width so no convolution is performed.
/// - Inputs must be type CV_32FC2
/// - Set center_zero_freq if you fftshifted the low frequencies into the center of the image.
inline void apply_raised_cosine_filter(cv::Mat &input, const float beta,
                                       bool center_zero_freq=false) {

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
    else // Low frequencies in the corners.
      y = std::min(r, input.rows-1-r);//;
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
      else // Low frequency in the corners
        x = std::min(c, input.cols-1-c);//
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

#endif // __VW_IMAGE_FOURIER_H__
