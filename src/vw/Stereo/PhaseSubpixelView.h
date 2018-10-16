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
#include <vw/Image/Fourier.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Stereo/Correlate.h>
#include <vw/Stereo/PreFilter.h>
#include <vw/Stereo/CostFunctions.h>
#include <vw/Stereo/Correlation.h>


/**

Implement the subpixel phase correlation algorithm from the following paper:

Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
"Efficient subpixel image registration algorithms," Opt. Lett. 33, 156-158 (2008).

Adapted from the MATLAB implementation available here:

https://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation

*/

namespace vw {
namespace stereo {

/// Use matrix multiplication to upsample a DFT in only a small region.
///  It is usually faster than doing the equivalent series of steps:
///   1) pad_fourier_transform(input, input.rows()*upscale, input.cols()*upscale)
///      dimension. fftshift(inverse) to bring the center of the image to (0,0).
///   2) dft(upsampled_image)
///   3) crop(dft_result, col_offset, row_offset, upsampled_width, upsampled_height)
inline
cv::Mat partial_upsample_dft(cv::Mat const& input, int upsampled_height, int upsampled_width, 
                             int upscale, int row_offset=0, int col_offset=0) {

  // Generate some linear row vectors
  const int COMPLEX_TYPE_CV = CV_32FC2;
  typedef std::complex<float> c_type;

  cv::Mat col_vector   (1, input.cols,       COMPLEX_TYPE_CV);
  cv::Mat row_vector   (1, input.rows,       COMPLEX_TYPE_CV);
  cv::Mat col_up_vector(1, upsampled_width,  COMPLEX_TYPE_CV);
  cv::Mat row_up_vector(1, upsampled_height, COMPLEX_TYPE_CV);

  for (int i=0; i<input.cols; ++i)
    col_vector.at<c_type>(0,i) = static_cast<c_type>(i);
  for (int i=0; i<input.rows; ++i)
    row_vector.at<c_type>(0,i) = static_cast<c_type>(i);
  for (int i=0; i<upsampled_width; ++i)
    col_up_vector.at<c_type>(0,i) = static_cast<c_type>(i);
  for (int i=0; i<upsampled_height; ++i)
    row_up_vector.at<c_type>(0,i) = static_cast<c_type>(i);

  const c_type neg_i(0, -1);
  float two_pi = 2.0*M_PI;

  cv::Mat col_kernel, row_kernel;
  cv::Mat temp_vec = fftshift(col_vector, true).t();
  cv::Mat v1 = temp_vec - floor(input.cols/2);
  cv::Mat v2 = col_up_vector - col_offset;
  
  cv::Mat dummy;
  cv::gemm(v1,v2,1.0, dummy, 0.0, col_kernel); //m  = v1*v2;
  
  c_type constant = neg_i*two_pi/static_cast<float>(input.cols*upscale);
  for (int i=0; i<col_kernel.rows; ++i)
    for (int j=0; j<col_kernel.cols; ++j)
      col_kernel.at<c_type>(i,j) = std::exp(col_kernel.at<c_type>(i,j)*constant);

  temp_vec = fftshift(row_vector, true);
  v1 = row_up_vector.t() - row_offset;
  v2 = temp_vec - floor(input.rows/2);
  
  cv::gemm(v1,v2,1.0, dummy, 0.0, row_kernel); //m  = v1*v2;
  
  constant = neg_i*two_pi/static_cast<float>(input.rows*upscale);
  for (int i=0; i<row_kernel.rows; ++i)
    for (int j=0; j<row_kernel.cols; ++j)
      row_kernel.at<c_type>(i,j) = std::exp(row_kernel.at<c_type>(i,j)*constant);

  //save_mag_from_ft(row_kernel,  "/home/smcmich1/data/subpixel/row_kernel.tif", false);
  //save_mag_from_ft(col_kernel,  "/home/smcmich1/data/subpixel/col_kernel.tif", false);

  // These muliplications are the slowest part of the phase correlation method.
  cv::Mat o1  = row_kernel*input;
  cv::Mat out = o1*col_kernel;

  //cv::Mat temp = complex_multiply(row_kernel, input);
  //cv::Mat out  = complex_multiply(temp, col_kernel);

  return out;
}

/// Compute the subpixel translation between two images using a two-pass frequency based method.
/// - The images must be the same size!
/// - Maximum accuracy is 1/subpixel_accuracy
template <class T1, class T2>
void phase_correlation_subpixel(ImageViewBase<T1> const& left_image,
                                ImageViewBase<T2> const& right_image,
                                Vector2f &offset,
                                int subpixel_accuracy=10,
                                bool debug = false
                               ) {

  if (left_image.get_size() != right_image.get_size()) {
    vw_throw( ArgumentErr() << "phase_correlation_subpixel requires images to be the same size!\n" );
  }

  // Fourier transform of the input images.
  // TODO: Use padding for an optimal DFT size?
  cv::Mat complexI_left, complexI_right;
  get_dft(left_image,  complexI_left);
  get_dft(right_image, complexI_right);

  // Some papers suggest filtering out high frequency image content prior to correlation
  //  but that does not seem to help in all cases.
  //const float BETA = 0.35;
  ////save_mag_from_ft(complexI_left,  "/home/smcmich1/data/subpixel/magI_left_prefilter.tif");
  ////save_mag_from_ft(complexI_right,  "/home/smcmich1/data/subpixel/magI_right_prefilter.tif");
  //apply_raised_cosine_filter(complexI_left,  BETA);
  //apply_raised_cosine_filter(complexI_right, BETA);
  
  // The first pass will try to find the best shift location at a low resolution, 
  // then the second pass will try to refine the result nearby that location.
  // By doing this we avoid doing full resolution computations over the entire image.

  // Compute convolution of the two images.
  cv::Mat initial_conj;
  cv::mulSpectrums(complexI_left, complexI_right, initial_conj, 0, true);

  // Pad the results to we get higher DFT accuracy.
  int pad_factor = subpixel_accuracy; // Controls maximum subpixel accuracy.

  const int INITIAL_PAD_FACTOR = 2;
  int padded_width  = INITIAL_PAD_FACTOR*initial_conj.cols;
  int padded_height = INITIAL_PAD_FACTOR*initial_conj.rows;
  cv::Mat padded_conj = pad_fourier_transform(initial_conj, padded_width, padded_height);

  // Inverse FFT to get back to image coordinates.
  cv::Mat conv;
  cv::dft(padded_conj, conv, cv::DFT_INVERSE + cv::DFT_REAL_OUTPUT+ cv::DFT_SCALE, 0);

  // Find the peak.
  int width  = conv.cols;
  int height = conv.rows;
  double maxVal;
  cv::Point maxLoc;
  cv::minMaxLoc(conv, NULL, &maxVal, NULL, &maxLoc);

  // Convert peak location back to the input pixel coordinates.
  float initial_shift_x = (maxLoc.x<width /2) ? (maxLoc.x) : (maxLoc.x-width );
  float initial_shift_y = (maxLoc.y<height/2) ? (maxLoc.y) : (maxLoc.y-height);
  initial_shift_x /= static_cast<float>(INITIAL_PAD_FACTOR);
  initial_shift_y /= static_cast<float>(INITIAL_PAD_FACTOR);

  if (debug) {

    std::cout << "padded_conj.size() = " << padded_conj.size() << std::endl;    
    std::cout << "maxLoc = " << maxLoc << std::endl;
    std::cout << "initial_shift_x = " << initial_shift_x << std::endl;
    std::cout << "initial_shift_y = " << initial_shift_y << std::endl;

    save_mag_from_ft(complexI_left,  "/home/smcmich1/data/subpixel/magI_left.tif");
    save_mag_from_ft(complexI_right, "/home/smcmich1/data/subpixel/magI_right.tif");
    save_mag_from_ft(initial_conj, "/home/smcmich1/data/subpixel/initial_conj.tif");
    save_mag_from_ft(padded_conj, "/home/smcmich1/data/subpixel/padded_conj.tif");

    boost::shared_ptr<cv::Mat> ocv_ptr(&conv, boost::null_deleter());
    ImageResourceView<float> ocv_view(new ImageResourceOpenCV(ocv_ptr));
    write_image( "/home/smcmich1/data/subpixel/conv.tif", ocv_view);
  }

  // End of the first pass, stop here if the output resolution is low.
  if (pad_factor <= 2) {
    offset[0] = initial_shift_x;
    offset[1] = initial_shift_y;
    return;
  }

  // No do a second pass to improve the answer resolution.

  // The size of the region (in units of input pixels) around the low-resolution
  // peak where we will search for the high resolution peak.
  const float UPSAMPLE_REGION_FACTOR = 1.5;

  // Compute the location of interest to be upsampled.
  float shift_x = round(initial_shift_x*pad_factor)/pad_factor; 
  float shift_y = round(initial_shift_y*pad_factor)/pad_factor; 
  float dft_shift = floor(ceil(pad_factor*UPSAMPLE_REGION_FACTOR)/2); //% Center of output array at dftshift+1

  // Use the matrix multiply trick to upsample the region of interest.
  int upsampled_height = ceil(pad_factor*UPSAMPLE_REGION_FACTOR);
  int upsampled_width  = ceil(pad_factor*UPSAMPLE_REGION_FACTOR);
  cv::Mat new_conj;
  cv::mulSpectrums(complexI_right, complexI_left, new_conj, 0, true); // Left/Right order is reversed here.
  cv::Mat partial_upsampled = partial_upsample_dft(new_conj, 
                                                   upsampled_height,
                                                   upsampled_width,
                                                   pad_factor,
                                                   dft_shift-shift_y*pad_factor,
                                                   dft_shift-shift_x*pad_factor);

  // Find the peak
  cv::Mat CC;
  get_magnitude(partial_upsampled, CC);
  cv::minMaxLoc(CC, NULL, &maxVal, NULL, &maxLoc); // TODO: magnitude(conv)??

  // Convert result into the final offset value
  maxLoc.y = maxLoc.y - dft_shift;
  maxLoc.x = maxLoc.x - dft_shift;
  shift_y = shift_y + static_cast<float>(maxLoc.y)/static_cast<float>(pad_factor);
  shift_x = shift_x + static_cast<float>(maxLoc.x)/static_cast<float>(pad_factor);

  offset[0] = shift_x;
  offset[1] = shift_y;

  if (debug) {
    std::cout << "partial_upsampled.size() = " << partial_upsampled.size() << std::endl;
    std::cout << "CC.type = " << CC.type() << std::endl;
    std::cout << "maxLoc = " << maxLoc << std::endl;
    std::cout << "shift_x = " << shift_x << std::endl;
    std::cout << "shift_y = " << shift_y << std::endl;

    save_mag_from_ft(partial_upsampled, "/home/smcmich1/data/subpixel/partial_upsampled.tif", false);
  }
}




/// Update the values in disparity_map according to region_of_interest.
/// - This function is set up to work with the PyramidSubpixelView class.
/// - use_second_refinement improves results by repeating the computation.
template<class ChannelT> void
subpixel_phase_2d(ImageView<PixelMask<Vector2f> > &disparity_map,
                  ImageView<ChannelT> const& left_image,
                  ImageView<ChannelT> const& right_image,
                  int32  kern_width, int32 kern_height,
                  BBox2i region_of_interest,
                  int  subpixel_accuracy = 20,
                  bool use_second_refinement = true) {

  VW_ASSERT( disparity_map.cols() == left_image.cols() &&
             disparity_map.rows() == left_image.rows(),
             ArgumentErr() << "subpixel_correlation: left image and "
                            << "disparity map do not have the same dimensions.");

  // Interpolated right image in case we go out of bounds.
  InterpolationView<EdgeExtensionView<ImageView<ChannelT>, ZeroEdgeExtension>, BilinearInterpolation> right_interp_image =
         interpolate(right_image, BilinearInterpolation(), ZeroEdgeExtension());

  // This is the maximum number of pixels that the solution can be
  // adjusted by subpixel refinement.
  // - Until we trust this method more, keep the value low!
  float SUBPIXEL_MAX_TRANSLATION = 3.0;//kern_width/2;

  const int32 kern_half_height    = kern_height/2;
  const int32 kern_half_width     = kern_width /2;

  // Iterate over all of the pixels in the disparity map except for the outer edges.
  for ( int32 y = std::max(region_of_interest.min().y()-1,kern_half_height);
              y < std::min(left_image.rows()-kern_half_height,
                           region_of_interest.max().y()+1); ++y) {

    for (int32 x = std::max(region_of_interest.min().x()-1,kern_half_width);
               x < std::min(left_image.cols()-kern_half_width,
                            region_of_interest.max().x()+1); ++x) {

      // Skip over pixels for which we have no initial disparity estimate
      // - There is no check that the nearby image regions are valid, hopefully
      //   if they are not then the disparity value is also invalid.
      if ( !is_valid(disparity_map(x,y)) )
        continue;

      // The window in the left image is the current pixel surrounded
      //  by the size of the kernel.  The window in the right image is
      //  the same region offset by the expected offset.
      BBox2i current_window(x-kern_half_width, y-kern_half_height,
                            kern_width, kern_height);
      BBox2i right_window = current_window + Vector2i(disparity_map(x,y)[0], disparity_map(x,y)[1]);

      // Compute the image patches
      ImageView<ChannelT> left_image_patch  = crop(left_image,         current_window);
      ImageView<ChannelT> right_image_patch = crop(right_interp_image, right_window);

      //write_image("/home/smcmich1/data/subpixel/left_patch.tif",  left_image_patch );
      //write_image("/home/smcmich1/data/subpixel/right_patch.tif", right_image_patch);

      // We are just solving for a simple translation vector
      Vector2f d;
      bool debug = false;
      int initial_subpixel_accuracy = subpixel_accuracy;
      if (use_second_refinement) // The first pass can be lower resolution.
        initial_subpixel_accuracy /= 2;
      phase_correlation_subpixel(left_image_patch, right_image_patch,
                                 d, initial_subpixel_accuracy, debug);

      if (use_second_refinement) {
        // Shift the right crop by the computed offset, then re-run
        // phase correlation to get a final offset.
        // - This improves the results at the cost of taking twice as long.

        ImageView<ChannelT> shift_right_crop = crop(translate( right_image,
                                                              d[0], d[1],
                                                              ZeroEdgeExtension(),
                                                              BicubicInterpolation()),
                                                    right_window
                                                  );
        //write_image("/home/smcmich1/data/subpixel/right_patch_refined.tif", shift_right_crop);

        Vector2f d2(0,0);
        phase_correlation_subpixel(left_image_patch, shift_right_crop,
                                  d2, subpixel_accuracy, debug);
        d += d2; // The second translation adds to the first one.
      }

      // If there is too much translation in our affine transform or we got NaNs, invalidate the pixel
      if ( norm_2(d) > SUBPIXEL_MAX_TRANSLATION ||
           std::isnan(d[0]) || std::isnan(d[1]) )
        invalidate(disparity_map(x,y));
      else
        remove_mask(disparity_map(x,y)) -= d; // TODO: Why is this subtracted?

      //vw_throw( NoImplErr() << "DEBUG!!!!" );
    } // X increment
  } // Y increment

}


}} // namespace vw::stereo

#endif // __VW_STEREO_PHASESUBPIXEL_VIEW__
