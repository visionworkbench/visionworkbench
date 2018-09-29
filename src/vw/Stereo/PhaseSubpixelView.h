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


namespace vw {
namespace stereo {

// TODO: Not needed???
/// Implement standard matrix multiplication for OpenCV complex matrices.
cv::Mat complex_multiply(cv::Mat const& a, cv::Mat const& b) {
  if (a.cols != b.rows) {
    vw_throw(ArgumentErr() << "complex_multiply dims error: a.cols = " << a.cols
                           << ", b.rows = " << b.rows);
  }
  
  typedef std::complex<float> c_type;
  cv::Mat out(a.rows, b.cols, a.type());
  for (int r=0; r<out.rows; ++r) {
    for (int c=0; c<out.cols; ++c) {
      c_type val = 0;
      for (int i=0; i<a.cols; ++i) {
        val += a.at<c_type>(r,i) * b.at<c_type>(i,c);
      }
      out.at<c_type>(r,c) = val;
    }
  }
  
  return out;
}

/// Upsample the input image and compute the DFT, but only 
  inline
cv::Mat partial_upsample_dft(cv::Mat const& input, int upsampled_height, int upsampled_width, 
                             int upscale, int row_offset=0, int col_offset=0) {
/*
function out=dftups(in,nor,noc,usfac,roff,coff)
% function out=dftups(in,nor,noc,usfac,roff,coff);
% Upsampled DFT by matrix multiplies, can compute an upsampled DFT in just
% a small region.
% usfac         Upsampling factor (default usfac = 1)
% [nor,noc]     Number of pixels in the output upsampled DFT, in
%               units of upsampled pixels (default = size(in))
% roff, coff    Row and column offsets, allow to shift the output array to
%               a region of interest on the DFT (default = 0)
% Recieves DC in upper left corner, image center must be in (1,1) 
% Manuel Guizar - Dec 13, 2007
% Modified from dftus, by J.R. Fienup 7/31/06

% This code is intended to provide the same result as if the following
% operations were performed
%   - Embed the array "in" in an array that is usfac times larger in each
%     dimension. ifftshift to bring the center of the image to (1,1).
%   - Take the FFT of the larger array
%   - Extract an [nor, noc] region of the result. Starting with the 
%     [roff+1 coff+1] element.

% It achieves this result by computing the DFT in the output array without
% the need to zeropad. Much faster and memory efficient than the
% zero-padded FFT approach if [nor noc] are much smaller than [nr*usfac nc*usfac]
*/
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
  
  std::cout << "upsampled_height = " << upsampled_height << std::endl;
  std::cout << "upsampled_width  = " << upsampled_width << std::endl;
  std::cout << "upscale  = " << upscale << std::endl;
  std::cout << "row_offset  = " << row_offset << std::endl;
  std::cout << "col_offset  = " << col_offset << std::endl;
  std::cout << "input.size() = " << input.size() << std::endl;
  std::cout << "input.rows = " << input.rows << std::endl;
  std::cout << "input.cols = " << input.cols << std::endl;
  
  cv::Mat col_kernel, row_kernel;
  //std::cout << "col_vector = " << col_vector << std::endl;
  cv::Mat temp_vec = fftshift(col_vector, true).t();
  //std::cout << "temp 1 = " << temp_vec << std::endl;
  cv::Mat v1 = temp_vec - floor(input.cols/2);
  cv::Mat v2 = col_up_vector - col_offset;
  cv::Mat m;
  std::cout << "v1.size() = " << v1.size() << std::endl;
  std::cout << "v2.size() = " << v2.size() << std::endl;
  cv::gemm(v1,v2,1.0, cv::Mat(), 0.0, m); //m  = v1*v2;
  std::cout << "m.size() = " << m.size() << std::endl;
  //std::cout << "m = \n" << m << std::endl;
  c_type constant = neg_i*two_pi/static_cast<float>(input.cols*upscale);
  std::cout << "constant = " << constant << std::endl;
  // TODO: Improve!
  for (int i=0; i<m.rows; ++i)
    for (int j=0; j<m.cols; ++j)
      m.at<c_type>(i,j) = std::exp(m.at<c_type>(i,j)*constant);
  //std::cout << "\n\nm2 = \n" << m << std::endl;
  //cv::exp( m, col_kernel ); // Does not work!
  col_kernel = m;
  std::cout << "col_kernel.size() = " << col_kernel.size() << std::endl;

  temp_vec = fftshift(row_vector, true);
  //std::cout << "temp 2 = " << temp_vec << std::endl;
  v1 = row_up_vector.t() - row_offset;
  v2 = temp_vec - floor(input.rows/2);
  std::cout << "v1.size() = " << v1.size() << std::endl;
  std::cout << "v2.size() = " << v2.size() << std::endl;
  cv::gemm(v1,v2,1.0, cv::Mat(), 0.0, m); //m  = v1*v2;
  std::cout << "m.size() = " << m.size() << std::endl;
  //std::cout << "m = " << m << std::endl;
  constant = neg_i*two_pi/static_cast<float>(input.rows*upscale);
  std::cout << "constant = " << constant << std::endl;
  std::cout << "m.type() = " << m.type() << std::endl;
  // TODO: Improve!
  for (int i=0; i<m.rows; ++i)
    for (int j=0; j<m.cols; ++j)
      m.at<c_type>(i,j) = std::exp(m.at<c_type>(i,j)*constant);
  //std::cout << "\n\nm2 = \n" << m << std::endl;
  //cv::exp( m, row_kernel ); // Does not work!
  row_kernel = m;
  std::cout << "row_kernel.size() = " << row_kernel.size() << std::endl;

  //std::cout << "row_kernel = \n" << row_kernel << std::endl;
  //std::cout << "\n\n\ncol_kernel = \n" << col_kernel << std::endl;
  //save_mag_from_ft(row_kernel,  "/home/smcmich1/data/subpixel/row_kernel.tif", false);
  //save_mag_from_ft(col_kernel,  "/home/smcmich1/data/subpixel/col_kernel.tif", false);
  
  //std::cout << "input = \n" << input << std::endl << std::endl;
  
  cv::Mat out = row_kernel*input*col_kernel;
  
  //cv::Mat temp = complex_multiply(row_kernel, input);
  //cv::Mat out  = complex_multiply(temp, col_kernel);
  
  std::cout << "out.size() = " << out.size() << std::endl;
  //std::cout << "out = \n" << out << std::endl;
  return out;
}



/// Compute the subpixel offset between two images.
/// - The images MUST be the same size!
/// - Maximum accuracy is 1/subpixel_accuracy
template <class T1, class T2>
void phase_correlation_subpixel(ImageViewBase<T1> const& left_image,
                                ImageViewBase<T2> const& right_image,
                                Vector2 &offset,
                                int subpixel_accuracy=10) {

  if (left_image.get_size() != right_image.get_size()) {
    vw_throw( ArgumentErr() 
      << "phase_correlation_subpixel requires images to be the same size!\n" );
  }
  
  // TODO: Make sure the result is not changed by padding!
  cv::Mat complexI_left, complexI_right;
  std::cout << "Getting DFT of left image...\n";
  get_dft(left_image,  complexI_left);
  std::cout << "Getting DFT of right image...\n";
  get_dft(right_image, complexI_right);


  std::cout << "Postprocessing...\n";
  save_mag_from_ft(complexI_left,  "/home/smcmich1/data/subpixel/magI_left.tif");
  save_mag_from_ft(complexI_right, "/home/smcmich1/data/subpixel/magI_right.tif");

  // Compute convolution of the two images.
  cv::Mat initial_conj;
  std::cout << "Multiply\n";
  cv::mulSpectrums(complexI_left, complexI_right, initial_conj, 0, true);
  
  save_mag_from_ft(initial_conj, "/home/smcmich1/data/subpixel/initial_conj.tif");
  
  int pad_factor = subpixel_accuracy; // Controls maximum subpixel accuracy.
  
  const int INITIAL_PAD_FACTOR = 2;
  int padded_width  = INITIAL_PAD_FACTOR*initial_conj.cols;
  int padded_height = INITIAL_PAD_FACTOR*initial_conj.rows;
  std::cout << "Padding FT\n";
  cv::Mat padded_conj = pad_fourier_transform(initial_conj, padded_width, padded_height);
  std::cout << "padded_conj.size() = " << padded_conj.size() << std::endl;
  save_mag_from_ft(padded_conj, "/home/smcmich1/data/subpixel/padded_conj.tif");
  
  // TODO: The scaling here is different than in Matlab!
  std::cout << "IFFT\n";
  cv::Mat conv;
  cv::dft(padded_conj, conv, cv::DFT_INVERSE + cv::DFT_REAL_OUTPUT+ cv::DFT_SCALE, 0); // TODO: set NonZeroRows?

  std::cout << "find peak\n";
  int width  = conv.cols;
  int height = conv.rows;
  double maxVal;
  cv::Point maxLoc;
  std::cout << "conv.type = " << conv.type() << std::endl;
  cv::minMaxLoc(conv, NULL, &maxVal, NULL, &maxLoc); // TODO: magnitude(conv)??
  std::cout << "maxLoc = " << maxLoc << std::endl;

  // Unscramble the peak location and get the final shift answer.
  float initial_shift_x = (maxLoc.x<width /2) ? (maxLoc.x) : (maxLoc.x-width );
  float initial_shift_y = (maxLoc.y<height/2) ? (maxLoc.y) : (maxLoc.y-height);
  initial_shift_x /= static_cast<float>(INITIAL_PAD_FACTOR);
  initial_shift_y /= static_cast<float>(INITIAL_PAD_FACTOR);
  std::cout << "initial_shift_x = " << initial_shift_x << std::endl;
  std::cout << "initial_shift_y = " << initial_shift_y << std::endl;

  boost::shared_ptr<cv::Mat> ocv_ptr(&conv, null_deleter);
  ImageResourceView<float> ocv_view(new ImageResourceOpenCV(ocv_ptr));
  write_image( "/home/smcmich1/data/subpixel/conv.tif", ocv_view);

  if (pad_factor <= 2) {
    offset[0] = initial_shift_x;
    offset[1] = initial_shift_y;
    return;
  }

  // Estimate refinement using a matrix multiply DFT

  // TODO: Redo the comments!

  // TODO: Make sure there is not a one pixel offset!
  
  //% Initial shift estimate in upsampled grid
  float shift_x = round(initial_shift_x*pad_factor)/pad_factor; 
  float shift_y = round(initial_shift_y*pad_factor)/pad_factor; 
  float dft_shift = floor(ceil(pad_factor*1.5)/2); //% Center of output array at dftshift+1
  //% Matrix multiply DFT around the current shift estimate
  const float UPSAMPLE_REGION_FACTOR = 1.5;
  int upsampled_height = ceil(pad_factor*UPSAMPLE_REGION_FACTOR);
  int upsampled_width  = ceil(pad_factor*UPSAMPLE_REGION_FACTOR);
  cv::Mat new_conj;
  cv::mulSpectrums(complexI_right, complexI_left, new_conj, 0, true); // Order is reversed here.
  cv::Mat partial_upsampled = partial_upsample_dft(new_conj, 
                                                   upsampled_height,
                                                   upsampled_width,
                                                   pad_factor,
                                                   dft_shift-shift_y*pad_factor,
                                                   dft_shift-shift_x*pad_factor);
  std::cout << "partial_upsampled.size() = " << partial_upsampled.size() << std::endl;
  save_mag_from_ft(partial_upsampled, "/home/smcmich1/data/subpixel/partial_upsampled.tif", false);
  //cv::Mat CC = conj(partial_upsampled); // TODO: Needed??
  cv::Mat CC;
  get_magnitude(partial_upsampled, CC);

  //% Locate maximum and map back to original pixel grid 
  std::cout << "CC.type = " << CC.type() << std::endl;
  cv::minMaxLoc(CC, NULL, &maxVal, NULL, &maxLoc); // TODO: magnitude(conv)??
  std::cout << "maxLoc = " << maxLoc << std::endl;
  //CCmax = CC(rloc,cloc); // Used for error calc
  maxLoc.y = maxLoc.y - dft_shift - 1;
  maxLoc.x = maxLoc.x - dft_shift - 1;
  shift_y = shift_y + static_cast<float>(maxLoc.y)/static_cast<float>(pad_factor);
  shift_x = shift_x + static_cast<float>(maxLoc.x)/static_cast<float>(pad_factor);
  std::cout << "shift_x = " << shift_x << std::endl;
  std::cout << "shift_y = " << shift_y << std::endl;

  offset[0] = shift_x;
  offset[1] = shift_y;

  std::cout << "Finished with subpixel calculation!\n";
}


  
  
  /*
  // TESTING: Subpixel FFT shift! --------------------------
  
  DiskImageView<PixelGray<float> > left_disk_image (left_file_name );
  DiskImageView<PixelGray<float> > right_disk_image(right_file_name );
  
  const int subpixel_accuracy = 10;
  float offset_cols, offset_rows;
  phase_correlation_subpixel(left_disk_image, right_disk_image,
                             offset_cols, offset_rows, subpixel_accuracy);
  
  std::cout << "Done with test!\n";
  return 0;
  // END TESTING ------------------------
  */
  
  



/// Subpixel view class based on the Phase Correlation method.
template <class DImageT, class Image1T, class Image2T>
class PhaseSubpixelView : public ImageViewBase<PhaseSubpixelView<DImageT,Image1T,Image2T> > {
  
/*
  /// Compute the subpixel disparity for each input integer disparity
  template <class FImage1T, class FImage2T>
  ImageView<PixelMask<Vector2f> >
  evaluate( ImageView<PixelMask<Vector2i> > const& integer_disparity, ///< Input disparity, cropped to disparity_region
            ImageViewBase<FImage1T>         const& left_filtered_image,
            ImageViewBase<FImage2T>         const& right_filtered_image,
            BBox2i const& left_region,      ///< The ROI of the left and right images that we will need
            BBox2i const& right_region,     ///  to use in order to do our computations.
            BBox2i const& disparity_region, ///< The ROI in the entire output image that we are computing
            BBox2i const& search_range ) const { ///< The range of input disparity values in disparity_region
    } // End of the evaluate() function
*/
public:
  typedef PixelMask<Vector2f> pixel_type;
  typedef pixel_type result_type;
  typedef ProceduralPixelAccessor<PhaseSubpixelView> pixel_accessor;

  PhaseSubpixelView(ImageViewBase<DImageT>    const& disparity,
                    ImageViewBase<Image1T>    const& left_image,
                    ImageViewBase<Image2T>    const& right_image,
                    PrefilterModeType prefilter_mode, float prefilter_width,
                    Vector2i const& kernel_size ) {}

  inline int32 cols  () const { return 0; }
  inline int32 rows  () const { return 0; }
  inline int32 planes() const { return 1; }

  inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }
  inline pixel_type operator() ( int32 /*i*/, int32 /*j*/, int32 /*p*/ = 0 ) const {
    vw_throw( NoImplErr() << "SubpixelView:operator() has not been implemented." );
    return pixel_type();
  }

  // Block rasterization section does the actual work
  typedef CropView<ImageView<pixel_type> > prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const {

  } // end function prerasterize

  template <class DestT>
  inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
    vw::rasterize(prerasterize(bbox), dest, bbox );
  }
}; // End class ParabolaSubpixelView

template <class DImageT, class Image1T, class Image2T>
PhaseSubpixelView<DImageT, Image1T, Image2T>
phase_subpixel( ImageViewBase<DImageT>    const& disparity,
                ImageViewBase<Image1T>    const& left_image,
                ImageViewBase<Image2T>    const& right_image,
                PrefilterModeType prefilter_mode, float prefilter_width,
                Vector2i const& kernel_size ) {
  typedef PhaseSubpixelView<DImageT, Image1T, Image2T> result_type;
  return result_type( disparity.impl(), left_image.impl(), right_image.impl(),
                      prefilter_mode, prefilter_width, kernel_size );
}

}} // namespace vw::stereo

#endif // __VW_STEREO_PHASESUBPIXEL_VIEW__
