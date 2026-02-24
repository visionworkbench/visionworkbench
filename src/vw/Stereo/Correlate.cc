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

#include <vw/Image/ImageMath.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/PixelAccessors.h>
#include <vw/Image/Filter.h>
#include <vw/Stereo/Correlate.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/LinearAlgebra.h>
#include <vw/Image/Interpolation.h>
#include <vw/Image/Manipulation.h>

using namespace vw;
using namespace vw::stereo;

// TODO(oalexan1): Many of these are unused
namespace vw { namespace stereo { namespace detail {

  ImageView<float> compute_spatial_weight_image(int32 kern_width, int32 kern_height,
                                                float two_sigma_sqr) {
    int32 center_pix_x = kern_width/2;
    int32 center_pix_y = kern_height/2;
    float sum;
    
    sum = 0.0;
    ImageView<float> weight(kern_width, kern_height);
    for (int32 j = 0; j < kern_height; ++j) {
      for (int32 i = 0; i < kern_width; ++i ) {
        weight(i,j) = exp(-1*((i-center_pix_x)*(i-center_pix_x) +
                              (j-center_pix_y)*(j-center_pix_y)) / two_sigma_sqr);
        sum += weight(i,j);
      }
    }
    
    weight /= sum;
    
    return weight;
  }
  
  inline double huber_robust_coefficient (double delta_norm, double b) {
    if (delta_norm < b)
      return delta_norm*delta_norm;
    else
      return 2*b*delta_norm - b*b;
  }
  
  inline float
  cauchy_robust_coefficient (float delta_norm, float b) {
    float b_sqr = b*b;
    return b_sqr*logf(1+delta_norm*delta_norm/b_sqr);
  }
  
  inline double
  blake_zisserman_robust_coefficient (double delta_norm, double b) {
    return -log(exp(-(delta_norm*delta_norm) ) + b);
  }
  
  struct HuberError {
    double m_b;
    HuberError(double b) : m_b(b) {}
    
    double operator() (double delta_norm) {
      if (delta_norm < m_b)
        return delta_norm*delta_norm;
      else
        return 2*m_b*delta_norm - m_b*m_b;
    }
  };
  
  inline ImageView<float>
  compute_gaussian_weight_image(int kern_width, int kern_height) {
    int center_pix_x = kern_width/2;
    int center_pix_y = kern_height/2;
    float two_sigma_sqr = 2.0*pow(float(kern_width)/7.0,2.0);
    
    ImageView<float> weight(kern_width, kern_height);
    for (int j = 0; j < kern_height; ++j) {
      for (int i = 0; i < kern_width; ++i ) {
        weight(i,j) = exp(-1*((i-center_pix_x)*(i-center_pix_x) +
                              (j-center_pix_y)*(j-center_pix_y)) / two_sigma_sqr);
      }
    }
    return weight;
  }

}}} // End namespace vw::stereo::detail

/// Research-y implementation of affine EM correlation.
void vw::stereo::subpixel_correlation_affine_2d_EM(ImageView<PixelMask<Vector2f> > &disparity_map,
                                                   ImageView<float> const& left_image,
                                                   ImageView<float> const& right_image,
                                                   int32 kern_width, int32 kern_height,
                                                   BBox2i region_of_interest,
                                                   bool do_horizontal_subpixel,
                                                   bool do_vertical_subpixel,
                                                   bool /*verbose*/ ) {
  typedef Vector<float,6> Vector6f;
  typedef Matrix<float,6,6> Matrix6x6f;
  typedef typename ImageView<float>::pixel_accessor ImageViewFAcc;
  typedef typename CropView<ImageView<float> >::pixel_accessor CropViewFAcc;
  typedef typename CropView<ImageView<float> >::pixel_accessor CropViewTAcc;

  // Bail out if no subpixel computation has been requested
  if (!do_horizontal_subpixel && !do_vertical_subpixel) return;

  // Fixed consts
  const unsigned M_MAX_EM_ITER = 2;
  //const float M_MIN_VAR2_PLANE = 1e-6;
  //const float M_MIN_VAR2_NOISE = 1e-6;
  const float two_sigma_sqr = 2.0*pow(float(kern_width)/5.0,2.0);

  VW_ASSERT( disparity_map.cols() == left_image.cols() &&
             disparity_map.rows() == left_image.rows(),
             ArgumentErr() << "subpixel_correlation: left image and "
             << "disparity map do not have the same dimensions.");

  // Currently unused confidence_image code
  //ImageView<float> confidence_image(left_image.cols(),
  //                                  left_image.rows());

  // Interpolated Input Images
  InterpolationView<EdgeExtensionView<ImageView<float>, ZeroEdgeExtension>, BilinearInterpolation> right_interp_image =
    interpolate(right_image, BilinearInterpolation(), ZeroEdgeExtension());

  // This is the maximum number of pixels that the solution can be
  // adjusted by affine subpixel refinement.
  float AFFINE_SUBPIXEL_MAX_TRANSLATION = kern_width/2;

  int32 kern_half_height = kern_height/2;
  int32 kern_half_width = kern_width/2;
  int32 kern_pixels = kern_height * kern_width;
  int32 weight_threshold = kern_pixels/2;

  ImageView<float> x_deriv = derivative_filter(left_image, 1, 0);
  ImageView<float> y_deriv = derivative_filter(left_image, 0, 1);
  ImageView<float> weight_template =
    vw::stereo::detail::compute_spatial_weight_image(kern_width, kern_height, two_sigma_sqr);

  // Workspace images are allocated up here out of the tight inner
  // loop.  We rasterize into these directly in the code below.
  ImageView<float> w(kern_width, kern_height);

  // Iterate over all of the pixels in the disparity map except for
  // the outer edges.
  for ( int32 y = std::max(region_of_interest.min().y()-1,kern_half_height);
        y < std::min(left_image.rows()-kern_half_height,
                     region_of_interest.max().y()+1); ++y) {

    for (int32 x=std::max(region_of_interest.min().x()-1,kern_half_width);
         x < std::min(left_image.cols()-kern_half_width,
                      region_of_interest.max().x()+1); ++x) {

      BBox2i current_window(x-kern_half_width, y-kern_half_height,
                            kern_width, kern_height);

      // Skip over pixels for which we have no initial disparity estimate
      if ( !is_valid(disparity_map(x,y)) )
        continue;

      // Define and initialize the model params
      // Initialize our affine transform with the identity.  The
      // entries of d are laid out in row major order:
      //
      //   | d(0) d(1) d(2) |
      //   | d(3) d(4) d(5) |
      //   |  0    0    1   |
      //
      Vector6f d;
      d(0) = 1.0; d(1) = 0.0; d(2) = 0.0;
      d(3) = 0.0; d(4) = 1.0; d(5) = 0.0;

      // Compute the derivative image patches
      CropView<ImageView<float> > left_image_patch = crop(left_image, current_window);
      CropView<ImageView<float> > I_x = crop(x_deriv, current_window);
      CropView<ImageView<float> > I_y = crop(y_deriv, current_window);

      // Compute the base weight image
      int32 good_pixels =
        adjust_weight_image(w, crop(disparity_map, current_window),
                            weight_template);

      // Skip over pixels for which there are very few good matches
      // in the neighborhood.
      if (good_pixels < weight_threshold) {
        invalidate(disparity_map(x,y));
        continue;
      }

      float curr_sum_I_e_val = 0.0;
      float prev_sum_I_e_val = 0.0;

      // Iterate until a solution is found or the max number of
      // iterations is reached.
      for (unsigned iter = 0; iter < 10; ++iter) {
        // First we check to see if our current subpixel translation
        // is less than one half of the window width.  If not, then
        // we are probably having trouble converging and we abort
        // this pixel!!
        if (norm_2( Vector2f(d[2],d[5]) ) > AFFINE_SUBPIXEL_MAX_TRANSLATION)
          break;

        float x_base = x + disparity_map(x,y)[0];
        float y_base = y + disparity_map(x,y)[1];

        Matrix6x6f rhs;
        Vector6f lhs, prev_lhs;

        //initial values here.
        ImageView<float> gamma_plane(kern_width, kern_height);
        ImageView<float> gamma_noise(kern_width, kern_height);

        //set init params - START
        float var2_plane = 1e-3;
        float mean_noise = 0.0;
        float var2_noise = 1e-2;
        float w_plane    = 0.8;
        float w_noise    = 0.2;
        //set init params - END

        float in_curr_sum_I_e_val = 0.0;
        float in_prev_sum_I_e_val = 1000000.0;
        Vector6f d_em;
        d_em = d;

        for (unsigned em_iter=0; em_iter < M_MAX_EM_ITER; em_iter++){
          float noise_norm_factor = 1.0/sqrt(2*M_PI*var2_noise);
          float plane_norm_factor = 1.0/sqrt(2*M_PI*var2_plane);

          //reset lhs and rhs
          std::fill( lhs.begin(), lhs.end(), 0.0f );
          std::fill( rhs.begin(), rhs.end(), 0.0f );

          in_curr_sum_I_e_val = 0.0;
          float mean_noise_tmp  = 0.0;
          float sum_gamma_noise = 0.0;
          float sum_gamma_plane = 0.0;

          // Set up pixel accessors
          CropViewFAcc I_x_row = I_x.origin(), I_y_row = I_y.origin();
          CropViewTAcc left_image_patch_row = left_image_patch.origin();
          ImageViewFAcc w_row = w.origin();

          // Perform loop that does Expectation and Maximization in one go
          for (int32 jj = -kern_half_height; jj <= kern_half_height; ++jj) {
            ImageViewFAcc w_ptr = w_row;
            CropViewFAcc I_x_ptr = I_x_row, I_y_ptr = I_y_row;
            CropViewTAcc left_image_patch_ptr = left_image_patch_row;
            int32 gamma_iy = jj + kern_half_height;
            float xx_partial = x_base + d[1] * jj + d[2];
            float yy_partial = y_base + d[4] * jj + d[5];
            float delta_x_partial = d_em[1] * jj + d_em[2];
            float delta_y_partial = d_em[4] * jj + d_em[5];

            for (int32 ii = -kern_half_width; ii <= kern_half_width; ++ii) {
              // First we compute the pixel offset for the right image
              // and the error for the current pixel.
              int32 gamma_ix = ii + kern_half_width;
              float xx = d[0] * ii + xx_partial;
              float yy = d[3] * ii + yy_partial;
              float delta_x = d_em[0] * ii + delta_x_partial;
              float delta_y = d_em[3] * ii + delta_y_partial;

              /// Expectation
              float interpreted_px = right_interp_image(xx,yy);
              float I_e_val = interpreted_px - (*left_image_patch_ptr);
              in_curr_sum_I_e_val += I_e_val;
              float temp_plane = I_e_val - delta_x*(*I_x_ptr) -
                delta_y*(*I_y_ptr);
              float temp_noise = interpreted_px - mean_noise;
              float plane_prob_exp = // precompute to avoid underflow
                -1*(temp_plane*temp_plane)/(2*var2_plane);
              float plane_prob =
                (plane_prob_exp < -75) ? 0.0f : plane_norm_factor *
                exp(plane_prob_exp);
              float noise_prob_exp =
                -1*(temp_noise*temp_noise)/(2*var2_noise);
              float noise_prob =
                (noise_prob_exp < -75) ? 0.0f : noise_norm_factor *
                exp(noise_prob_exp);

              float sum = plane_prob*w_plane + noise_prob*w_noise;
              gamma_plane(gamma_ix, gamma_iy) = plane_prob*w_plane/sum;
              gamma_noise(gamma_ix, gamma_iy) = noise_prob*w_noise/sum;
              // End Expectation

              // Maximization (computing the d_em vector)
              mean_noise_tmp +=
                interpreted_px * gamma_noise(gamma_ix, gamma_iy);
              sum_gamma_plane += gamma_plane(gamma_ix, gamma_iy);
              sum_gamma_noise += gamma_noise(gamma_ix, gamma_iy);

              float robust_weight = gamma_plane(gamma_ix, gamma_iy);

              // We combine the error value with the derivative and
              // add this to the update equation.
              float weight  = robust_weight*(*w_ptr);
              if ( weight < 1e-26 ) {
                // avoid underflow
                I_x_ptr.next_col();
                I_y_ptr.next_col();
                left_image_patch_ptr.next_col();
                continue;
              }
              float I_x_val = weight * (*I_x_ptr);
              float I_y_val = weight * (*I_y_ptr);
              float I_x_sqr = I_x_val * (*I_x_ptr);
              float I_y_sqr = I_y_val * (*I_y_ptr);
              float I_x_I_y = I_x_val * (*I_y_ptr);

              // Left hand side
              lhs(0) += ii * I_x_val * I_e_val;
              lhs(1) += jj * I_x_val * I_e_val;
              lhs(2) +=      I_x_val * I_e_val;
              lhs(3) += ii * I_y_val * I_e_val;
              lhs(4) += jj * I_y_val * I_e_val;
              lhs(5) +=      I_y_val * I_e_val;

              // Right Hand Side UL
              rhs(0,0) += ii*ii * I_x_sqr;
              rhs(0,1) += ii*jj * I_x_sqr;
              rhs(0,2) += ii    * I_x_sqr;
              rhs(1,1) += jj*jj * I_x_sqr;
              rhs(1,2) += jj    * I_x_sqr;
              rhs(2,2) +=         I_x_sqr;

              // Right Hand Side UR
              rhs(0,3) += ii*ii * I_x_I_y;
              rhs(0,4) += ii*jj * I_x_I_y;
              rhs(0,5) += ii    * I_x_I_y;
              rhs(1,4) += jj*jj * I_x_I_y;
              rhs(1,5) += jj    * I_x_I_y;
              rhs(2,5) +=         I_x_I_y;

              // Right Hand Side LR
              rhs(3,3) += ii*ii * I_y_sqr;
              rhs(3,4) += ii*jj * I_y_sqr;
              rhs(3,5) += ii    * I_y_sqr;
              rhs(4,4) += jj*jj * I_y_sqr;
              rhs(4,5) += jj    * I_y_sqr;
              rhs(5,5) +=         I_y_sqr;
              // End Maximization

              I_x_ptr.next_col();
              I_y_ptr.next_col();
              left_image_patch_ptr.next_col();
            }
            I_x_row.next_row();
            I_y_row.next_row();
            left_image_patch_row.next_row();
          }
          lhs *= -1;

          // Fill in symmetric entries
          rhs(1,0) = rhs(0,1);
          rhs(2,0) = rhs(0,2);
          rhs(2,1) = rhs(1,2);
          rhs(1,3) = rhs(0,4);
          rhs(2,3) = rhs(0,5);
          rhs(2,4) = rhs(1,5);
          rhs(3,0) = rhs(0,3);
          rhs(3,1) = rhs(1,3);
          rhs(3,2) = rhs(2,3);
          rhs(4,0) = rhs(0,4);
          rhs(4,1) = rhs(1,4);
          rhs(4,2) = rhs(2,4);
          rhs(4,3) = rhs(3,4);
          rhs(5,0) = rhs(0,5);
          rhs(5,1) = rhs(1,5);
          rhs(5,2) = rhs(2,5);
          rhs(5,3) = rhs(3,5);
          rhs(5,4) = rhs(4,5);

          // Solves lhs = rhs * x, and stores the result in-place in lhs.
          try {
            solve_symmetric_modify(rhs,lhs);
          } catch (const ArgumentErr& /*e*/) {} // Do Nothing

          //normalize the mean of the noise
          mean_noise = mean_noise_tmp/sum_gamma_noise;

          /* CURRENTLY UNUSED VARIANCE CALCULATION CODE
          //compute the variance for noise and plane
          float var2_noise_tmp  = 0.0;
          float var2_plane_tmp  = 0.0;

          // Set up pixel accessors
          I_x_row = I_x.origin();
          I_y_row = I_y.origin();
          left_image_patch_row = left_image_patch.origin();

          for (int32 jj = -kern_half_height; jj <= kern_half_height; ++jj) {
          CropViewFAcc I_x_ptr = I_x_row, I_y_ptr = I_y_row;
          CropViewTAcc left_image_patch_ptr = left_image_patch_row;
          int32 gamma_iy = jj + kern_half_height;
          float xx_partial = x_base + d[1] * jj + d[2];
          float yy_partial = y_base + d[4] * jj + d[5];
          float delta_x_partial = d_em[1] * jj + d_em[2];
          float delta_y_partial = d_em[4] * jj + d_em[5];

          for (int32 ii = -kern_half_width; ii <= kern_half_width; ++ii) {
          // First we compute the pixel offset for the right image
          // and the error for the current pixel.
          int32 gamma_ix = ii + kern_half_width;
          float xx = d[0] * ii + xx_partial;
          float yy = d[3] * ii + yy_partial;
          float delta_x = d_em[0] * ii + delta_x_partial;
          float delta_y = d_em[3] * ii + delta_y_partial;

          float interpreted_px = right_interp_image(xx,yy);
          float I_e_val = interpreted_px - (*left_image_patch_ptr);
          float temp_plane = I_e_val - delta_x*(*I_x_ptr) -
          delta_y*(*I_y_ptr);
          float temp_noise = interpreted_px - mean_noise;

          var2_noise_tmp += temp_noise*temp_noise*
          gamma_noise(gamma_ix, gamma_iy);
          var2_plane_tmp += temp_plane*temp_plane*
          gamma_plane(gamma_ix, gamma_iy);

          I_x_ptr.next_col();
          I_y_ptr.next_col();
          left_image_patch_ptr.next_col();
          }
          I_x_row.next_row();
          I_y_row.next_row();
          left_image_patch_row.next_row();
          }

          confidence_image(x,y) = in_curr_sum_I_e_val;
          var2_noise = var2_noise_tmp/sum_gamma_noise;
          var2_plane = var2_plane_tmp/sum_gamma_plane;

          if (var2_noise < M_MIN_VAR2_NOISE) var2_noise = M_MIN_VAR2_NOISE;
          if (var2_plane < M_MIN_VAR2_PLANE) var2_plane = M_MIN_VAR2_PLANE;
          */

          w_plane = sum_gamma_plane/(float)(kern_pixels);
          w_noise = sum_gamma_noise/(float)(kern_pixels);

          //Termination
          float conv_error = norm_2(prev_lhs - lhs);
          d_em = d + lhs;
          if (in_curr_sum_I_e_val < 0)
            in_curr_sum_I_e_val = - in_curr_sum_I_e_val;

          curr_sum_I_e_val = in_curr_sum_I_e_val;
          prev_lhs = lhs;

          // Termination condition
          if ((conv_error < 1E-3) && (em_iter > 0))
            break;
          else
            in_prev_sum_I_e_val = in_curr_sum_I_e_val;

        } // for em_iter end

        d += lhs;
        if (curr_sum_I_e_val < 0)
          curr_sum_I_e_val = - curr_sum_I_e_val;

        // Termination condition
        if ((prev_sum_I_e_val < curr_sum_I_e_val) && (iter > 0))
          break;
        else
          prev_sum_I_e_val = curr_sum_I_e_val;

      }

      if ( norm_2( Vector2f(d[2],d[5]) ) >
           AFFINE_SUBPIXEL_MAX_TRANSLATION ||
           std::isnan(d[2]) || std::isnan(d[5]) )
        invalidate(disparity_map(x,y));
      else
        remove_mask(disparity_map(x,y)) += Vector2f(d[2],d[5]);
    } // X increment
  } // Y increment
}

// Speed at all cost implementation
//
// In this version we don't keep around future research ideas
// since they are slow.
void vw::stereo::subpixel_optimized_affine_2d_EM(ImageView<PixelMask<Vector2f>> & disparity_map,
                                                 ImageView<float> const& left_image,
                                                 ImageView<float> const& right_image,
                                                 int32  kern_width, int32 kern_height,
                                                 BBox2i region_of_interest,
                                                 bool   do_horizontal_subpixel,
                                                 bool   do_vertical_subpixel,
                                                 bool /*verbose*/ ) {
  typedef Vector<float,6>   Vector6f;
  typedef Matrix<float,6,6> Matrix6x6f;
  typedef typename ImageView<float>::pixel_accessor ImageViewFAcc;
  typedef typename CropView<ImageView<float>>::pixel_accessor CropViewFAcc;
  typedef typename CropView<ImageView<float>>::pixel_accessor CropViewTAcc;

  // Bail out if no subpixel computation has been requested
  if (!do_horizontal_subpixel && !do_vertical_subpixel) return;

  // Fixed consts
  const unsigned M_MAX_EM_ITER = 2;
  const float two_sigma_sqr = 2.0*pow(float(kern_width)/5.0,2.0);

  VW_ASSERT( disparity_map.cols() == left_image.cols() &&
             disparity_map.rows() == left_image.rows(),
             ArgumentErr() << "subpixel_correlation: left image and "
             << "disparity map do not have the same dimensions.");

  // Interpolated Input Images
  InterpolationView<EdgeExtensionView<ImageView<float>, ZeroEdgeExtension>, BilinearInterpolation> right_interp_image =
    interpolate(right_image, BilinearInterpolation(), ZeroEdgeExtension());

  // This is the maximum number of pixels that the solution can be
  // adjusted by affine subpixel refinement.
  float AFFINE_SUBPIXEL_MAX_TRANSLATION = kern_width/2;

  const int32 kern_half_height = kern_height/2;
  const int32 kern_half_width  = kern_width/2;
  const int32 kern_pixels      = kern_height * kern_width;
  const int32 weight_threshold = kern_pixels/2;

  ImageView<float> x_deriv = derivative_filter(left_image, 1, 0).impl();
  ImageView<float> y_deriv = derivative_filter(left_image, 0, 1).impl();
  ImageView<float> weight_template =
    vw::stereo::detail::compute_spatial_weight_image(kern_width, kern_height, two_sigma_sqr);

  // Workspace images are allocated up here out of the tight inner
  // loop.  We rasterize into these directly in the code below.
  ImageView<float> w(kern_width, kern_height);

  // Iterate over all of the pixels in the disparity map except for
  // the outer edges.
  for ( int32 y = std::max(region_of_interest.min().y()-1,kern_half_height);
        y < std::min(left_image.rows()-kern_half_height,
                     region_of_interest.max().y()+1); ++y) {

    for (int32 x=std::max(region_of_interest.min().x()-1,kern_half_width);
         x < std::min(left_image.cols()-kern_half_width,
                      region_of_interest.max().x()+1); ++x) {

      BBox2i current_window(x-kern_half_width, y-kern_half_height,
                            kern_width, kern_height);

      // Skip over pixels for which we have no initial disparity estimate
      if ( !is_valid(disparity_map(x,y)) )
        continue;

      // Define and initialize the model params
      // Initialize our affine transform with the identity.  The
      // entries of d are laid out in row major order:
      //
      //   | d(0) d(1) d(2) |
      //   | d(3) d(4) d(5) |
      //   |  0    0    1   |
      //
      Vector6f d;
      d(0) = 1.0; d(1) = 0.0; d(2) = 0.0;
      d(3) = 0.0; d(4) = 1.0; d(5) = 0.0;
      float *dPtr = &(d[0]);

      // Compute the derivative image patches
      CropView<ImageView<float> > left_image_patch = crop(left_image, current_window);
      CropView<ImageView<float> > I_x = crop(x_deriv, current_window);
      CropView<ImageView<float> > I_y = crop(y_deriv, current_window);

      // Compute the base weight image
      int32 good_pixels = adjust_weight_image(w, crop(disparity_map, current_window), weight_template);

      // Skip over pixels for which there are very few good matches
      // in the neighborhood.
      if (good_pixels < weight_threshold) {
        invalidate(disparity_map(x,y));
        continue;
      }

      float curr_sum_I_e_val = 0.0;
      float prev_sum_I_e_val = 0.0;

      // Iterate until a solution is found or the max number of
      // iterations is reached.
      for (unsigned iter = 0; iter < 10; ++iter) {
        // First we check to see if our current subpixel translation
        // is less than one half of the window width.  If not, then
        // we are probably having trouble converging and we abort
        // this pixel!!
        if (norm_2( Vector2f(d[2],d[5]) ) > AFFINE_SUBPIXEL_MAX_TRANSLATION)
          break;

        float x_base = x + disparity_map(x,y)[0];
        float y_base = y + disparity_map(x,y)[1];

        Matrix6x6f rhs;
        Vector6f lhs, prev_lhs;

        //set init params - START
        float var2_plane = 1e-3;
        float mean_noise = 0.0;
        float var2_noise = 1e-2;
        float w_plane    = 0.8;
        float w_noise    = 0.2;
        //set init params - END

        float in_curr_sum_I_e_val = 0.0;
        Vector6f d_em;
        d_em = d;
        float *d_emPtr = &(d_em[0]);

        for (unsigned em_iter=0; em_iter < M_MAX_EM_ITER; em_iter++){
          float noise_norm_factor = 1.0/sqrt(2*M_PI*var2_noise);
          float plane_norm_factor = 1.0/sqrt(2*M_PI*var2_plane);

          //reset lhs and rhs
          std::fill( lhs.begin(), lhs.end(), 0.0f );
          std::fill( rhs.begin(), rhs.end(), 0.0f );

          in_curr_sum_I_e_val = 0.0;
          float mean_noise_tmp  = 0.0;
          float sum_gamma_noise = 0.0;
          float sum_gamma_plane = 0.0;

          // Set up pixel accessors
          CropViewFAcc I_x_row = I_x.origin(), I_y_row = I_y.origin();
          CropViewTAcc left_image_patch_row = left_image_patch.origin();
          ImageViewFAcc w_row = w.origin();

          // Counts pixels skipped on evaluation.
          int32 skip = 0;

          // Perform loop that does Expectation and Maximization in one go
          for (int32 jj = -kern_half_height; jj <= kern_half_height; ++jj) {
            ImageViewFAcc w_ptr = w_row;
            CropViewFAcc I_x_ptr = I_x_row, I_y_ptr = I_y_row;
            CropViewTAcc left_image_patch_ptr = left_image_patch_row;
            float xx_partial      = x_base + d[1] * jj + d[2];
            float yy_partial      = y_base + d[4] * jj + d[5];
            float delta_x_partial = d_em[1] * jj + d_em[2];
            float delta_y_partial = d_em[4] * jj + d_em[5];

            for (int32 ii = -kern_half_width; ii <= kern_half_width; ++ii) {
              // First we compute the pixel offset for the right image
              // and the error for the current pixel.
              float xx      = dPtr   [0] * ii + xx_partial;
              float yy      = dPtr   [3] * ii + yy_partial;
              float delta_x = d_emPtr[0] * ii + delta_x_partial;
              float delta_y = d_emPtr[3] * ii + delta_y_partial;

              /// Expectation
              float interpreted_px = right_interp_image(xx,yy);  // 10% of function time is spent in this call, 6% performing bounds checking.
              float I_e_val = interpreted_px - (*left_image_patch_ptr);
              in_curr_sum_I_e_val += I_e_val;
              float temp_plane     = I_e_val - delta_x*(*I_x_ptr) - delta_y*(*I_y_ptr);
              float temp_noise     = interpreted_px - mean_noise;
              float plane_prob_exp = -1*(temp_plane*temp_plane)/(2*var2_plane); // precompute to avoid underflow
              float plane_prob     = (plane_prob_exp < -75) ? 0.0f : plane_norm_factor * exp(plane_prob_exp);
              float noise_prob_exp = -1*(temp_noise*temp_noise)/(2*var2_noise);
              float noise_prob     = (noise_prob_exp < -75) ? 0.0f : noise_norm_factor * exp(noise_prob_exp);
              // 40% of function time is spent on the two exp() calls

              float sum         = plane_prob*w_plane + noise_prob*w_noise;
              float gamma_plane = plane_prob*w_plane/sum;
              float gamma_noise = noise_prob*w_noise/sum;
              // End Expectation

              // Maximization (computing the d_em vector)
              mean_noise_tmp  += interpreted_px * gamma_noise;
              sum_gamma_plane += gamma_plane;
              sum_gamma_noise += gamma_noise;

              // We combine the error value with the derivative and
              // add this to the update equation.
              float weight  = gamma_plane*(*w_ptr);
              if ( weight < 1e-26 ) {
                // avoid underflow
                I_x_ptr.next_col();
                I_y_ptr.next_col();
                left_image_patch_ptr.next_col();
                skip++;
                continue;
              }
              float I_x_val = weight  * (*I_x_ptr);
              float I_y_val = weight  * (*I_y_ptr);
              float I_x_sqr = I_x_val * (*I_x_ptr);
              float I_y_sqr = I_y_val * (*I_y_ptr);
              float I_x_I_y = I_x_val * (*I_y_ptr);

              // Left hand side
              lhs(0) -= ii * I_x_val * I_e_val;
              lhs(1) -= jj * I_x_val * I_e_val;
              lhs(2) -=      I_x_val * I_e_val;
              lhs(3) -= ii * I_y_val * I_e_val;
              lhs(4) -= jj * I_y_val * I_e_val;
              lhs(5) -=      I_y_val * I_e_val;

              float multipliers[3];
              multipliers[0] = ii*ii;
              multipliers[1] = ii*jj;
              multipliers[2] = jj*jj;



              // Right Hand Side UL
              float* rhsData = rhs.data(); // TEST: Access RHS by raw data pointer
              rhsData[ 0] += multipliers[0] * I_x_sqr;
              rhsData[ 1] += multipliers[1] * I_x_sqr;
              rhsData[ 2] += ii    * I_x_sqr;
              rhsData[ 7] += multipliers[2] * I_x_sqr;
              rhsData[ 8] += jj    * I_x_sqr;
              rhsData[14] +=         I_x_sqr;

              // Right Hand Side UR
              rhsData[ 3] += multipliers[0] * I_x_I_y;
              rhsData[ 4] += multipliers[1] * I_x_I_y;
              rhsData[ 5] += ii    * I_x_I_y;
              rhsData[10] += multipliers[2] * I_x_I_y;
              rhsData[11] += jj    * I_x_I_y;
              rhsData[17] +=         I_x_I_y;

              // Right Hand Side LR
              rhsData[21] += multipliers[0] * I_y_sqr;
              rhsData[22] += multipliers[1] * I_y_sqr;
              rhsData[23] += ii    * I_y_sqr;
              rhsData[28] += multipliers[2] * I_y_sqr;
              rhsData[29] += jj    * I_y_sqr;
              rhsData[35] +=         I_y_sqr;
              
              
/*
              // Right Hand Side UL
              rhs(0,0) += multipliers[0] * I_x_sqr;
              rhs(0,1) += multipliers[1] * I_x_sqr;
              rhs(0,2) += ii    * I_x_sqr;
              rhs(1,1) += multipliers[2] * I_x_sqr;
              rhs(1,2) += jj    * I_x_sqr;
              rhs(2,2) +=         I_x_sqr;

              // Right Hand Side UR
              rhs(0,3) += multipliers[0] * I_x_I_y;
              rhs(0,4) += multipliers[1] * I_x_I_y;
              rhs(0,5) += ii    * I_x_I_y;
              rhs(1,4) += multipliers[2] * I_x_I_y;
              rhs(1,5) += jj    * I_x_I_y;
              rhs(2,5) +=         I_x_I_y;

              // Right Hand Side LR
              rhs(3,3) += multipliers[0] * I_y_sqr;
              rhs(3,4) += multipliers[1] * I_y_sqr;
              rhs(3,5) += ii    * I_y_sqr;
              rhs(4,4) += multipliers[2] * I_y_sqr;
              rhs(4,5) += jj    * I_y_sqr;
              rhs(5,5) +=         I_y_sqr;
              // End Maximization
*/
              I_x_ptr.next_col();
              I_y_ptr.next_col();
              left_image_patch_ptr.next_col();
            }
            I_x_row.next_row();
            I_y_row.next_row();
            left_image_patch_row.next_row();
          }

          // Checking for early termination
          if ( skip == kern_pixels )
            break;

          // Fill in symmetric entries
          rhs(1,0) = rhs(0,1);
          rhs(2,0) = rhs(0,2);
          rhs(2,1) = rhs(1,2);
          rhs(3,0) = rhs(0,3);
          rhs(1,3) = rhs(3,1) = rhs(4,0) = rhs(0,4);
          rhs(2,3) = rhs(3,2) = rhs(5,0) = rhs(0,5);
          rhs(4,1) = rhs(1,4);
          rhs(2,4) = rhs(4,2) = rhs(5,1) = rhs(1,5);
          rhs(5,2) = rhs(2,5);
          rhs(4,3) = rhs(3,4);
          rhs(5,3) = rhs(3,5);
          rhs(5,4) = rhs(4,5);

          // Solves lhs = rhs * x, and stores the result in-place in lhs.
          const math::f77_int n = 6;
          const math::f77_int nrhs = 1;
          math::f77_int info;
          math::posv('L',n,nrhs,&(rhs(0,0)), n, &(lhs(0)), n, &info);

          //normalize the mean of the noise
          mean_noise = mean_noise_tmp/sum_gamma_noise;

          w_plane = sum_gamma_plane/(float)(kern_pixels);
          w_noise = sum_gamma_noise/(float)(kern_pixels);

          //Termination
          float conv_error = norm_2(prev_lhs - lhs);
          d_em = d + lhs;
          if (in_curr_sum_I_e_val < 0)
            in_curr_sum_I_e_val = - in_curr_sum_I_e_val;

          curr_sum_I_e_val = in_curr_sum_I_e_val;
          prev_lhs = lhs;

          // Termination condition
          if ((conv_error < 1E-3) && (em_iter > 0))
            break;

        } // for em_iter end

        d += lhs;
        if (curr_sum_I_e_val < 0)
          curr_sum_I_e_val = - curr_sum_I_e_val;

        // Termination condition
        if ((prev_sum_I_e_val < curr_sum_I_e_val) && (iter > 0))
          break;
        else
          prev_sum_I_e_val = curr_sum_I_e_val;

      }

      if ( norm_2( Vector2f(d[2],d[5]) ) >
           AFFINE_SUBPIXEL_MAX_TRANSLATION ||
           std::isnan(d[2]) || std::isnan(d[5]) )
        invalidate(disparity_map(x,y));
      else
        remove_mask(disparity_map(x,y)) += Vector2f(d[2],d[5]);
    } // X increment
  } // Y increment
}

/// Affine subpixel correlation function
/// - Similar to the _em function but about 5X faster and less accurate.
void vw::stereo::subpixel_optimized_affine_2d(ImageView<PixelMask<Vector2f> > &disparity_map,
                                              ImageView<float> const& left_image,
                                              ImageView<float> const& right_image,
                                              int32  kern_width, int32 kern_height,
                                              BBox2i region_of_interest,
                                              bool   do_horizontal_subpixel,
                                              bool   do_vertical_subpixel,
                                              bool /*verbose*/ ) {
  
  typedef Vector<float,6  > Vector6f;
  typedef Matrix<float,6,6> Matrix6x6f;
  typedef typename ImageView<float>::pixel_accessor ImageViewFAcc;
  typedef typename CropView<ImageView<float   > >::pixel_accessor CropViewFAcc;
  typedef typename CropView<ImageView<float> >::pixel_accessor CropViewTAcc;

  // Bail out if no subpixel computation has been requested
  if (!do_horizontal_subpixel && !do_vertical_subpixel) return;

  // Fixed consts
  const unsigned MAX_NUM_ITERATIONS = 10;  // <-- Does not seem to make a big difference around 8ish
  const float    two_sigma_sqr = 2.0*pow(float(kern_width)/5.0,2.0);

  VW_ASSERT( disparity_map.cols() == left_image.cols() &&
             disparity_map.rows() == left_image.rows(),
             ArgumentErr() << "subpixel_correlation: left image and "
             << "disparity map do not have the same dimensions.");

  // It seems like Bicubic interpolation should help here, but when tested it produced undetectable
  //  changes in the results and significantly increased the run time.  This is true for both the pyramid
  //  and non-pyramid version of this algorithm.

  // Interpolated Input Images
  InterpolationView<EdgeExtensionView<ImageView<float>, ZeroEdgeExtension>, BilinearInterpolation> right_interp_image =
         interpolate(right_image, BilinearInterpolation(), ZeroEdgeExtension());
  // An alternate interpolation view with no bounds checking!
  InterpolationView<EdgeExtensionView<ImageView<float>, NoEdgeExtension>, BilinearInterpolation> right_interp_image_unsafe =
         interpolate(right_image, BilinearInterpolation(), NoEdgeExtension());



  // This is the maximum number of pixels that the solution can be
  // adjusted by affine subpixel refinement.
  float AFFINE_SUBPIXEL_MAX_TRANSLATION = kern_width/2;

  const int32 kern_half_height    = kern_height/2;
  const int32 kern_half_width     = kern_width /2;
  const int32 kern_pixels         = kern_height * kern_width;
  const int32 min_num_good_pixels = kern_pixels/2;

  // Get X and Y derivatives of the input images
  ImageView<float> x_deriv = derivative_filter(left_image, 1, 0);
  ImageView<float> y_deriv = derivative_filter(left_image, 0, 1);
  ImageView<float> weight_template =
    vw::stereo::detail::compute_spatial_weight_image(kern_width, kern_height, two_sigma_sqr);

  // Workspace images are allocated up here out of the tight inner
  // loop.  We rasterize into these directly in the code below.
  ImageView<float> w(kern_width, kern_height);

  // Iterate over all of the pixels in the disparity map except for the outer edges.
  for ( int32 y = std::max(region_of_interest.min().y()-1,kern_half_height);
              y < std::min(left_image.rows()-kern_half_height,
                           region_of_interest.max().y()+1); ++y) {

    for (int32 x = std::max(region_of_interest.min().x()-1,kern_half_width);
               x < std::min(left_image.cols()-kern_half_width,
                            region_of_interest.max().x()+1); ++x) {

      BBox2i current_window(x-kern_half_width, y-kern_half_height,
                            kern_width, kern_height);

      // Skip over pixels for which we have no initial disparity estimate
      if ( !is_valid(disparity_map(x,y)) )
        continue;

      // Define and initialize the model params
      // Initialize our affine transform with the identity.  The
      // entries of d are laid out in row major order:
      //
      //   | d(0) d(1) d(2) |
      //   | d(3) d(4) d(5) |
      //   |  0    0    1   |
      //
      Vector6f d;
      d(0) = 1.0; d(1) = 0.0; d(2) = 0.0;
      d(3) = 0.0; d(4) = 1.0; d(5) = 0.0;
      float *dPtr = &(d[0]); // Raw data pointer access to avoid inlining failure

      // Compute the derivative image patches
      CropView<ImageView<float> > left_image_patch = crop(left_image, current_window);
      CropView<ImageView<float> > I_x = crop(x_deriv, current_window);
      CropView<ImageView<float> > I_y = crop(y_deriv, current_window);

      // Compute the base weight image
      int32 good_pixels = adjust_weight_image(w, crop(disparity_map, current_window), weight_template);

      // Skip over pixels for which there are very few good matches
      // in the neighborhood.
      if (good_pixels < min_num_good_pixels) {
        invalidate(disparity_map(x,y));
        continue;
      }

      //float curr_sum_I_e_val = 0.0;
      //float prev_sum_I_e_val = 0.0;

      // Iterate until a solution is found or the max number of
      // iterations is reached.
      for (unsigned iter = 0; iter < MAX_NUM_ITERATIONS; ++iter) {
        // First we check to see if our current subpixel translation
        // is less than one half of the window width.  If not, then
        // we are probably having trouble converging and we abort
        // this pixel!!
        if (norm_2( Vector2f(d[2],d[5]) ) > AFFINE_SUBPIXEL_MAX_TRANSLATION)
          break;

        float x_base = x + disparity_map(x,y)[0];
        float y_base = y + disparity_map(x,y)[1];

        //curr_sum_I_e_val = 0.0;
        
        Matrix6x6f rhs;
        Vector6f lhs;


        // Compute the outer range of xx and yy values and determine if everything will
        // fall within the bounds of right_interp_image
        // - If everything is safely in bounds, we can skip bounds checking in the main pixel loop below.
        float xVals[4], yVals[4];
        xVals[0] = dPtr[0]*(-kern_half_width) + dPtr[1]*  kern_half_height  + dPtr[2];
        yVals[0] = dPtr[3]*(-kern_half_width) + dPtr[4]*  kern_half_height  + dPtr[5];
        xVals[1] = dPtr[0]*(-kern_half_width) + dPtr[1]*(-kern_half_height) + dPtr[2];
        yVals[1] = dPtr[3]*(-kern_half_width) + dPtr[4]*(-kern_half_height) + dPtr[5];
        xVals[2] = dPtr[0]*  kern_half_width  + dPtr[1]*(-kern_half_height) + dPtr[2];
        yVals[2] = dPtr[3]*  kern_half_width  + dPtr[4]*(-kern_half_height) + dPtr[5];
        xVals[3] = dPtr[0]*  kern_half_width  + dPtr[1]*  kern_half_height  + dPtr[2];
        yVals[3] = dPtr[3]*  kern_half_width  + dPtr[4]*  kern_half_height  + dPtr[5];
        bool use_unsafe_interp = true;
        const float SAFETY_BUFFER = kern_half_width + kern_half_height; // Extra buffer to make sure we stay in bounds
        float minCol = -x_base + SAFETY_BUFFER;
        float maxCol = right_image.cols() - x_base - SAFETY_BUFFER;
        float minRow = -y_base + SAFETY_BUFFER;
        float maxRow = right_image.rows() - y_base - SAFETY_BUFFER;
        for (int i=0; i<4; ++i) { // If each corner of the kernel is in the safe bounds we should be safe
          if ( (xVals[i] >= maxCol) || (xVals[i] <= minCol) ||
               (yVals[i] >= maxRow) || (yVals[i] <= minRow)   ) {
            use_unsafe_interp = false;
            break;
          }
        }

        // Set up pixel accessors
        CropViewFAcc I_x_row = I_x.origin(), I_y_row = I_y.origin();
        CropViewTAcc left_image_patch_row = left_image_patch.origin();
        ImageViewFAcc w_row = w.origin();

        for (int32 jj = -kern_half_height; jj <= kern_half_height; ++jj) {
        
          ImageViewFAcc w_ptr = w_row;
          CropViewFAcc  I_x_ptr = I_x_row, I_y_ptr = I_y_row;
          CropViewTAcc  left_image_patch_ptr = left_image_patch_row;
          float xx_partial      = x_base + dPtr[1] * jj + dPtr[2]; // Compute outside inner loop for speed
          float yy_partial      = y_base + dPtr[4] * jj + dPtr[5];


          for (int32 ii = -kern_half_width; ii <= kern_half_width; ++ii) {
            // First we compute the pixel offset for the right image
            // and the error for the current pixel.
            float xx = dPtr[0] * ii + xx_partial;
            float yy = dPtr[3] * ii + yy_partial;

            // Avoid using the edge-extension view when possible.
            float interpreted_px;
            if (use_unsafe_interp)
              interpreted_px = right_interp_image_unsafe(xx,yy);
            else
              interpreted_px = right_interp_image(xx,yy);
            float I_e_val = interpreted_px - (*left_image_patch_ptr);
            
            
            // Apply the robust cost function.  We use a cauchy
            // function to gently remove outliers for small errors.
            //float thresh = 1e-4;

            // Cauchy seems to work well with thresh ~= 1e-4
            //float error_value   = fabsf(I_e_val);
            //float robust_weight = sqrtf(detail::cauchy_robust_coefficient(error_value,thresh))/error_value;

            //curr_sum_I_e_val += error_value; // Accumulated error value over pixels

            // Huber seems to work well with thresh >= 1e-5
            //        float robust_weight = sqrt(detail::huber_robust_coefficient(fabs(I_e_val),thresh))/fabs(I_e_val);

            // Disable robust cost function altogether
            float robust_weight = 1;
        
            // We combine the error value with the derivative and
            // add this to the update equation.
            float weight = robust_weight *(*w_ptr);
            
            //in_curr_sum_I_e_val += I_e_val;
           
            
            float I_x_val = weight  * (*I_x_ptr);
            float I_y_val = weight  * (*I_y_ptr);
            float I_x_sqr = I_x_val * (*I_x_ptr);
            float I_y_sqr = I_y_val * (*I_y_ptr);
            float I_x_I_y = I_x_val * (*I_y_ptr);

            // Left hand side
            float IxIe = I_x_val * I_e_val;
            float IyIe = I_y_val * I_e_val;
            lhs(0) -= ii * IxIe;
            lhs(1) -= jj * IxIe;
            lhs(2) -=      IxIe;
            lhs(3) -= ii * IyIe;
            lhs(4) -= jj * IyIe;
            lhs(5) -=      IyIe;

            float multipliers[3];
            multipliers[0] = ii*ii;
            multipliers[1] = ii*jj;
            multipliers[2] = jj*jj;

            float* rhsData = rhs.data(); // Access RHS by raw data pointer to avoid inlining failure
            
            // Right Hand Side UL
            rhsData[ 0] += multipliers[0] * I_x_sqr;
            rhsData[ 1] += multipliers[1] * I_x_sqr;
            rhsData[ 2] += ii    * I_x_sqr;
            rhsData[ 7] += multipliers[2] * I_x_sqr;
            rhsData[ 8] += jj    * I_x_sqr;
            rhsData[14] +=         I_x_sqr;

            // Right Hand Side UR
            rhsData[ 3] += multipliers[0] * I_x_I_y;
            rhsData[ 4] += multipliers[1] * I_x_I_y;
            rhsData[ 5] += ii    * I_x_I_y;
            rhsData[10] += multipliers[2] * I_x_I_y;
            rhsData[11] += jj    * I_x_I_y;
            rhsData[17] +=         I_x_I_y;

            // Right Hand Side LR
            rhsData[21] += multipliers[0] * I_y_sqr;
            rhsData[22] += multipliers[1] * I_y_sqr;
            rhsData[23] += ii    * I_y_sqr;
            rhsData[28] += multipliers[2] * I_y_sqr;
            rhsData[29] += jj    * I_y_sqr;
            rhsData[35] +=         I_y_sqr;

/*          For some reason these function calls are not being inlined!!!
            // Right Hand Side UL
            rhs(0,0) += multipliers[0] * I_x_sqr;
            rhs(0,1) += multipliers[1] * I_x_sqr;
            rhs(0,2) += ii    * I_x_sqr;
            rhs(1,1) += multipliers[2] * I_x_sqr;
            rhs(1,2) += jj    * I_x_sqr;
            rhs(2,2) +=         I_x_sqr;

            // Right Hand Side UR
            rhs(0,3) += multipliers[0] * I_x_I_y;
            rhs(0,4) += multipliers[1] * I_x_I_y;
            rhs(0,5) += ii    * I_x_I_y;
            rhs(1,4) += multipliers[2] * I_x_I_y;
            rhs(1,5) += jj    * I_x_I_y;
            rhs(2,5) +=         I_x_I_y;

            // Right Hand Side LR
            rhs(3,3) += multipliers[0] * I_y_sqr;
            rhs(3,4) += multipliers[1] * I_y_sqr;
            rhs(3,5) += ii    * I_y_sqr;
            rhs(4,4) += multipliers[2] * I_y_sqr;
            rhs(4,5) += jj    * I_y_sqr;
            rhs(5,5) +=         I_y_sqr;
*/

            I_x_ptr.next_col();
            I_y_ptr.next_col();
            left_image_patch_ptr.next_col();
          }
          I_x_row.next_row();
          I_y_row.next_row();
          left_image_patch_row.next_row();
        }

        // Fill in symmetric entries
        rhs(1,0) = rhs(0,1);
        rhs(2,0) = rhs(0,2);
        rhs(2,1) = rhs(1,2);
        rhs(3,0) = rhs(0,3);
        rhs(1,3) = rhs(3,1) = rhs(4,0) = rhs(0,4);
        rhs(2,3) = rhs(3,2) = rhs(5,0) = rhs(0,5);
        rhs(4,1) = rhs(1,4);
        rhs(2,4) = rhs(4,2) = rhs(5,1) = rhs(1,5);
        rhs(5,2) = rhs(2,5);
        rhs(4,3) = rhs(3,4);
        rhs(5,3) = rhs(3,5);
        rhs(5,4) = rhs(4,5);

        // Solves lhs = rhs * x, and stores the result in-place in lhs.
        const math::f77_int n = 6;
        const math::f77_int nrhs = 1;
        math::f77_int info;
        math::posv('L',n,nrhs,&(rhs(0,0)), n, &(lhs(0)), n, &info);

        //Termination
        //if (in_curr_sum_I_e_val < 0)
        //  in_curr_sum_I_e_val = - in_curr_sum_I_e_val;

        //curr_sum_I_e_val = in_curr_sum_I_e_val;

        d += lhs; // Update the affine transform
        //if (curr_sum_I_e_val < 0)
        //  curr_sum_I_e_val = - curr_sum_I_e_val;

        // Termination condition
        // - Quit if the change in the affine transform is tiny
        //if (norm_2(lhs) < 0.05) // If change in affine transform is small, quit the iteration loop
        //  break;                // - The value here strongly affects the results
        Vector6f weighted_lhs(lhs);
        const int32 kern_quarter_height = kern_half_height/2;
        const int32 kern_quarter_width  = kern_half_width /2;
        weighted_lhs[0] *= kern_quarter_width;
        weighted_lhs[1] *= kern_quarter_height;
        weighted_lhs[3] *= kern_quarter_width;
        weighted_lhs[4] *= kern_quarter_height;
        if (norm_2(weighted_lhs) < 0.05)
          break;

        /*
        // Termination condition
        // - Quit if the fit error starts increasing
        //if ((prev_sum_I_e_val < curr_sum_I_e_val) && (iter > 0))
        if ( (fabs(prev_sum_I_e_val - curr_sum_I_e_val) < 0.1)  && (iter > 0))
          break;
        else
          prev_sum_I_e_val = curr_sum_I_e_val;
        */
        
        
      } // End multiple iteration loop
      
      // If there is too much translation in our affine transform or we got NaNs, invalidate the pixel
      if ( norm_2( Vector2f(d[2],d[5]) ) >
           AFFINE_SUBPIXEL_MAX_TRANSLATION ||
           std::isnan(d[2]) || std::isnan(d[5]) )
        invalidate(disparity_map(x,y));
      else
        remove_mask(disparity_map(x,y)) += Vector2f(d[2],d[5]);
    } // X increment
  } // Y increment
}

/// Lucas-Kanade subpixel correlation function
/// - So far this one does not work that well.
void vw::stereo::subpixel_optimized_LK_2d(ImageView<PixelMask<Vector2f>> & disparity_map,
                                          ImageView<float> const& left_image,
                                          ImageView<float> const& right_image,
                                          int32  kern_width, int32 kern_height,
                                          BBox2i region_of_interest,
                                          bool   do_horizontal_subpixel,
                                          bool   do_vertical_subpixel) {

  typedef Vector<float, 2>    Vector2f;
  typedef Matrix<float, 2, 2> Matrix2x2f;
  typedef typename ImageView<float>::pixel_accessor ImageViewFAcc;
  typedef typename CropView<ImageView<float>>::pixel_accessor CropViewFAcc;
  typedef typename CropView<ImageView<float>>::pixel_accessor CropViewTAcc;

  // Bail out if no subpixel computation has been requested
  if (!do_horizontal_subpixel && !do_vertical_subpixel) return;

  // Fixed consts
  const unsigned MAX_NUM_ITERATIONS = 10;
  const float    two_sigma_sqr = 2.0*pow(float(kern_width)/5.0,2.0);

  VW_ASSERT( disparity_map.cols() == left_image.cols() &&
             disparity_map.rows() == left_image.rows(),
             ArgumentErr() << "subpixel_correlation: left image and "
             << "disparity map do not have the same dimensions.");

  // Interpolated Input Images
  InterpolationView<EdgeExtensionView<ImageView<float>, ZeroEdgeExtension>, BilinearInterpolation> right_interp_image =
         interpolate(right_image, BilinearInterpolation(), ZeroEdgeExtension());

  // This is the maximum number of pixels that the solution can be
  // adjusted by subpixel refinement.
  float SUBPIXEL_MAX_TRANSLATION = kern_width/2;

  const int32 kern_half_height    = kern_height/2;
  const int32 kern_half_width     = kern_width /2;
  const int32 kern_pixels         = kern_height * kern_width;
  const int32 min_num_good_pixels = kern_pixels/2;

  // Get X and Y derivatives of the input images
  ImageView<float> x_deriv = derivative_filter(left_image, 1, 0);
  ImageView<float> y_deriv = derivative_filter(left_image, 0, 1);
  ImageView<float> weight_template = vw::stereo::detail::compute_spatial_weight_image(kern_width, kern_height, two_sigma_sqr);

  // Workspace images are allocated up here out of the tight inner
  // loop.  We rasterize into these directly in the code below.
  ImageView<float> w(kern_width, kern_height);

  // Iterate over all of the pixels in the disparity map except for the outer edges.
  for ( int32 y = std::max(region_of_interest.min().y()-1,kern_half_height);
              y < std::min(left_image.rows()-kern_half_height,
                           region_of_interest.max().y()+1); ++y) {

    for (int32 x = std::max(region_of_interest.min().x()-1,kern_half_width);
               x < std::min(left_image.cols()-kern_half_width,
                            region_of_interest.max().x()+1); ++x) {

      BBox2i current_window(x-kern_half_width, y-kern_half_height,
                            kern_width, kern_height);

      // Skip over pixels for which we have no initial disparity estimate
      if ( !is_valid(disparity_map(x,y)) )
        continue;

      // We are just solving for a simple translation vector
      Vector2f d;
      d(0) = 0.0; d(1) = 0.0;

      // Compute the derivative image patches
      CropView<ImageView<float> > left_image_patch = crop(left_image, current_window);
      CropView<ImageView<float> > I_x = crop(x_deriv, current_window);
      CropView<ImageView<float> > I_y = crop(y_deriv, current_window);

      // Compute the base weight image
      int32 good_pixels = adjust_weight_image(w, crop(disparity_map, current_window), weight_template);

      // Skip over pixels for which there are very few good matches
      // in the neighborhood.
      if (good_pixels < min_num_good_pixels) {
        invalidate(disparity_map(x,y));
        continue;
      }

      // Iterate until a solution is found or the max number of
      // iterations is reached.
      for (unsigned iter = 0; iter < MAX_NUM_ITERATIONS; ++iter) {
        // First we check to see if our current subpixel translation
        // is less than one half of the window width.  If not, then
        // we are probably having trouble converging and we abort
        // this pixel!!
        if (norm_2(d) > SUBPIXEL_MAX_TRANSLATION)
          break;

        float x_base = x + disparity_map(x,y)[0];
        float y_base = y + disparity_map(x,y)[1];

        Matrix2x2f rhs;
        Vector2f lhs;
   
        // Set up pixel accessors
        CropViewFAcc I_x_row = I_x.origin(), I_y_row = I_y.origin();
        CropViewTAcc left_image_patch_row = left_image_patch.origin();
        ImageViewFAcc w_row = w.origin();

        // Perform loop that does Expectation and Maximization in one go
        for (int32 jj = -kern_half_height; jj <= kern_half_height; ++jj) {
        
          ImageViewFAcc w_ptr = w_row;
          CropViewFAcc  I_x_ptr = I_x_row, I_y_ptr = I_y_row;
          CropViewTAcc  left_image_patch_ptr = left_image_patch_row;
          
          float xx_partial = x_base + d[0]; // Compute outside inner loop for speed
          
          float yy = y_base + jj + d[1];

          for (int32 ii = -kern_half_width; ii <= kern_half_width; ++ii) {
            // First we compute the pixel offset for the right image
            // and the error for the current pixel.
            float xx = ii + xx_partial;

            /// Expectation
            float interpreted_px = right_interp_image(xx,yy);
            float I_e_val = interpreted_px - (*left_image_patch_ptr);
            
            // Disable robust cost function altogether
            float robust_weight = 1;
        
            // We combine the error value with the derivative and
            // add this to the update equation.
            float weight = robust_weight *(*w_ptr);
        
            
            float I_x_val = weight  * (*I_x_ptr);
            float I_y_val = weight  * (*I_y_ptr);
            float I_x_sqr = I_x_val * (*I_x_ptr);
            float I_y_sqr = I_y_val * (*I_y_ptr);
            float I_x_I_y = I_x_val * (*I_y_ptr);

            // Left hand side
            lhs(0) -= I_x_val * I_e_val;
            lhs(1) -= I_y_val * I_e_val;

            float* rhsData = rhs.data(); // Access RHS by raw data pointer to avoid inlining failure

            // Right Hand Side UL
            rhsData[0] += I_x_sqr;

            // Right Hand Side UR
            rhsData[1] += I_x_I_y;

            // Right Hand Side LR
            rhsData[3] += I_y_sqr;

            I_x_ptr.next_col();
            I_y_ptr.next_col();
            left_image_patch_ptr.next_col();
          }
          I_x_row.next_row();
          I_y_row.next_row();
          left_image_patch_row.next_row();
        }

        // Fill in symmetric entries
        rhs(1,0) = rhs(0,1);

        // Solves lhs = rhs * x, and stores the result in-place in lhs.
        const math::f77_int n    = 2;
        const math::f77_int nrhs = 1;
        math::f77_int info;
        math::posv('L',n,nrhs,&(rhs(0,0)), n, &(lhs(0)), n, &info);

        d += lhs; // Update the affine transform

        // Termination condition
        // - Quit if the change in the affine transform is tiny
        if (norm_2(lhs) < 0.05) // If change in affine transform is small, quit the iteration loop
          break;
        
      } // End multiple iteration loop
      
      // If there is too much translation in our affine transform or we got NaNs, invalidate the pixel
      if ( norm_2(d) > SUBPIXEL_MAX_TRANSLATION ||
           std::isnan(d[0]) || std::isnan(d[1]) )
        invalidate(disparity_map(x,y));
      else
        remove_mask(disparity_map(x,y)) += d;
    } // X increment
  } // Y increment
}

int vw::stereo::adjust_weight_image(ImageView<float> &weight,
                                    ImageView<PixelMask<Vector2f>> const& disparity_map_patch,
                                    ImageView<float> const& weight_template) {
  float sum = 0;
  int32 num_good_pix = 0;
  typedef ImageView<float>::pixel_accessor IViewFAcc;
  typedef ImageView<PixelMask<Vector2f> >::pixel_accessor IViewDAcc;
  IViewFAcc weight_row_acc = weight.origin();
  IViewFAcc template_row_acc = weight_template.origin();
  IViewDAcc disp_row_acc = disparity_map_patch.origin();
  for (int32 j = 0; j < weight_template.rows(); ++j) {
    IViewFAcc weight_col_acc = weight_row_acc;
    IViewFAcc template_col_acc = template_row_acc;
    IViewDAcc disp_col_acc = disp_row_acc;
    for (int32 i = 0; i < weight_template.cols(); ++i ) {

      // Mask is zero if the disparity map's pixel is missing...
      if ( !is_valid(*disp_col_acc) )
        *weight_col_acc = 0;

      // ... otherwise we use the weight from the weight template
      else {
        *weight_col_acc = *template_col_acc;
        sum += *weight_col_acc;
        ++num_good_pix;
      }

      disp_col_acc.next_col();
      weight_col_acc.next_col();
      template_col_acc.next_col();
    }
    disp_row_acc.next_row();
    weight_row_acc.next_row();
    template_row_acc.next_row();
  }

  // Normalize the weight image
  if (sum == 0)
    vw_throw(LogicErr() << "subpixel_weight: Sum of weight image was zero. This isn't supposed to happen!");
  else
    weight /= sum;
  return num_good_pix;
}

template <class ImageT1, class ImageT2>
void vw::stereo::cross_corr_consistency_check(
    ImageViewBase<ImageT1> const& l2r,
    ImageViewBase<ImageT2> const& r2l,
    float cross_corr_threshold,
    ImageView<PixelMask<float>> * lr_disp_diff,
    Vector2i ul_corner_offset,
    bool verbose) {
  int32 l2r_rows = l2r.impl().rows(), l2r_cols = l2r.impl().cols(),
    r2l_rows = r2l.impl().rows(), r2l_cols = r2l.impl().cols();
  size_t count = 0, match_count = 0;

  if (verbose)
    vw_out(VerboseDebugMessage, "stereo") << "\tCrosscorr threshold: "
                                          << cross_corr_threshold << "\n";
  VW_DEBUG_ASSERT(cross_corr_threshold >= 0.0,
                  ArgumentErr()
                    << "cross_corr_consistency_check: "
                    << "the threshold is less than 0.");

  typename ImageT1::pixel_accessor l2r_row = l2r.impl().origin();
  for (int32 r = 0; r < l2r_rows; r++) {
    typename ImageT1::pixel_accessor l2r_col = l2r_row;
    for (int32 c = 0; c < l2r_cols; c++) {

      int32 r2l_x = c + (*l2r_col)[0];
      int32 r2l_y = r + (*l2r_col)[1];

      if (r2l_x < 0 || r2l_x >= r2l_cols ||
          r2l_y < 0 || r2l_y >= r2l_rows) {
        invalidate(*l2r_col);
      } else if (!is_valid(*l2r_col) ||
                 !is_valid(r2l.impl()(r2l_x, r2l_y))) {
        invalidate(*l2r_col);
      } else {
        float disp_diff =
          std::max(fabs((*l2r_col)[0] + r2l.impl()(r2l_x, r2l_y)[0]),
                   fabs((*l2r_col)[1] + r2l.impl()(r2l_x, r2l_y)[1]));
        if (cross_corr_threshold >= disp_diff) {
          count++;
          match_count++;
          if (lr_disp_diff != NULL)
            (*lr_disp_diff)(c + ul_corner_offset[0],
                            r + ul_corner_offset[1])
              = PixelMask<float>(disp_diff);
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
    vw_out(VerboseDebugMessage, "stereo")
      << "\tCross-correlation retained " << count
      << " / " << match_count << " matches ("
      << ((float)count / match_count * 100) << " percent).\n";
}

// Explicit instantiations for the types used in production and tests
typedef PixelMask<Vector2i> PMV2i;
typedef PixelMask<Vector2f> PMV2f;

template void vw::stereo::cross_corr_consistency_check(
    ImageViewBase<CropView<ImageView<PMV2i>>> const&,
    ImageViewBase<ImageView<PMV2i>> const&,
    float, ImageView<PixelMask<float>>*, Vector2i, bool);

template void vw::stereo::cross_corr_consistency_check(
    ImageViewBase<ImageView<PMV2f>> const&,
    ImageViewBase<ImageView<PMV2f>> const&,
    float, ImageView<PixelMask<float>>*, Vector2i, bool);
