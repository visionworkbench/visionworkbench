// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Stereo/CorrelateResearch.h>
#include <vw/Math/Matrix.h>

#include <vw/Core/Stopwatch.h>
#include <vw/Math/LinearAlgebra.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/Interpolation.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/Filter.h>
#include <cmath>
#include <ctime>
#include <vw/FileIO.h>

namespace vw {
namespace stereo {

  namespace detail {

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
  }

  template<class ChannelT>
  void subpixel_correlation_affine_2d(ImageView<PixelMask<Vector2f> > &disparity_map,
                                      ImageView<ChannelT> const& left_input_image,
                                      ImageView<ChannelT> const& right_input_image,
                                      int kern_width, int kern_height,
                                      bool do_horizontal_subpixel,
                                      bool do_vertical_subpixel,
                                      bool verbose = false) {

    VW_ASSERT( disparity_map.cols() == left_input_image.cols() &&
               disparity_map.rows() == left_input_image.rows(),
               ArgumentErr() << "subpixel_correlation: left image and disparity map do not have the same dimensions.");


    for (float blur_sigma = 3; blur_sigma >= 1.0; blur_sigma /= 2.0) {
      ImageView<ChannelT> left_image = LogStereoPreprocessingFilter(blur_sigma)(left_input_image);
      ImageView<ChannelT> right_image = LogStereoPreprocessingFilter(blur_sigma)(right_input_image);

      // This is the maximum number of pixels that the solution can be
      // adjusted by affine subpixel refinement.
      float AFFINE_SUBPIXEL_MAX_TRANSLATION = kern_width/2;
      int kern_half_height = kern_height/2;
      int kern_half_width = kern_width/2;

      // Robust cost function settings
      const float thresh = 0.01;
      detail::HuberError robust_cost_fn(thresh);

      int kern_pixels = kern_height * kern_width;
      int weight_threshold = kern_pixels / 2;

      // Bail out if no subpixel computation has been requested
      if (!do_horizontal_subpixel && !do_vertical_subpixel) return;

      ImageView<float> x_deriv = derivative_filter(left_image, 1, 0);
      ImageView<float> y_deriv = derivative_filter(left_image, 0, 1);
      ImageView<float> weight_template = detail::compute_gaussian_weight_image(kern_width, kern_height);

      // Workspace images are allocated up here out of the tight inner
      // loop.  We rasterize into these directly in the code below.
      ImageView<float> w(kern_width, kern_height);

      // Iterate over all of the pixels in the disparity map except for
      // the outer edges.
      Stopwatch sw;
      sw.start();
      double last_time = 0;

      for ( int y = kern_half_height;
            y < left_image.rows()-kern_half_height; ++y) {
        if (verbose && y % 10 == 0) {
          sw.stop();
          vw_out(InfoMessage, "stereo") << "\tProcessing subpixel line: " << y << " / " << left_image.rows()
                                        << "    (" << (10 * left_image.cols() / (sw.elapsed_seconds() - last_time))
                                        << " pixels/s, "<< sw.elapsed_seconds() << " s total )      \r" << std::flush;
          last_time = sw.elapsed_seconds();
          sw.start();
        }

        for ( int x = kern_half_width;
              x < left_image.cols()-kern_half_width; ++x) {

          BBox2i current_window(x-kern_half_width, y-kern_half_height, kern_width, kern_height);
          Vector2f base_offset(-disparity_map(x,y).child());

          // Skip over pixels for which we have no initial disparity estimate
          if ( !is_valid(disparity_map(x,y)) )
            continue;

          // Initialize our affine transform with the identity.  The
          // entries of d are laid out in row major order:
          //
          //   | d(0) d(1) d(2) |
          //   | d(3) d(4) d(5) |
          //   |  0    0    1   |
          //
          Vector<float,6> d;
          d(0) = 1.0;
          d(4) = 1.0;

          // Compute the derivative image patches
          CropView<ImageView<ChannelT> > left_image_patch = crop(left_image, current_window);
          CropView<ImageView<float> > I_x = crop(x_deriv, current_window);
          CropView<ImageView<float> > I_y = crop(y_deriv, current_window);

          // Compute the base weight image
          int good_pixels = adjust_weight_image(w, crop(disparity_map, current_window), weight_template);

          // Skip over pixels for which there are very few good matches
          // in the neighborhood.
          if (good_pixels < weight_threshold) {
            invalidate( disparity_map(x,y) );
            continue;
          }

          // Iterate until a solution is found or the max number of
          // iterations is reached.
          for (unsigned iter = 0; iter < 10; ++iter) {
            // First we check to see if our current subpixel translation
            // is less than one half of the window width.  If not, then
            // we are probably having trouble converging and we abort
            // this pixel!!
            if (norm_2( Vector<float,2>(d[2],d[5]) ) > AFFINE_SUBPIXEL_MAX_TRANSLATION)
              break;

            InterpolationView<EdgeExtensionView<ImageView<ChannelT>, ZeroEdgeExtension>, BilinearInterpolation> right_interp_image =
              interpolate(right_image, BilinearInterpolation(), ZeroEdgeExtension());

            float x_base = x + disparity_map(x,y)[0];
            float y_base = y + disparity_map(x,y)[1];

            Matrix<float,6,6> rhs;
            Vector<float,6> lhs;

            // Set up pixel accessors
            typename ImageView<float>::pixel_accessor w_row = w.origin();
            typename CropView<ImageView<float> >::pixel_accessor I_x_row = I_x.origin();
            typename CropView<ImageView<float> >::pixel_accessor I_y_row = I_y.origin();
            typename CropView<ImageView<ChannelT> >::pixel_accessor left_image_patch_row = left_image_patch.origin();

            for (int jj = -kern_half_height; jj <= kern_half_height; ++jj) {
              typename ImageView<float>::pixel_accessor w_ptr = w_row;
              typename CropView<ImageView<float> >::pixel_accessor I_x_ptr = I_x_row;
              typename CropView<ImageView<float> >::pixel_accessor I_y_ptr = I_y_row;
              typename CropView<ImageView<ChannelT> >::pixel_accessor left_image_patch_ptr = left_image_patch_row;

              for (int ii = -kern_half_width; ii <= kern_half_width; ++ii) {

                // First we compute the pixel offset for the right image
                // and the error for the current pixel.
                float xx = x_base + d[0] * ii + d[1] * jj + d[2];
                float yy = y_base + d[3] * ii + d[4] * jj + d[5];
                float I_e_val = right_interp_image(xx,yy) - (*left_image_patch_ptr) + 1e-16;
                //              error_total += pow(I_e_val,2);

                // Apply the robust cost function.  We use a cauchy
                // function to gently remove outliers for small errors.
                float thresh = 1e-3;

                // Cauchy seems to work well with thresh ~= 1e-4
                float error_value = fabsf(I_e_val);
                float robust_weight = sqrtf(detail::cauchy_robust_coefficient(error_value,thresh))/error_value;

                // Huber seems to work well with thresh >= 1e-5
                //        float robust_weight = sqrt(detail::huber_robust_coefficient(fabs(I_e_val),thresh))/fabs(I_e_val);

                // Disable robust cost function altogether
                //        float robust_weight = 1;

                // We combine the error value with the derivative and
                // add this to the update equation.
                float weight = robust_weight *(*w_ptr);
                //float weight = robust_weight;// *(*w_ptr);
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

                w_ptr.next_col();
                I_x_ptr.next_col();
                I_y_ptr.next_col();
                left_image_patch_ptr.next_col();
              }
              w_row.next_row();
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
            //           Matrix<double,6,6> pre_rhs = rhs;
            //           Vector<double,6> pre_lhs = lhs;
            try {
              solve_symmetric_nocopy(rhs,lhs);
            } catch (ArgumentErr &/*e*/) {} // Do nothing
            d += lhs;

            // Termination condition
            if (norm_2(lhs) < 0.01)
              break;
          }
          if ( norm_2( Vector2f(d[2],d[5]) ) >
               AFFINE_SUBPIXEL_MAX_TRANSLATION ||
               std::isnan(d[2]) || std::isnan(d[5]) )
              invalidate(disparity_map(x,y));
          else
            remove_mask(disparity_map(x,y)) += Vector2f(d[2],d[5]);
        }
      }

    }

    if (verbose)
      vw_out(InfoMessage, "stereo") << "\tProcessing subpixel line: done.                                         \n";
  }

  template<class ChannelT> void
  subpixel_correlation_affine_2d_bayesian(ImageView<PixelMask<Vector2f> > &disparity_map,
                                          ImageView<ChannelT> const& left_input_image,
                                          ImageView<ChannelT> const& right_input_image,
                                          int kern_width, int kern_height,
                                          bool do_horizontal_subpixel,
                                          bool do_vertical_subpixel,
                                          bool verbose) {

    VW_ASSERT( disparity_map.cols() == left_input_image.cols() &&
               disparity_map.rows() == left_input_image.rows(),
               ArgumentErr() << "subpixel_correlation: left image and disparity map do not have the same dimensions.");

    for (float blur_sigma = 3; blur_sigma >= 1.0; blur_sigma /= 2.0) {
      ImageView<ChannelT> left_image = LogStereoPreprocessingFilter(blur_sigma)(left_input_image);
      ImageView<ChannelT> right_image = LogStereoPreprocessingFilter(blur_sigma)(right_input_image);

      // This is the maximum number of pixels that the solution can be
      // adjusted by affine subpixel refinement.
      float AFFINE_SUBPIXEL_MAX_TRANSLATION = kern_width/2;
      //float AFFINE_SUBPIXEL_MAX_TRANSLATION = kern_width;
      int kern_half_height = kern_height/2;
      int kern_half_width = kern_width/2;

      // Robust cost function settings
      //const float thresh = 0.01;
      //detail::HuberError robust_cost_fn(thresh);

      int kern_pixels = kern_height * kern_width;
      int weight_threshold = kern_pixels / 2;

      // Bail out if no subpixel computation has been requested
      if (!do_horizontal_subpixel && !do_vertical_subpixel) return;

      ImageView<float> x_deriv = derivative_filter(left_image, 1, 0);
      ImageView<float> y_deriv = derivative_filter(left_image, 0, 1);
      ImageView<float> weight_template = detail::compute_gaussian_weight_image(kern_width, kern_height);

      // Workspace images are allocated up here out of the tight inner
      // loop.  We rasterize into these directly in the code below.
      ImageView<float> w(kern_width, kern_height);

      // Iterate over all of the pixels in the disparity map except for
      // the outer edges.
      Stopwatch sw;
      sw.start();
      double last_time = 0;

      for (int y=kern_half_height; y<left_image.rows()-kern_half_height; ++y) {
        if (verbose && y % 10 == 0) {
          sw.stop();
          vw_out(InfoMessage, "stereo") << "\tProcessing subpixel line: " << y << " / " << left_image.rows() << "    (" << (10 * left_image.cols() / (sw.elapsed_seconds() - last_time)) << " pixels/s, "<< sw.elapsed_seconds() << " s total )      \r" << std::flush;
          last_time = sw.elapsed_seconds();
          sw.start();
        }

        for (int x=kern_half_width; x<left_image.cols()-kern_half_width; ++x) {
          BBox2i current_window(x-kern_half_width, y-kern_half_height, kern_width, kern_height);
          Vector2 base_offset(-disparity_map(x,y).child());

          // Skip over pixels for which we have no initial disparity estimate
          if ( !is_valid(disparity_map(x,y)) )
            continue;

          // Initialize our affine transform with the identity.  The
          // entries of d are laid out in row major order:
          //
          //   | d(0) d(1) d(2) |
          //   | d(3) d(4) d(5) |
          //   |  0    0    1   |
          //
          Vector<float,6> d;
          d(0) = 1.0;
          d(4) = 1.0;

          // Compute the derivative image patches
          CropView<ImageView<ChannelT> > left_image_patch = crop(left_image, current_window);
          CropView<ImageView<float> > I_x = crop(x_deriv, current_window);
          CropView<ImageView<float> > I_y = crop(y_deriv, current_window);

          // Compute the base weight image
          int good_pixels = adjust_weight_image(w, crop(disparity_map, current_window), weight_template);

          // Skip over pixels for which there are very few good matches
          // in the neighborhood.
          if (good_pixels < weight_threshold) {
            invalidate( disparity_map(x,y) );
            continue;
          }

          float curr_sum_I_e_val = 0.0;
          float prev_sum_I_e_val = 0.0;
          unsigned iter;
          // Iterate until a solution is found or the max number of
          // iterations is reached.
          for (iter = 0; iter < 10; ++iter) {
            // First we check to see if our current subpixel translation
            // is less than one half of the window width.  If not, then
            // we are probably having trouble converging and we abort
            // this pixel!!
            if (norm_2( Vector<float,2>(d[2],d[5]) ) > AFFINE_SUBPIXEL_MAX_TRANSLATION)
              break;

            InterpolationView<EdgeExtensionView<ImageView<ChannelT>, ZeroEdgeExtension>, BilinearInterpolation> right_interp_image =
              interpolate(right_image, BilinearInterpolation(), ZeroEdgeExtension());

            float x_base = x + disparity_map(x,y)[0];
            float y_base = y + disparity_map(x,y)[1];
            //          float error_total = 0;

            Matrix<float,6,6> rhs;
            Vector<float,6> lhs;

            ImageView<float> ll_value(kern_width, kern_height);
            float sum_error_value = 0;
            float mean_l = 0.0;
            float mean_r = 0.0;

            // Set up pixel accessors
            typename ImageView<float>::pixel_accessor w_row = w.origin();
            typename CropView<ImageView<float> >::pixel_accessor I_x_row = I_x.origin();
            typename CropView<ImageView<float> >::pixel_accessor I_y_row = I_y.origin();
            typename CropView<ImageView<ChannelT> >::pixel_accessor left_image_patch_row = left_image_patch.origin();


            for (int jj = -kern_half_height; jj <= kern_half_height; ++jj) {
              typename ImageView<float>::pixel_accessor w_ptr = w_row;
              typename CropView<ImageView<float> >::pixel_accessor I_x_ptr = I_x_row;
              typename CropView<ImageView<float> >::pixel_accessor I_y_ptr = I_y_row;
              typename CropView<ImageView<ChannelT> >::pixel_accessor left_image_patch_ptr = left_image_patch_row;

              for (int ii = -kern_half_width; ii <= kern_half_width; ++ii) {

                // First we compute the pixel offset for the right image
                // and the error for the current pixel.
                float xx = x_base + d[0] * ii + d[1] * jj + d[2];
                float yy = y_base + d[3] * ii + d[4] * jj + d[5];
                float I_e_val = right_interp_image(xx,yy) - (*left_image_patch_ptr);// + 1e-16;
                //              error_total += pow(I_e_val,2);

                // Apply the robust cost function.  We use a huber
                // function to gently remove outliers for small errors,
                // but we set a hard limit a 5 times the cost threshold
                // to remove major (salt&pepper) noise.
                //float thresh = 1e-3;

                float two_sigma_2 = 1e-4;//1e-3;//1e-4;//1e-5;//1e-6;

                ll_value(jj+kern_half_height, ii+kern_half_width) =  exp(-1*(I_e_val*I_e_val)/two_sigma_2);
                sum_error_value = sum_error_value + ll_value(jj+kern_half_height, ii+kern_half_width);

                mean_l = mean_l + (*left_image_patch_ptr);
                mean_r = mean_r + right_interp_image(xx,yy);

                w_ptr.next_col();
                I_x_ptr.next_col();
                I_y_ptr.next_col();
                left_image_patch_ptr.next_col();
              }
              w_row.next_row();
              I_x_row.next_row();
              I_y_row.next_row();
              left_image_patch_row.next_row();
            }

            mean_r = mean_r/(kern_width*kern_height);
            mean_l = mean_l/(kern_width*kern_height);

            // Set up pixel accessors
            w_row = w.origin();
            I_x_row = I_x.origin();
            I_y_row = I_y.origin();
            left_image_patch_row = left_image_patch.origin();

            curr_sum_I_e_val = 0.0;

            for (int jj = -kern_half_height; jj <= kern_half_height; ++jj) {
              typename ImageView<float>::pixel_accessor w_ptr = w_row;
              typename CropView<ImageView<float> >::pixel_accessor I_x_ptr = I_x_row;
              typename CropView<ImageView<float> >::pixel_accessor I_y_ptr = I_y_row;
              typename CropView<ImageView<ChannelT> >::pixel_accessor left_image_patch_ptr = left_image_patch_row;

              for (int ii = -kern_half_width; ii <= kern_half_width; ++ii) {

                // First we compute the pixel offset for the right image
                // and the error for the current pixel.
                float xx = x_base + d[0] * ii + d[1] * jj + d[2];
                float yy = y_base + d[3] * ii + d[4] * jj + d[5];
                float I_e_val = right_interp_image(xx,yy) - (*left_image_patch_ptr);
                curr_sum_I_e_val = curr_sum_I_e_val + I_e_val;

                float robust_weight = ll_value(jj+kern_half_height, ii+kern_half_width)/sum_error_value;

                // We combine the error value with the derivative and
                // add this to the update equation.
                float weight = robust_weight*(*w_ptr);//*ll;
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

                w_ptr.next_col();
                I_x_ptr.next_col();
                I_y_ptr.next_col();
                left_image_patch_ptr.next_col();
              }
              w_row.next_row();
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
              solve_symmetric_nocopy(rhs,lhs);
            } catch (ArgumentErr &/*e*/) {} // Do nothing
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
            invalidate( disparity_map(x,y) );
          else
            remove_mask(disparity_map(x,y)) += Vector2f(d[2],d[5]);
        }
      }
    }
  }

}}
