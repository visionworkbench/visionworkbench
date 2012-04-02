// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


namespace detail {
  // is_equal_elements, tells if all the elements are equal to each other
  template <class T>
  struct IsEqualElements {
    bool result;

    IsEqualElements() : result(true), m_comparing(false) {}

    void operator()( T const& i ) {
      if ( m_comparing ) {
        result = result && (*m_first_element == i);
      } else {
        m_comparing = true;
        m_first_element = &i;
      }
    }
  protected:
    bool m_comparing;
    T const* m_first_element;
  };

  ///-------------------------------------------------------------------------

  inline ImageView<float>
  compute_spatial_weight_image(int32 kern_width, int32 kern_height,
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

}

template<class ChannelT> void
subpixel_correlation_affine_2d_EM(ImageView<PixelMask<Vector2f> > &disparity_map,
                                  ImageView<ChannelT> const& left_image,
                                  ImageView<ChannelT> const& right_image,
                                  int32 kern_width, int32 kern_height,
                                  BBox2i region_of_interest,
                                  bool do_horizontal_subpixel,
                                  bool do_vertical_subpixel,
                                  bool /*verbose*/ ) {
  typedef Vector<float,6> Vector6f;
  typedef Matrix<float,6,6> Matrix6x6f;
  typedef typename ImageView<float>::pixel_accessor ImageViewFAcc;
  typedef typename CropView<ImageView<float> >::pixel_accessor CropViewFAcc;
  typedef typename CropView<ImageView<ChannelT> >::pixel_accessor CropViewTAcc;

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
  InterpolationView<EdgeExtensionView<ImageView<ChannelT>, ZeroEdgeExtension>, BilinearInterpolation> right_interp_image =
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
    detail::compute_spatial_weight_image(kern_width, kern_height, two_sigma_sqr);

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
      CropView<ImageView<ChannelT> > left_image_patch = crop(left_image, current_window);
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
        float w_plane = 0.8;
        float w_noise = 0.2;
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
              ChannelT interpreted_px = right_interp_image(xx,yy);
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
            solve_symmetric_nocopy(rhs,lhs);
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

          ChannelT interpreted_px = right_interp_image(xx,yy);
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
template<class ChannelT> void
subpixel_optimized_affine_2d_EM(ImageView<PixelMask<Vector2f> > &disparity_map,
                                ImageView<ChannelT> const& left_image,
                                ImageView<ChannelT> const& right_image,
                                int32 kern_width, int32 kern_height,
                                BBox2i region_of_interest,
                                bool do_horizontal_subpixel,
                                bool do_vertical_subpixel,
                                bool /*verbose*/ ) {
  typedef Vector<float,6> Vector6f;
  typedef Matrix<float,6,6> Matrix6x6f;
  typedef typename ImageView<float>::pixel_accessor ImageViewFAcc;
  typedef typename CropView<ImageView<float> >::pixel_accessor CropViewFAcc;
  typedef typename CropView<ImageView<ChannelT> >::pixel_accessor CropViewTAcc;

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
  InterpolationView<EdgeExtensionView<ImageView<ChannelT>, ZeroEdgeExtension>, BilinearInterpolation> right_interp_image =
    interpolate(right_image, BilinearInterpolation(), ZeroEdgeExtension());

  // This is the maximum number of pixels that the solution can be
  // adjusted by affine subpixel refinement.
  float AFFINE_SUBPIXEL_MAX_TRANSLATION = kern_width/2;

  const int32 kern_half_height = kern_height/2;
  const int32 kern_half_width = kern_width/2;
  const int32 kern_pixels = kern_height * kern_width;
  const int32 weight_threshold = kern_pixels/2;

  ImageView<float> x_deriv = derivative_filter(left_image, 1, 0);
  ImageView<float> y_deriv = derivative_filter(left_image, 0, 1);
  ImageView<float> weight_template =
    detail::compute_spatial_weight_image(kern_width, kern_height, two_sigma_sqr);

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
      CropView<ImageView<ChannelT> > left_image_patch = crop(left_image, current_window);
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

        //set init params - START
        float var2_plane = 1e-3;
        float mean_noise = 0.0;
        float var2_noise = 1e-2;
        float w_plane = 0.8;
        float w_noise = 0.2;
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

          // Counts pixels skipped on evaluation.
          int32 skip = 0;

          // Perform loop that does Expectation and Maximization in one go
          for (int32 jj = -kern_half_height; jj <= kern_half_height; ++jj) {
            ImageViewFAcc w_ptr = w_row;
            CropViewFAcc I_x_ptr = I_x_row, I_y_ptr = I_y_row;
            CropViewTAcc left_image_patch_ptr = left_image_patch_row;
            float xx_partial = x_base + d[1] * jj + d[2];
            float yy_partial = y_base + d[4] * jj + d[5];
            float delta_x_partial = d_em[1] * jj + d_em[2];
            float delta_y_partial = d_em[4] * jj + d_em[5];

            for (int32 ii = -kern_half_width; ii <= kern_half_width; ++ii) {
              // First we compute the pixel offset for the right image
              // and the error for the current pixel.
              float xx = d[0] * ii + xx_partial;
              float yy = d[3] * ii + yy_partial;
              float delta_x = d_em[0] * ii + delta_x_partial;
              float delta_y = d_em[3] * ii + delta_y_partial;

              /// Expectation
              ChannelT interpreted_px = right_interp_image(xx,yy);
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
              float gamma_plane = plane_prob*w_plane/sum;
              float gamma_noise = noise_prob*w_noise/sum;
              // End Expectation

              // Maximization (computing the d_em vector)
              mean_noise_tmp +=
                interpreted_px * gamma_noise;
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
              float I_x_val = weight * (*I_x_ptr);
              float I_y_val = weight * (*I_y_ptr);
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
