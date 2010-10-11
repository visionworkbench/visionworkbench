// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_STEREO_AFFINE_MIXTURE_COMPONENT__
#define __VW_STEREO_AFFINE_MIXTURE_COMPONENT__

#include "MixtureComponent.h"

#include <vw/Image.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Math.h>
#include <iostream>
#include <cmath>

namespace vw {
  namespace stereo {
    class AffineTransformOrigin : public TransformHelper<AffineTransformOrigin,ConvexFunction,ConvexFunction> {
    protected:
      double a, b, c, d;
      double ai, bi, ci, di;
      double x, y;
      double x0, y0;

    public:

    AffineTransformOrigin() : a(1), b(0), c(0), d(1), x(0), y(0), x0(0), y0(0) {  // defualt is the identity transform
        ai = 1;
        bi = 0;
        ci = 0;
        di = 1;
      }

    AffineTransformOrigin( Matrix2x2 const& matrix, Vector2 const& offset, Vector2 const& origin) : a( matrix(0,0) ), b( matrix(0,1) ), c( matrix(1,0) ), d( matrix(1,1))
        {
          Matrix2x2 inv = inverse( matrix );
          ai = inv(0,0);
          bi = inv(0,1);
          ci = inv(1,0);
          di = inv(1,1);

          x = offset(0);
          y = offset(1);
          x0 = origin(0);
          y0 = origin(1);
        }
      // forward
      inline Vector2 reverse( Vector2 const& p ) const {
        return Vector2(a*(p[0]-x0)+b*(p[1]-y0) + x + x0,
                       c*(p[0]-x0)+d*(p[1]-y0) + y + y0);
      }
      // reverse
      inline Vector2 forward( Vector2 const& p ) const {
        double px = p[0]-x -x0;
        double py = p[1]-y -y0;
        return Vector2(ai*px+bi*py+x0,
                       ci*px+di*py+y0);
      }
    };



    template<class PixelT, class PrecisionT>
    class AffineMixtureComponent : public MixtureComponentBase<AffineMixtureComponent<PixelT, PrecisionT> > {
    public:

      template<class ImageT>
      AffineMixtureComponent(ImageViewBase<ImageT> const& left,
                             ImageViewBase<ImageT> const& right,
                             Vector2i window_size, PrecisionT sigma_0, PrecisionT sigma_min) :
      left_r(left.impl()), right_r(right.impl()), s(sigma_0), s0(sigma_0), s_min(sigma_min),
        M_transform_linear(matrix_data_linear), M_transform_offset(matrix_data_offset),
        pos_linearize(0, 0), w(window_size(0), window_size(1)), p(window_size(0), window_size(1)), w_gaussian(window_size(0), window_size(1)),
        r_window(window_size(0), window_size(1)), l_window(window_size(0), window_size(1)), err(window_size(0), window_size(1))
        {
          matrix_data_linear[0] = 1.; matrix_data_linear[1] = 0.; // linear component of affine transform
          matrix_data_linear[2] = 0; matrix_data_linear[3] = 1; // linear component of affine transform
          matrix_data_offset[0] = 0;
          matrix_data_offset[1] = 0;
          T = AffineTransformOrigin(M_transform_linear, M_transform_offset, pos_linearize);

          l_set = false;

          sqrt_2_pi = sqrt(2*M_PI);

          fill(w, 1.);
          fill(p, .5);

          // make derivative kernels
          ImageView<PrecisionT> dx_kernel(3,3), dy_kernel(3,3);
          fill(dx_kernel, 0.);
          dx_kernel(0,1) = .5; dx_kernel(1,1) = 0; dx_kernel(2,1) = -.5;
          fill(dy_kernel, 0.);
          dy_kernel(1,0) = .5; dy_kernel(1,1) = 0; dy_kernel(1,2) = -.5;

          // precompute and cache right image derivatives
          r_image_dx = convolution_filter(right_r, dx_kernel);
          r_image_dy = convolution_filter(right_r, dy_kernel);

          l_likelihood = INFINITY;

          debug = false;

          // numerical optimization parameter default values
          epsilon_inner = 1e-8;
          inner_loop_iter_max = 8; //20;
          min_det = .1; // minimum determinant of the affine transformation; this determines the ratio that the affine transform is allowed to skew the area of the window by.
          max_det =  1.9; // maximum determinant of the affine transformation; this determines the ratio that the affine transform is allowed to skew the area of the window by.

          window.min() = Vector2(0, 0);
          window.max() = Vector2(window_size(0), window_size(1));
        }

      void reset(BBox2i left_window,
                 PrecisionT disparity_x, PrecisionT disparity_y,
                 Matrix2x2 warp, ImageView<PrecisionT> const& g) {
        fill(w, 1.);
        fill(p, .5);
        s = s0;


        //ut << "reseting warp to: " << warp << std::endl;

        l_set = false;

        window = left_window;
        pos_linearize = window.center();

        matrix_data_linear[0] = warp(0, 0); matrix_data_linear[1] = warp(0, 1); // linear component of affine transform
        matrix_data_linear[2] = warp(1, 0); matrix_data_linear[3] = warp(1, 1); // linear component of affine transform
        matrix_data_offset[0] = disparity_x;
        matrix_data_offset[1] = disparity_y;

        T = AffineTransformOrigin(M_transform_linear, M_transform_offset, pos_linearize);

        r_window = crop(transform(right_r, T), window);
        l_window = crop(left_r, window);
        err = r_window - l_window;

        w_gaussian = g;
      }

      PrecisionT log_likelihood() {
        if(!l_set) {
          PrecisionT sum_weights = sum_of_pixel_values(w); // this is the normalization constant for the weights
          l_likelihood = sum_of_pixel_values(.5*w*pow(err/s,2) - w*log(sqrt(w_gaussian))) + sum_weights*(log(s) + log(sqrt_2_pi));
          l_set = true;
        }
        return l_likelihood;
      }

      void update_posterior() {
        // our model distribution is N(J_{Tx}, s/w_gaussian)
        p = exp(-.5*w_gaussian*pow(err/s, 2))*sqrt(w_gaussian)/(s*sqrt_2_pi); // gaussian posterior density
        // <=> to exp(-.5*pow(err/(s/sqrt(w_g)), 2))/(sqrt_2_pi*s/sqrt(w_g)); // gaussian posterior density
      }

      PrecisionT sigma() {
        return s;
      }

      void print_status(std::string const& prefix) {
        std::cout << prefix << "sigma = " << s << std::endl;
      }

      Matrix<PrecisionT, 6, 6> const& hessian() const { return hess; }


      MatrixProxy<PrecisionT, 2, 2> const& affine_transform_mat() {
        return M_transform_linear;
      }

      VectorProxy<PrecisionT, 2> const& affine_transform_offset() {
        return M_transform_offset;
      }

      AffineTransformOrigin const& affine_transform() const { return T; }

      ImageView<PrecisionT> const& errors() { return err; }

      ImageView<PrecisionT> & weights() { return w; }

      ImageView<PrecisionT> const& weights() const { return w; }

      ImageView<PrecisionT> const& prob() const { return p; }

      void set_debug(bool p_debug) { debug = p_debug; }

      inline void m_compute_gradient_hessian(Vector<PrecisionT, 6> &gradient, Matrix<double, 6, 6> &hessian,
                                             int x, int y,
                                             ImageView<PrecisionT> const &weights, ImageView<PixelT> const &err,
                                             ImageView<PixelT> const &r_window_dx, ImageView<PixelT> const &r_window_dy) const
      {
        int kernel_height = err.rows();
        int kernel_width = err.cols();

        typename ImageView<PrecisionT>::pixel_accessor weights_iter = weights.origin();
        typename ImageView<PixelT>::pixel_accessor err_iter = err.origin();

        typename ImageView<PixelT>::pixel_accessor dx_iter = r_window_dx.origin();
        typename ImageView<PixelT>::pixel_accessor dy_iter = r_window_dy.origin();

        fill(gradient, 0);
        hessian.set_zero();
        int i, j;
        for(i = 0; i < kernel_height; i++) {
          for(j = 0; j < kernel_width; j++) {
            PrecisionT ty = y + i - (kernel_height-1)/2.;
            PrecisionT tx = x + j - (kernel_width-1)/2.;

            PrecisionT txtx = tx*tx;
            PrecisionT tyty = ty*ty;
            PrecisionT txty = tx*ty;


            PrecisionT cache_w = *weights_iter;
            PrecisionT cache_e = *err_iter;
            PrecisionT cache_dx = *dx_iter;
            PrecisionT cache_dy = *dy_iter;

            // these are implemented directly below, but commented out here as a reference
            //PrecisionT cache_a = cache_dx*cache_dx;
            //PrecisionT cache_b = cache_dx*cache_dy;
            //PrecisionT cache_c = cache_dy*cache_dy;

            PrecisionT cache_w_cache_e = cache_w*cache_e;
            PrecisionT cache_w_cache_a = cache_w*cache_dx*cache_dx;
            PrecisionT cache_w_cache_b = cache_w*cache_dx*cache_dy;
            PrecisionT cache_w_cache_c = cache_w*cache_dy*cache_dy;

            gradient(0) += cache_w_cache_e*tx*cache_dx; // lhs(0) += tx * dx_val * de_val; // a
            gradient(1) += cache_w_cache_e*tx*cache_dy; // lhs(3) += tx * dy_val * e_val;  // c
            gradient(2) += cache_w_cache_e*ty*cache_dx; // lhs(1) += ty * dx_val * e_val;  // b
            gradient(3) += cache_w_cache_e*ty*cache_dy; // lhs(4) += ty * dy_val * e_val;  // d
            gradient(4) += cache_w_cache_e*cache_dx;    // lhs(2) +=      dx_val * e_val;  // horizontal
            gradient(5) += cache_w_cache_e*cache_dy;    // lhs(5) +=      dy_val * e_val;  // vertical

            hessian(0, 0) += cache_w_cache_a*txtx;
            hessian(0, 1) += cache_w_cache_b*txtx;
            hessian(0, 2) += cache_w_cache_a*txty;

            hessian(0, 3) += cache_w_cache_b*txty;
            hessian(0, 4) += cache_w_cache_a*tx;
            hessian(0, 5) += cache_w_cache_b*tx;

            hessian(1, 1) += cache_w_cache_c*txtx;
            hessian(1, 2) += cache_w_cache_b*txty;

            hessian(1, 3) += cache_w_cache_c*txty;
            hessian(1, 4) += cache_w_cache_b*tx;
            hessian(1, 5) += cache_w_cache_c*tx;

            hessian(2, 2) += cache_w_cache_a*tyty;
            hessian(2, 3) += cache_w_cache_b*tyty;
            hessian(2, 4) += cache_w_cache_a*ty;
            hessian(2, 5) += cache_w_cache_b*ty;

            hessian(3, 3) += cache_w_cache_c*tyty;
            hessian(3, 4) += cache_w_cache_b*ty;
            hessian(3, 5) += cache_w_cache_c*ty;

            hessian(4, 4) += cache_w_cache_a;
            hessian(4, 5) += cache_w_cache_b;

            hessian(5, 5) += cache_w_cache_c;

            weights_iter.next_col();
            err_iter.next_col();
            dx_iter.next_col();
            dy_iter.next_col();
          }

        }
        int k, l;
        for(k = 0; k < 6; k++) {
          for(l = k; l < 6; l++) {
            hessian(l, k) = hessian(k, l);
          }
        }
      }

      void set_max_iter(int max_iter) {
        inner_loop_iter_max = max_iter;
      }
      void set_epsilon_convergence(PrecisionT epsilon) {
        epsilon_inner = epsilon;
      }
      void set_min_determinant(PrecisionT p_min_det) {
        min_det = p_min_det;
      }
      void set_max_determinant(PrecisionT p_max_det) {
        max_det = p_max_det;
      }

      void fit_parameters() {
        PrecisionT x = pos_linearize(0);
        PrecisionT y = pos_linearize(1);

        ImageViewRef<PrecisionT> w_adj = w*w_gaussian;

        ImageView<PixelT> r_window_dx(window.width(), window.height());
        ImageView<PixelT> r_window_dy(window.width(), window.height());

        Vector<PrecisionT, 6> gradient;
        Vector<PrecisionT, 6> soln;
        Matrix<PrecisionT, 6, 6> hessian_temp;


        Vector<PrecisionT, 6> gradient_hist[inner_loop_iter_max];
        Vector<PrecisionT, 6> state_hist[inner_loop_iter_max];
        Matrix<PrecisionT, 6, 6> hessian_hist[inner_loop_iter_max];

        PrecisionT initial_norm_grad = 0;
        if(debug) {
          std::cout << "initializing..." << std::endl;
        }
        // initialized cropped transformations of the right window under the current value of T; and err
        T = AffineTransformOrigin(M_transform_linear, M_transform_offset, Vector2(x, y));
        r_window = crop(transform(right_r, T), window);
        l_window = crop(left_r, window);
        err = r_window - l_window;

        double l_likelihood_old = sum_of_pixel_values(.5*w_adj*pow(err/s,2) - w*log(sqrt(w_gaussian))) + sum_of_pixel_values(w)*(log(s) + log(sqrt_2_pi)); //TODO: remove

        // the transformation parameter, 'T', is estimated iteratively in this inner loop
        int inner_loop_iter;
        PrecisionT f_value = INFINITY;
        Stopwatch inner_loop_timer;
        inner_loop_timer.start();

        for(inner_loop_iter = 0; inner_loop_iter < inner_loop_iter_max; inner_loop_iter++) {
          Stopwatch inner_loop_window_timer;
          inner_loop_window_timer.start();

          // compute the transformed derivatives
          r_window_dx = crop(transform(r_image_dx, T), window);
          r_window_dy = crop(transform(r_image_dy, T), window);

          err = r_window - l_window;

          inner_loop_window_timer.stop();
          //std::cout << "inner loop window setup " << 1000*inner_loop_window_timer.elapsed_seconds() << std::endl;


          Stopwatch gradient_timer;
          gradient_timer.start();


          if(debug) {
            std::cout << "computing gradient" << std::endl;
          }


          m_compute_gradient_hessian(gradient, hess, 0, 0,
                                     w_adj, err,
                                     r_window_dx, r_window_dy);

          if(inner_loop_iter == 0) {
            initial_norm_grad = norm_2(gradient);
          }

          if(debug) {
            hessian_hist[inner_loop_iter] = hess;

            state_hist[inner_loop_iter][0] = matrix_data_linear[0];
            state_hist[inner_loop_iter][1] = matrix_data_linear[2];
            state_hist[inner_loop_iter][2] = matrix_data_linear[1];
            state_hist[inner_loop_iter][3] = matrix_data_linear[3];

            state_hist[inner_loop_iter][4] = matrix_data_offset[0];
            state_hist[inner_loop_iter][5] = matrix_data_offset[1];
          }
          //std::cout << "gradient = " << gradient << std::endl;

          //std::cout << "||gradient|| = " << norm_2(gradient) << std::endl;

          //std::cout << "hessian computed as " << hess << std::endl;
          gradient_timer.stop();

          //std::cout << "gradient/hessian took " << 1000*gradient_timer.elapsed_seconds() << std::endl;

          Stopwatch lm_timer;
          // perform the update here; do line search
          PrecisionT f_value_last = sum_of_pixel_values(w_adj*pow(err, 2));
          bool min_step_size_hit = false;
          lm_timer.start();



          if(debug) {
            std::cout << "solving..." << std::endl;
          }
          soln = gradient;
          hessian_temp = hess;
          try {
            solve_symmetric_nocopy(hessian_temp, soln);

            for(int i = 0; i < 6; i++) {
              if(soln[i] != soln[i]) { // test for NaNs
                vw_out() << "NaN in soln!" << std::endl;
                soln = gradient;
                break;
              }
            }

          }
          catch(vw::Exception &e) {
            soln = gradient;
            //std::cout << "using gradient " << std::endl;
          }

          if(debug) {
            gradient_hist[inner_loop_iter] = soln;
          }

          if(debug) {
            std::cout << "line searching..." << std::endl;
          }
          PrecisionT alpha = .1;
          //std::cout << "f_value_last = " << f_value_last << std::endl;
          //std::cout << M_transform_linear << std::endl;
          PrecisionT matrix_data_linear_0[4];
          PrecisionT matrix_data_offset_0[2];

          while(true) {
            if(debug) {
              std::cout << "soln = " << (soln) << std::endl;
            }
            //std::cout << "M: " << M_transform_linear << std::endl;
            if(norm_2_sqr(soln) < epsilon_inner*epsilon_inner) { // if we've gotten to a step size this small, there's no more improvement to be made
              //std::cout << "updating" << std::endl;
              //PrecisionT determinant = matrix_data_linear[0]*matrix_data_linear[3] - matrix_data_linear[1]*matrix_data_linear[2]; // unused
              T = AffineTransformOrigin(M_transform_linear, M_transform_offset, Vector2(x, y));
              r_window = crop(transform(right_r, T), window);
              err = r_window - l_window;
              min_step_size_hit = true;
              //std::cout << "min_step_size_hit" << std::endl;
              break;
            }

            //std::cout << "trying step" << std::endl;

            // back up the current parameters
            matrix_data_linear_0[0] = matrix_data_linear[0];
            matrix_data_linear_0[1] = matrix_data_linear[1];
            matrix_data_linear_0[2] = matrix_data_linear[2];
            matrix_data_linear_0[3] = matrix_data_linear[3];
            matrix_data_offset_0[0] = matrix_data_offset[0];
            matrix_data_offset_0[1] = matrix_data_offset[1];

            //make a step
            // the weird indexing shift below is because soln is vectorized in column-major order, where as VW stores things row-major
            matrix_data_linear[0] -= soln[0];
            matrix_data_linear[1] -= soln[2];
            matrix_data_linear[2] -= soln[1];
            matrix_data_linear[3] -= soln[3];
            matrix_data_offset[0] -= soln[4];
            matrix_data_offset[1] -= soln[5];

            // make sure determinant is positive (there no flips involved)
            //std::cout << det(M_transform_linear) << std::endl;
            PrecisionT determinant = matrix_data_linear[0]*matrix_data_linear[3] - matrix_data_linear[1]*matrix_data_linear[2];
            if(debug) {
              std::cout << "det = " << determinant << std::endl;
            }
            if(determinant < min_det || determinant > max_det) {
              // restore last good parameters
              matrix_data_linear[0] = matrix_data_linear_0[0];
              matrix_data_linear[1] = matrix_data_linear_0[1];
              matrix_data_linear[2] = matrix_data_linear_0[2];
              matrix_data_linear[3] = matrix_data_linear_0[3];
              matrix_data_offset[0] = matrix_data_offset_0[0];
              matrix_data_offset[1] = matrix_data_offset_0[1];

              //std::cout << "min_det hit" << std::endl;
              //std::cout << soln << std::endl;
              soln *= alpha;
              continue;
            }

            // update T, the right window view, and the err
            T = AffineTransformOrigin(M_transform_linear, M_transform_offset, Vector2(x, y));

            //std::cout << "updating values" << std::endl;
            r_window = crop(transform(right_r, T), window);
            err = r_window - l_window;
            f_value = sum_of_pixel_values(w_adj*pow(err, 2));
            //std::cout << "solved" << std::endl;
            //std::cout << "testing: " << matrix_data_linear[0] << " " << matrix_data_linear[1] << " " << matrix_data_linear[2] << " " << matrix_data_linear[3] << " "
            // << matrix_data_offset[0] << " " << matrix_data_offset[1] << std::endl;

            //std::cout << "\t w/ f_value = " << f_value << std::endl;
            if(f_value < f_value_last) { // if there has been an improvement
              if(debug) {
                std::cout << "improved: " << f_value_last << " -> " << f_value << std::endl;
              }
              break;
            }
            else { // there is no improvement
              f_value = f_value_last;
              //std::cout << "backtracking" << std::endl;
              // restore last good parameters
              matrix_data_linear[0] = matrix_data_linear_0[0];
              matrix_data_linear[1] = matrix_data_linear_0[1];
              matrix_data_linear[2] = matrix_data_linear_0[2];
              matrix_data_linear[3] = matrix_data_linear_0[3];
              matrix_data_offset[0] = matrix_data_offset_0[0];
              matrix_data_offset[1] = matrix_data_offset_0[1];

              soln *= alpha;
              continue;
            }
          }
          //std::cout << "f_value = " << f_value << std::endl;

          lm_timer.stop();
          //std::cout << "linesearch took " << 1000*lm_timer.elapsed_seconds() << std::endl;

          // check termination condition for inner loop
          //std::cout << "delta is " << fabs(f_value_last - f_value) << std::endl;
          if(min_step_size_hit || fabs(f_value_last - f_value) < epsilon_inner) {// if we could not improve anything after degrading to gradient descent, or improvement is too small
            if(debug) {
              std::cout << "converged to " << f_value << std::endl;
              std::cout << "\tw/ ||g||/||g0|| = " << norm_2(gradient)/initial_norm_grad << std::endl;
            }

            break;
          }
        }

        // recompute these with the final converged values
        T = AffineTransformOrigin(M_transform_linear, M_transform_offset, Vector2(x, y));
        r_window = crop(transform(right_r, T), window);
        err = r_window - l_window;

        PrecisionT sum_weights = sum_of_pixel_values(w); // this is the normalization constant for the weights

        s = sqrt(sum_of_pixel_values(w_adj*pow(err, 2))/sum_weights);
        s = std::max<double>(s, s_min); // threshold the variance so it never gets too small

        l_likelihood = sum_of_pixel_values(.5*w_adj*pow(err/s,2) - w*log(sqrt(w_gaussian))) + sum_weights*(log(s) + log(sqrt_2_pi));
        l_set = true;

        if(l_likelihood_old + 1e-8 < l_likelihood) {
          std::cout << "LIKELIHOOD WENT UP IN Affine by "
                    <<  l_likelihood-l_likelihood_old << std::endl;
        }

        hess = hess/(sum_weights*s*s);
      }

    private:
      bool debug;

      ImageViewRef<PixelT> left_r;
      ImageViewRef<PixelT> right_r;

      PrecisionT s, s0;
      PrecisionT s_min;
      PrecisionT sqrt_2_pi;
      PrecisionT l_likelihood;

      ImageView<PixelT> r_image_dx, r_image_dy;

      Matrix<PrecisionT, 6, 6> hess;
      PrecisionT matrix_data_linear[4];
      PrecisionT matrix_data_offset[2];

      MatrixProxy<PrecisionT, 2, 2> M_transform_linear;
      VectorProxy<PrecisionT, 2> M_transform_offset;

      BBox2i window;
      Vector2 pos_linearize; // linearization point of the minimizer

      ImageView<PrecisionT> w;  // weights
      ImageView<PrecisionT> p; // posterior probabilities of data in window
      ImageView<PrecisionT> w_gaussian;  // gaussian window

      ImageView<PixelT> r_window;
      ImageView<PixelT> l_window;
      ImageView<PrecisionT> err;

      AffineTransformOrigin T;
      bool l_set;

      // numerical optimization parameters
      PrecisionT epsilon_inner;
      int inner_loop_iter_max;
      PrecisionT min_det;
      PrecisionT max_det;
    };


  }
}

#endif//__VW_STEREO_AFFINE_MIXTURE_COMPONENT__
