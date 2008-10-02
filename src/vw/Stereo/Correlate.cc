// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__
#include <vw/Stereo/Correlate.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/Interpolation.h>

//#include <vDSP.h>

namespace vw {
namespace stereo {  

  /// This routine cross checks L2R and R2L, placing the final version
  /// of the disparity map in L2R.
  void cross_corr_consistency_check(ImageView<PixelDisparity<float> > &L2R, 
                                    ImageView<PixelDisparity<float> > const& R2L,
                                    double cross_corr_threshold, bool verbose) {
    int32 xx,yy;
    int count = 0, match_count = 0;
  
    if (verbose)
      vw_out(InfoMessage, "stereo") << "\tCrosscorr threshold: " << cross_corr_threshold << "\n";
    if (cross_corr_threshold < 0) 
      vw_throw( vw::ArgumentErr() << "CrossCorrConsistencyCheck2D: the crosscorr threshold was less than 0." );
  
    for(xx = 0; xx < L2R.cols(); xx++) {     
      for(yy = 0; yy < L2R.rows(); yy++) {
      
        int xOffset = (int)L2R(xx,yy).h();
        int yOffset = (int)L2R(xx,yy).v();
      
        // Check to make sure we are within the image bounds
        if(xx+xOffset < 0 || yy+yOffset < 0 ||
           xx+xOffset >= R2L.cols() || yy+yOffset >= R2L.rows()) {
          L2R(xx,yy) = PixelDisparity<float>();  // Default constructor is missing pixel.
        }

        // Check for missing pixels
        else if ( L2R(xx,yy).missing() ||
                  R2L(xx+xOffset, yy+yOffset).missing() ) {
          L2R(xx,yy) = PixelDisparity<float>();  // Default constructor is missing pixel.
        }
      
        // Check for correlation consistency
        //
        // Since the hdisp for the R2L and L2R buffers will be opposite 
        // in sign, we determine their similarity by *summing* them, rather
        // than differencing them as you might expect.
        else if (cross_corr_threshold >= fabs(L2R(xx,yy).h() + R2L(xx+xOffset,yy+yOffset).h()) &&
                 cross_corr_threshold >= fabs(L2R(xx,yy).v() + R2L(xx+xOffset,yy+yOffset).v())) {
          count++;
          match_count++;
        }
      
        // Otherwise, the pixel is bad.
        else {
          match_count++;
          L2R(xx,yy) = PixelDisparity<float>();  // Default constructor is missing pixel.
        } 
      }
    } 
    if (verbose) 
      vw_out(InfoMessage, "stereo") << "\tCross-correlation retained " << count << " / " << match_count << " matches (" << ((float)count/match_count*100) <<" percent).\n";
  }

  static double find_minimum(double lt, double mid, double rt) {
    double a = (rt+lt)*0.5-mid;
    double b = (rt-lt)*0.5;
    return -b/(2.0*a);
  }

  /* 
   * Find the minimun of a 2d hyperbolic surface that is fit to the nine points 
   * around and including the peak in the disparity map.  This gives better 
   * subpixel resolution when both horizontal and vertical subpixel is requested.
   * 
   * The equation of the surface we are fitting is:
   *    z = ax^2 + by^2 + cxy + dx + ey + f
   */
  template <class VectorT, class MatrixT>
  static vw::Vector2 find_minimum_2d(vw::VectorBase<VectorT> &points, vw::MatrixBase<MatrixT> &pinvA) {

    vw::Vector2 offset;

    /* 
     * First, compute the parameters of the hyperbolic surface by fitting the nine points in 'points'
     * using a linear least squares fit.  This process is fairly fast, since we have already pre-computed
     * the inverse of the A matrix in Ax = b.
     */
    vw::Vector<double> x = pinvA * points;
  
    /* 
     * With these parameters, we have a closed form expression for the surface.  We compute the 
     * derivative, and find the point where the slope is zero.  This is our maximum.
     *
     * Max is at [x,y] where:
     *
     *   dz/dx = 2ax + cy + d = 0
     *   dz/dy = 2by + cx + e = 0
     * 
     * Of course, we optimize this computation a bit by unrolling it by hand beforehand.
     */
    double denom = 4 * x(0) * x(1) - (x(2) * x(2));
  
    offset(0) = ( x(2) * x(4) - 2 * x(1) * x(3) ) / denom;
    offset(1) = ( x(2) * x(3) - 2 * x(0) * x(4) ) / denom;

    return offset;
  }

  ///-------------------------------------------------------------------------

  inline ImageView<float> compute_gaussian_weight_image(int kern_width, int kern_height) {

    int center_pix_x = kern_width/2;
    int center_pix_y = kern_height/2;
    int two_sigma_sqr = 2*pow(kern_width/3,2);

    ImageView<float> weight(kern_width, kern_height);
    for (int j = 0; j < kern_height; ++j) {
      for (int i = 0; i < kern_width; ++i ) {
        weight(i,j) = exp(-1 * (pow(i-center_pix_x,2) + pow(j-center_pix_y,2)) / two_sigma_sqr);
      }
    }
    return weight;
  }

  inline int adjust_weight_image(ImageView<float> &weight,
                                 ImageView<PixelDisparity<float> > const& disparity_map_patch,
                                 ImageView<float> const& weight_template) {
    
    //    const float continuity_threshold_squared = 64;  // T = 8
    int center_pix_x = weight_template.cols()/2;
    int center_pix_y = weight_template.rows()/2;
    PixelDisparity<float> center_pix = disparity_map_patch(center_pix_x, center_pix_y);

    float sum = 0;
    int num_good_pix = 0;
    ImageView<float>::pixel_accessor weight_row_acc = weight.origin();
    ImageView<float>::pixel_accessor template_row_acc = weight_template.origin();
    ImageView<PixelDisparity<float> >::pixel_accessor disp_row_acc = disparity_map_patch.origin();
    for (int j = 0; j < weight_template.rows(); ++j) {
      ImageView<float>::pixel_accessor weight_col_acc = weight_row_acc;
      ImageView<float>::pixel_accessor template_col_acc = template_row_acc;
      ImageView<PixelDisparity<float> >::pixel_accessor disp_col_acc = disp_row_acc;
      for (int i = 0; i < weight_template.cols(); ++i ) {

        // Mask is zero if the disparity map's pixel is missing...
        if ( (*disp_col_acc).missing()) 
          *weight_col_acc = 0;

//         // ... or if there is a large discontinuity ...
//         if (pow( (*disp_col_acc).h()-center_pix.h(),2) + pow( (*disp_col_acc).v()-center_pix.v(),2) >= continuity_threshold_squared)
//           *weight_col_acc = 0;

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
      vw_throw(LogicErr() << "subpixel_weight: Sum of weight image was zero.  This isn't supposed to happen!");
    else 
      weight /= sum;
    return num_good_pix;
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

  template<class ChannelT>
  void subpixel_correlation_affine_2d(ImageView<PixelDisparity<float> > &disparity_map,
                                      ImageView<ChannelT> const& left_image,
                                      ImageView<ChannelT> const& right_image,
                                      int kern_width, int kern_height,
                                      bool do_horizontal_subpixel,
                                      bool do_vertical_subpixel,
                                      bool verbose) {
    
    VW_ASSERT( disparity_map.cols() == left_image.cols() &&
               disparity_map.rows() == left_image.rows(),
               ArgumentErr() << "subpixel_correlation: left image and disparity map do not have the same dimensions.");

    // Robust cost function settings
    const float thresh = 0.01;
    HuberError robust_cost_fn(thresh);

    int kern_pixels = kern_height * kern_width;
    int weight_threshold = kern_pixels / 4;

    // Bail out if no subpixel computation has been requested 
    if (!do_horizontal_subpixel && !do_vertical_subpixel) return;

    ImageView<float> x_deriv = derivative_filter(left_image, 1, 0);
    ImageView<float> y_deriv = derivative_filter(left_image, 0, 1);
    ImageView<float> weight_template = compute_gaussian_weight_image(kern_width, kern_height);

    // Workspace images are allocated up here out of the tight inner
    // loop.  We rasterize into these directly in the code below.
    ImageView<float> w(kern_width, kern_height);
    
    Vector2 offset;

    // Iterate over all of the pixels in the disparity map except for
    // the outer edges.
    Stopwatch sw;
    sw.start();
    double last_time = 0;
    for (int y=kern_height/2; y<left_image.rows()-kern_height/2; ++y) {
      if (verbose && y % 10 == 0) {
        sw.stop();
        vw_out(InfoMessage, "stereo") << "\tProcessing subpixel line: " << y << " / " << left_image.rows() << "    (" << (10 * left_image.cols() / (sw.elapsed_seconds() - last_time)) << " pixels/s, "<< sw.elapsed_seconds() << " s total )      \r" << std::flush;
        last_time = sw.elapsed_seconds();
        sw.start();
      }
      for (int x=kern_width/2; x<left_image.cols()-kern_width/2; ++x) {
        BBox2i current_window(x-kern_width/2, y-kern_height/2, kern_width, kern_height);
        Vector2 base_offset( -disparity_map(x,y).h() , -disparity_map(x,y).v() );          

        // Skip over pixels for which we have no initial disparity estimate
        if (disparity_map(x,y).missing())
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
          disparity_map(x,y) = PixelDisparity<float>();
          continue;
        }
                
        // Iterate until a solution is found or the max number of
        // iterations is reached.
        for (unsigned iter = 0; iter < 10; ++iter) {
          offset(0) = d[2];
          offset(1) = d[5];
          
          // First we check to see if our current subpixel translation
          // is less than one half of the window width.  If not, then
          // we are probably having trouble converging and we abort
          // this pixel!!
          if (norm_2(offset) > kern_width/2) 
            break;
          
          InterpolationView<EdgeExtensionView<ImageView<ChannelT>, ZeroEdgeExtension>, BilinearInterpolation> right_interp_image =
            interpolate(right_image, BilinearInterpolation(), ZeroEdgeExtension());
          
          float x_base = x + disparity_map(x,y).h();
          float y_base = y + disparity_map(x,y).v();
          //          float error_total = 0;

          Matrix<float,6,6> rhs;
          Vector<float,6> lhs;
          for (int jj = -kern_height/2; jj <= kern_height/2; ++jj) {
            for (int ii = -kern_width/2; ii <= kern_width/2; ++ii) {
              int i = ii + kern_width/2;
              int j = jj + kern_height/2;

              // First we compute the pixel offset for the right image
              // and the error for the current pixel.
              float xx = x_base + d[0] * ii + d[1] * jj + offset(0);
              float yy = y_base + d[3] * ii + d[4] * jj + offset(1);
              float I_e_val = right_interp_image(xx,yy) - left_image_patch(i,j) + 1e-16; 

              // Apply the robust cost function.  We use a huber
              // function to gently remove outliers for small errors,
              // but we set a hard limit a 5 times the cost threshold
              // to remove major (salt&pepper) noise.
              //              if (fabs(I_e_val) > thresh*5 || I_e_val == 0.0)
              float robust_weight = sqrt(robust_cost_fn(fabs(I_e_val)))/fabs(I_e_val);
              //              error_total += pow(I_e_val,2);

              // We combine the error value with the derivative and
              // add this to the update equation.
              float I_x_val = robust_weight * w(i,j) * I_x(i,j);
              float I_y_val = robust_weight * w(i,j) * I_y(i,j);
              float I_x_sqr = robust_weight * w(i,j) * I_x(i,j) * I_x(i,j);
              float I_y_sqr = robust_weight * w(i,j) * I_y(i,j) * I_y(i,j);
              float I_x_I_y = robust_weight * w(i,j) * I_x(i,j) * I_y(i,j);

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
            }
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

          //           {          
          //             ImageView<ChannelT> right_image_patch(kern_width, kern_height);
          //             for (int jj = -kern_height/2; jj <= kern_height/2; ++jj) {
          //               for (int ii = -kern_width/2; ii <= kern_width/2; ++ii) {
          //                 float xx = x_base + d[0] * ii + d[1] * jj + offset(0);
          //                 float yy = y_base + d[3] * ii + d[4] * jj + offset(1);
          //                 right_image_patch(ii+kern_width/2, jj+kern_width/2) = right_interp_image(xx,yy);
          //               }
          //             }
          //             std::ostringstream ostr;
          //             ostr << x << "_" << y << "-" << iter;
          //             write_image("small/left-"+ostr.str()+".tif", left_image_patch);
          //             write_image("small/right-"+ostr.str()+".tif", right_image_patch);
          //             write_image("small/weight-"+ostr.str()+".tif", w);
          //           }


          // Solves lhs = rhs * x, and stores the result in-place in lhs.
          Matrix<double,6,6> pre_rhs = rhs;
          Vector<double,6> pre_lhs = lhs;
          try { solve_symmetric_nocopy(rhs,lhs);
          } catch (ArgumentErr &e) {
            std::cout << "Error @ " << x << " " << y << "\n";
//             std::cout << "Exception caught: " << e.what() << "\n";
//             std::cout << "PRERHS: " << pre_rhs << "\n";
//             std::cout << "PRELHS: " << pre_lhs << "\n\n";
//             std::cout << "RHS: " << rhs << "\n";
//             std::cout << "LHS: " << lhs << "\n\n";
//             std::cout << "DEBUG: " << rhs(0,1) << "   " << rhs(1,0) << "\n\n";
//             exit(0);
          }
          d += lhs;

          //          std::cout << "Update: " << lhs << "     " << d << "     " << sqrt(error_total) << "    " << (sqrt(lhs[2]*lhs[2]+lhs[5]*lhs[5])) << "\n";

          // Termination condition
          if (norm_2(lhs) < 0.01) 
            break;
        }
        //        std::cout << "----> " << d << "\n\n";
        offset(0) = d[2];
        offset(1) = d[5];
        
        if (norm_2(offset) > 1.5 || 
            offset(0) != offset(0) ||  // Check to make sure the offset is not NaN...
            offset(1) != offset(1) ) { // ... ditto.
          disparity_map(x,y) = PixelDisparity<float>();
        } else {
          disparity_map(x,y).h() += offset(0);
          disparity_map(x,y).v() += offset(1);
        }
      }
    }
    if (verbose) 
      vw_out(InfoMessage, "stereo") << "\tProcessing subpixel line: done.                                         \n";
  }


  template<class ChannelT>
  void subpixel_correlation_parabola(ImageView<PixelDisparity<float> > &disparity_map,
                                        ImageView<ChannelT> const& left_image,
                                        ImageView<ChannelT> const& right_image,
                                        int kern_width, int kern_height,
                                        bool do_horizontal_subpixel,
                                        bool do_vertical_subpixel,
                                        bool verbose) {

    VW_ASSERT(left_image.cols() == right_image.cols() && left_image.cols() == disparity_map.cols() &&
              left_image.rows() == right_image.rows() && left_image.rows() == disparity_map.rows(),
              ArgumentErr() << "subpixel_correlation: input image dimensions do not agree.\n");

    int height = disparity_map.rows();
    int width = disparity_map.cols();

    ChannelT *new_img0 = &(left_image(0,0));
    ChannelT *new_img1 = &(right_image(0,0));
  
    // Bail out if no subpixel computation has been requested 
    if (!do_horizontal_subpixel && !do_vertical_subpixel) return;
  
    // We get a considerable speedup in our 2d subpixel correlation if
    // we go ahead and compute the pseudoinverse of the A matrix (where
    // each row in A is [ x^2 y^2 xy x y 1] (our 2d parabolic surface)
    // for the range of x = [-1:1] and y = [-1:1].
    static double pinvA_data[] = { 1.0/6,  1.0/6,  1.0/6, -1.0/3, -1.0/3, -1.0/3,  1.0/6,  1.0/6,  1.0/6,
                                   1.0/6, -1.0/3,  1.0/6,  1.0/6, -1.0/3,  1.0/6,  1.0/6, -1.0/3,  1.0/6,
                                   1.0/4,    0.0, -1.0/4,    0.0,    0.0,    0.0, -1.0/4,    0.0,  1.0/4,
                                   -1.0/6, -1.0/6, -1.0/6,    0.0,    0.0,   0.0,  1.0/6,  1.0/6,  1.0/6,
                                   -1.0/6,    0.0,  1.0/6, -1.0/6,    0.0, 1.0/6, -1.0/6,    0.0,  1.0/6,
                                   -1.0/9,  2.0/9, -1.0/9,  2.0/9,  5.0/9, 2.0/9, -1.0/9,  2.0/9, -1.0/9 }; 
    vw::MatrixProxy<double,6,9> pinvA(pinvA_data);
    for (int r = 0; r < height; r++) {
      if (r%100 == 0) 
        if (verbose) vw_out(InfoMessage, "stereo") << "\tPerforming sub-pixel correlation... "<< (double(r)/height * 100) << "%        \r" << std::flush;
    
      for (int c = 0; c < width; c++) {
      
        if ( !disparity_map(c,r).missing() ) {
          int hdisp= (int)disparity_map(c,r).h();
          int vdisp= (int)disparity_map(c,r).v();

          double mid = compute_soad(new_img0, new_img1,
                                    r, c,
                                    hdisp,   vdisp,
                                    kern_width, kern_height,
                                    width, height);
        
          // If only horizontal subpixel resolution is requested 
          if (do_horizontal_subpixel && !do_vertical_subpixel) {
            double lt= compute_soad(new_img0, new_img1,
                                    r, c,
                                    hdisp-1, vdisp,
                                    kern_width, kern_height,
                                    width, height);
            double rt= compute_soad(new_img0, new_img1,
                                    r, c,
                                    hdisp+1, vdisp,
                                    kern_width, kern_height,
                                    width, height);
          
            if ((mid <= lt && mid < rt) || (mid <= rt && mid < lt)) {
              disparity_map(c,r).h() += find_minimum(lt, mid, rt);
            } else {
              disparity_map(c,r) = PixelDisparity<float>();
            }
          }
      
          // If only vertical subpixel resolution is requested 
          if (do_vertical_subpixel && !do_horizontal_subpixel) {
            double up= compute_soad(new_img0, new_img1,
                                    r, c,
                                    hdisp, vdisp-1,
                                    kern_width, kern_height,
                                    width, height);
            double dn= compute_soad(new_img0, new_img1,
                                    r, c,
                                    hdisp, vdisp+1,
                                    kern_width, kern_height,
                                    width, height);
          
            if ((mid <= up && mid < dn) || (mid <= dn && mid < up)) {
              disparity_map(c,r).v() += find_minimum(up, mid, dn);
            } else {
              disparity_map(c,r) = PixelDisparity<float>();
            }
          }
        
        
          // If both vertical and horizontal subpixel resolution is requested,
          // we try to fit a 2d hyperbolic surface using the 9 points surrounding the
          // peak SOAD value.  
          //
          // We place the soad values into a vector using the following indices
          // (i.e. index 4 is the max disparity value)
          // 
          //     0  3  6
          //     1  4  7
          //     2  5  8
          //
          if (do_vertical_subpixel && do_horizontal_subpixel) {
            vw::Vector<double,9> points;
          
            points(0) = (double)compute_soad(new_img0, new_img1,
                                             r, c,
                                             hdisp-1, vdisp-1,
                                             kern_width, kern_height,
                                             width, height);
            points(1) = (double)compute_soad(new_img0, new_img1,
                                             r, c,
                                             hdisp-1, vdisp,
                                             kern_width, kern_height,
                                             width, height);
            points(2) = (double)compute_soad(new_img0, new_img1,
                                             r, c,
                                             hdisp-1, vdisp+1,
                                             kern_width, kern_height,
                                             width, height);
            points(3) = (double)compute_soad(new_img0, new_img1,
                                             r, c,
                                             hdisp, vdisp-1,
                                             kern_width, kern_height,
                                             width, height);
            points(4) = (double)mid;
            points(5) = (double)compute_soad(new_img0, new_img1,
                                             r, c,
                                             hdisp, vdisp+1,
                                             kern_width, kern_height,
                                             width, height);
            points(6) = (double)compute_soad(new_img0, new_img1,
                                             r, c,
                                             hdisp+1, vdisp-1,
                                             kern_width, kern_height,
                                             width, height);
            points(7) = (double)compute_soad(new_img0, new_img1,
                                             r, c,
                                             hdisp+1, vdisp,
                                             kern_width, kern_height,
                                             width, height);
            points(8) = (double)compute_soad(new_img0, new_img1,
                                             r, c,
                                             hdisp+1, vdisp+1,
                                             kern_width, kern_height,
                                             width, height);
            
            vw::Vector2 offset = find_minimum_2d(points, pinvA);
          
            // This prevents us from adding in large offsets for
            // poorly fit data.
            if (fabs(offset(0)) < 2.0 && fabs(offset(1)) < 2.0) {
              disparity_map(c,r).h() += offset(0);
              disparity_map(c,r).v() += offset(1);
            } else {
              disparity_map(c,r) = PixelDisparity<float>();
            }
          }
        } 
      }
    }
    if (verbose) vw_out(InfoMessage, "stereo") << "\tPerforming sub-pixel correlation... done.                 \n";
  }

  template void subpixel_correlation_affine_2d(ImageView<PixelDisparity<float> > &disparity_map,
                                     ImageView<uint8> const& left_image,
                                     ImageView<uint8> const& right_image,
                                     int kern_width, int kern_height,
                                     bool do_horizontal_subpixel,
                                     bool do_vertical_subpixel,
                                     bool verbose);

  template void subpixel_correlation_affine_2d(ImageView<PixelDisparity<float> > &disparity_map,
                                     ImageView<float> const& left_image,
                                     ImageView<float> const& right_image,
                                     int kern_width, int kern_height,
                                     bool do_horizontal_subpixel,
                                     bool do_vertical_subpixel,
                                     bool verbose);

  template void subpixel_correlation_parabola(ImageView<PixelDisparity<float> > &disparity_map,
                                     ImageView<uint8> const& left_image,
                                     ImageView<uint8> const& right_image,
                                     int kern_width, int kern_height,
                                     bool do_horizontal_subpixel,
                                     bool do_vertical_subpixel,
                                     bool verbose);

  template void subpixel_correlation_parabola(ImageView<PixelDisparity<float> > &disparity_map,
                                     ImageView<float> const& left_image,
                                     ImageView<float> const& right_image,
                                     int kern_width, int kern_height,
                                     bool do_horizontal_subpixel,
                                     bool do_vertical_subpixel,
                                     bool verbose);

}} // namespace vw::stereo
