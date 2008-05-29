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
#ifndef __VW_STEREO_CORRELATE_H__
#define __VW_STEREO_CORRELATE_H__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/LinearAlgebra.h>
#include <vw/Stereo/DisparityMap.h>

#define VW_STEREO_MISSING_PIXEL -32000

namespace vw {
namespace stereo {

VW_DEFINE_EXCEPTION(CorrelatorErr, vw::Exception);

  // Timing (profiling) utilites
  static inline double Time(void) {
    struct timeval tv;
    gettimeofday(&tv, 0);
    return ((double) tv.tv_sec + (double) tv.tv_usec / 1.0e6);
  }


  // Sign of the Laplacian of the Gaussian pre-processing
  // 
  // Default gaussian blur standard deviation is 1.5 pixels.
  class SlogStereoPreprocessingFilter {
    float m_slog_width;

  public:
    typedef ImageView<uint8> result_type;

    SlogStereoPreprocessingFilter(float slog_width = 1.5) : m_slog_width(slog_width) {} 

    template <class ViewT>
    result_type operator()(ImageViewBase<ViewT> const& view) const {
      return channel_cast<uint8>(threshold(laplacian_filter(gaussian_filter(channel_cast<float>(view.impl()),m_slog_width)),0.0));
    }

    static bool use_bit_image() { return true; }
  };
  
  // Laplacian of Gaussian pre-processing
  // 
  // Default gaussian blur standard deviation is 1.5 pixels.
  class LogStereoPreprocessingFilter {
    float m_log_width;
    
  public:
    typedef ImageView<float> result_type;

    LogStereoPreprocessingFilter(float log_width = 1.5) : m_log_width(log_width) {}
    
    template <class ViewT>
    result_type operator()(ImageViewBase<ViewT> const& view) const {
      return laplacian_filter(gaussian_filter(channel_cast<float>(view.impl()),m_log_width));
    }
    
    static bool use_bit_image() { return false; }
  };
  
  // No pre-processing
  class NullStereoPreprocessingFilter {
  public:
    typedef ImageView<uint8> result_type;

    template <class ViewT>
    result_type operator()(ImageViewBase<ViewT> const& view) const {
      return channel_cast_rescale<uint8>(view.impl());
    }
    static bool use_bit_image() { return false; }    
  };



  // Return absolute difference of two bit images
  inline unsigned char correlator_absolute_difference(unsigned char val1, unsigned char val2) {
    return val1 ^ val2;
  }
  
  // Return absolute difference of two bit images
  inline float correlator_absolute_difference(float val1, float val2) {
    float diff = val1 - val2;
    if (diff > 0) return diff; else return -diff;
  }


  /// Compute the sum of the absolute difference between a template
  /// region taken from img1 and the window centered at (c,r) in img0.
  template <class ChannelT>
  inline double compute_soad(ChannelT *img0, ChannelT *img1,
                             int r, int c,                   // row and column in img0
                             int hdisp, int vdisp,           // Current disparity offset from (c,r) for img1
                             int kern_width, int kern_height,// Kernel dimensions
                             int width, int height) {        // Image dimensions
  
    r -= kern_height/2;
    c -= kern_width/2;
    if (r<0         || c<0       || r+kern_height>=height       || c+kern_width>=width ||
        r+vdisp < 0 || c+hdisp<0 || r+vdisp+kern_height>=height || c+hdisp+kern_width>=width) {
      return VW_STEREO_MISSING_PIXEL;
    }

    ChannelT *new_img0 = img0;
    ChannelT *new_img1 = img1;

    new_img0 += c + r*width;
    new_img1 += (c+hdisp) + (r+vdisp)*width;
  
    typename vw::AccumulatorType<ChannelT>::type ret = 0;
    for (int rr= 0; rr< kern_height; rr++) {
      for (int cc= 0; cc< kern_width; cc++) {
        ret += correlator_absolute_difference(new_img0[cc], new_img1[cc]);
      }
      new_img0 += width;
      new_img1 += width;
    }
    return double(ret);
  }

  /// For a given set of images, compute the optimal disparity (minimum
  /// SOAD) at position left_image(i,j) for the given correlation window
  /// settings.
  /// 
  /// The left_image and right_image must have the same dimensions, but
  /// this is only checked here if debugging is enabled.
  template <class ChannelT>
  inline PixelDisparity<float> compute_disparity(ImageView<ChannelT> &left_image,
                                                 ImageView<ChannelT> &right_image,
                                                 int i, int j,
                                                 int kern_width, int kern_height,
                                                 int min_h_disp, int max_h_disp,
                                                 int min_v_disp, int max_v_disp) {

    const double default_soad = 1.0e10;     // Impossibly large value
    double min_soad = default_soad;
    PixelDisparity<float> best_disparity; // Starts as a missing pixel
    for (int ii = min_h_disp; ii <= max_h_disp; ++ii) {
      for (int jj = min_v_disp; jj <= max_v_disp; ++jj) {
        double soad = compute_soad(&(left_image(0,0)), &(right_image(0,0)),
                                   j, i, ii, jj,kern_width, kern_height, 
                                   left_image.cols(), left_image.rows());
        if (soad != VW_STEREO_MISSING_PIXEL && soad < min_soad) {
          min_soad = soad;
          best_disparity = PixelDisparity<float>(ii, jj);
        }
      }
    }
    return best_disparity;
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
  
//     offset(0) = ( x(2) * x(4) - 2 * x(1) * x(3) ) / denom;
//     offset(1) = ( x(2) * x(3) - 2 * x(1) * x(4) ) / denom;
  
    offset(0) = ( x(2) * x(4) - 2 * x(1) * x(3) ) / denom;
    offset(1) = ( x(2) * x(3) - 2 * x(0) * x(4) ) / denom;

    return offset;
  }

  template <class ChannelT> 
  void subpixel_correlation_new(ImageView<PixelDisparity<float> > &disparity_map,
                            ImageView<ChannelT> const& left_input_image,
                            ImageView<ChannelT> const& right_input_image,
                            int kern_width, int kern_height,
                            bool do_horizontal_subpixel = true,
                            bool do_vertical_subpixel = true,
                            bool verbose = false) {
    VW_ASSERT(left_input_image.cols() == right_input_image.cols() && left_input_image.cols() == disparity_map.cols() &&
              left_input_image.rows() == right_input_image.rows() && left_input_image.rows() == disparity_map.rows(),
              ArgumentErr() << "subpixel_correlation: input image dimensions do not agree.\n");

     // Bail out if no subpixel computation has been requested 
    if (!do_horizontal_subpixel && !do_vertical_subpixel) return;

    EdgeExtensionView<ImageView<ChannelT>, ZeroEdgeExtension > left_image = edge_extend(left_input_image, ZeroEdgeExtension());
    EdgeExtensionView<ImageView<ChannelT>, ZeroEdgeExtension > right_image = edge_extend(right_input_image, ZeroEdgeExtension());

    int kern_half_height = kern_height/2;
    int kern_half_width = kern_width/2;
    
    Vector<double> I_1((kern_width+1)*(kern_height+1));        
    Matrix<double> I_2((kern_width+1)*(kern_height+1), 10); 

    for (int j = 0; j < disparity_map.rows(); j++) {
      if (j%100 == 0) 
        if (verbose) vw_out(InfoMessage, "stereo") << "\tPerforming sub-pixel correlation... "<< (double(j)/disparity_map.rows() * 100) << "%              \r" << std::flush;
      for (int i = 0; i < disparity_map.cols(); i++) {
        if ( !disparity_map(i,j).missing() ) {

          // What is the starting (integer) disparity offset?
          int hdisp= (int)disparity_map(i,j).h();
          int vdisp= (int)disparity_map(i,j).v();

//           if (i == 237 && j == 237)
//             std::cout << "\n\nPixel 237x237: " << hdisp << " " << vdisp << "      ";

          // Copy the image data into I_1 and I_2
          int idx = 0;
          for (int r = -kern_half_height; r <= kern_half_height; ++r) {
            for (int c = -kern_half_width; c <= kern_half_width; ++c) { 
              I_1(idx) = left_image(i+c, j+r);
//               I_2(idx,0) = right_image(i+c+hdisp-1, j+r+vdisp  );
//               I_2(idx,1) = right_image(i+c+hdisp  , j+r+vdisp-1);
//               I_2(idx,2) = right_image(i+c+hdisp  , j+r+vdisp  );
//               I_2(idx,3) = right_image(i+c+hdisp  , j+r+vdisp+1);
//               I_2(idx,4) = right_image(i+c+hdisp+1, j+r+vdisp  );
//               I_2(idx,5) = 1;
              I_2(idx,0) = right_image(i+c+hdisp-1, j+r+vdisp-1);
              I_2(idx,1) = right_image(i+c+hdisp  , j+r+vdisp-1);
              I_2(idx,2) = right_image(i+c+hdisp+1, j+r+vdisp-1);
              I_2(idx,3) = right_image(i+c+hdisp-1, j+r+vdisp  );
              I_2(idx,4) = right_image(i+c+hdisp  , j+r+vdisp  );
              I_2(idx,5) = right_image(i+c+hdisp+1, j+r+vdisp  );
              I_2(idx,6) = right_image(i+c+hdisp-1, j+r+vdisp+1);
              I_2(idx,7) = right_image(i+c+hdisp  , j+r+vdisp+1);
              I_2(idx,8) = right_image(i+c+hdisp+1, j+r+vdisp+1);
              I_2(idx,9) = 1;
              ++idx;
            }
          }

          // Solve for the coefficients
          Matrix<double> L = inverse(transpose(I_2)*I_2)*transpose(I_2);
          Vector<double> A = L*I_1;
          double S1 = sum(subvector(A,0,9));  
          Vector2 dX_1 ( ( A(5)-A(3) ) / S1, ( A(7)-A(1) ) / S1 );
//           double S1 = sum(subvector(A,0,5));  
//           Vector2 dX_1 ( ( A(4)-A(0) ) / S1, ( A(3)-A(1) ) / S1 );
          Matrix2x2 CV_1;

          Vector<double> sum_vec = I_1;
          for (unsigned z = 0; z < sum_vec.size(); ++z)
            sum_vec[z] -= A[9];  
          for (unsigned z = 0; z < 9; ++z)
            sum_vec -= (A[z]*select_col(I_2,z));
          Vector<double> square_sum_vec = elem_prod(sum_vec, sum_vec) / (2*kern_width+1)*(2*kern_height+1);
          double sigma_1 = sqrt(sum(square_sum_vec));

          double scaled_sigma_1 = (sigma_1*sigma_1)/(S1*S1);;
          Vector<double> L1 = select_row(L,1);
          Vector<double> L3 = select_row(L,3);
          Vector<double> L5 = select_row(L,5);
          Vector<double> L7 = select_row(L,7);
          CV_1(0,0) = scaled_sigma_1 * dot_prod(L5-L3,L5-L3);   // var(dx)
          CV_1(0,1) = scaled_sigma_1 * dot_prod(L5-L3,L7-L1);   // cov(dx,dy)
          CV_1(1,0) = CV_1(0,1);
          CV_1(1,1) = scaled_sigma_1 * dot_prod(L7-L1,L7-L1);   // var(dy)
//           Vector<double> L0 = select_row(L,0);
//           Vector<double> L4 = select_row(L,4);
//           Vector<double> L1 = select_row(L,1);
//           Vector<double> L3 = select_row(L,3);
//           CV_1(0,0) = scaled_sigma_1 * dot_prod(L4-L0,L4-L0);   // var(dx)
//           CV_1(0,1) = scaled_sigma_1 * dot_prod(L4-L0,L3-L1);   // cov(dx,dy)
//           CV_1(1,0) = CV_1(0,1);
//           CV_1(1,1) = scaled_sigma_1 * dot_prod(L3-L1,L3-L1);   // var(dy)

          // Copy the image data into I_1 and I_2
          idx = 0;
          for (int r = -kern_half_height; r <= kern_half_height; ++r) {
            for (int c = -kern_half_width; c <= kern_half_width; ++c) { 
              I_1(idx) = right_image(i+hdisp+c, j+vdisp+r);
              I_2(idx,0) = left_image(i+c-1, j+r-1);
              I_2(idx,1) = left_image(i+c  , j+r-1);
              I_2(idx,2) = left_image(i+c+1, j+r-1);
              I_2(idx,3) = left_image(i+c-1, j+r  );
              I_2(idx,4) = left_image(i+c  , j+r  );
              I_2(idx,5) = left_image(i+c+1, j+r  );
              I_2(idx,6) = left_image(i+c-1, j+r+1);
              I_2(idx,7) = left_image(i+c  , j+r+1);
              I_2(idx,8) = left_image(i+c+1, j+r+1);
              I_2(idx,9) = 1;
//               I_2(idx,0) = left_image(i+c-1, j+r  );
//               I_2(idx,1) = left_image(i+c  , j+r-1);
//               I_2(idx,2) = left_image(i+c  , j+r  );
//               I_2(idx,3) = left_image(i+c  , j+r+1);
//               I_2(idx,4) = left_image(i+c+1, j+r  );
//               I_2(idx,5) = 1;
              ++idx;
            }
          }

          // Solve for the coefficients
          L = inverse(transpose(I_2)*I_2)*transpose(I_2);
          Vector<double> B = L * I_1;
          double S2 = sum(subvector(B,0,9));  
          Vector2 dX_2 ( ( B(3)-B(5) ) / S2, ( B(1)-B(7) ) / S2 );
//           double S2 = sum(subvector(B,0,5));  
//           Vector2 dX_2 ( ( B(0)-B(4) ) / S2, ( B(1)-B(3) ) / S2 );
          Matrix2x2 CV_2;

          sum_vec = I_1;
          for (unsigned z = 0; z < sum_vec.size(); ++z)
            sum_vec[z] -= B[9];  
          for (unsigned z = 0; z < 9; ++z)
            sum_vec -= (B[z]*select_col(I_2,z));
          square_sum_vec = elem_prod(sum_vec, sum_vec) / (2*kern_width+1)*(2*kern_height+1);
          double sigma_2 = sqrt(sum(square_sum_vec));

          double scaled_sigma_2 = (sigma_2*sigma_2)/(S2*S2);
//           L4 = select_row(L,4);
//           L0 = select_row(L,0);
//           L3 = select_row(L,3);
//           L1 = select_row(L,1);
//           CV_2(0,0) = scaled_sigma_2 * dot_prod(L0-L4,L0-L4);   // var(dx)
//           CV_2(0,1) = scaled_sigma_2 * dot_prod(L0-L4,L1-L3);   // cov(dx,dy)
//           CV_2(1,0) = CV_2(0,1);
//           CV_2(1,1) = scaled_sigma_2 * dot_prod(L1-L3,L1-L3);   // var(dy)
          L1 = select_row(L,1);
          L3 = select_row(L,3);
          L5 = select_row(L,5);
          L7 = select_row(L,7);
          CV_2(0,0) = scaled_sigma_2 * dot_prod(L3-L5,L3-L5);   // var(dx)
          CV_2(0,1) = scaled_sigma_2 * dot_prod(L3-L5,L1-L7);   // cov(dx,dy)
          CV_2(1,0) = CV_2(0,1);
          CV_2(1,1) = scaled_sigma_2 * dot_prod(L1-L7,L1-L7);   // var(dy)

          Matrix2x2 CV_1i = inverse(CV_1);
          Matrix2x2 CV_2i = inverse(CV_2);
          Vector2 update = inverse(CV_1i + CV_2i) * (CV_1i*dX_1 + CV_2i*dX_2);
  
          disparity_map(i,j).h() += update[0]; // dx
          disparity_map(i,j).v() += update[1]; // dy

//           if (i == 237 && j == 237) {
//             std::cout << "\n\t\t" << dX_1[0] << " " << dX_1[1] << "   " << A << " " << sum(subvector(A,0,9)) << "  " << sigma_1 << "  \n";
//             std::cout << "\t\t" << dX_2[0] << " " << dX_2[1] << "   " << B << " " << sum(subvector(B,0,9)) << "  " << sigma_2 << "  \n";
//             std::cout << "\t\t" << update[0] << " " << update[1] << "   " << CV_1 << "   " << CV_2 << "\n";
//             std::cout << "Pixel 237x237: " << disparity_map(i,j).h() << " " << disparity_map(i,j).v() << " \n\n";
//           }
        }
      }
    }

    if (verbose) vw_out(InfoMessage, "stereo") << "\tPerforming sub-pixel correlation... done.                 \n";
  }


  template <class ChannelT> 
  void subpixel_correlation(ImageView<PixelDisparity<float> > &disparity_map,
                            ImageView<ChannelT> const& left_image,
                            ImageView<ChannelT> const& right_image,
                            int kern_width, int kern_height,
                            bool do_horizontal_subpixel = true,
                            bool do_vertical_subpixel = true,
                            bool verbose = false) {

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
    // each row in A is [ x^2 y^2 xy x y 1] (our 2d hyperbolic surface)
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

//           if (r == 237 && c == 237)
//             std::cout << "Pixel 237x237: " << hdisp << " " << vdisp << "      ";
        
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

//               if (r == 237 && c == 237) {
//                 std::cout << offset(0) << " " << offset(1) << "      \n";
//                 std::cout << "Pixel 237x237: " << disparity_map(c,r).h() << " " << disparity_map(c,r).v() << "      \n\n";
//               }

            } else {
              disparity_map(c,r) = PixelDisparity<float>();
              //            std::cout << "Bad offset: " << offset(0) << " " << offset(1) << "\n";
            }
          }
        } 
      }
    }
    if (verbose) vw_out(InfoMessage, "stereo") << "\tPerforming sub-pixel correlation... done.                 \n";
  }

  /// This routine cross checks L2R and R2L, placing the final version
  /// of the disparity map in L2R.
  void cross_corr_consistency_check(ImageView<PixelDisparity<float> > &L2R, 
                                    ImageView<PixelDisparity<float> > &R2L,
                                    double cross_corr_threshold, bool verbose = false);

}} // namespace vw::stereo

#endif // __VW_STEREO_CORRELATOR_H__
