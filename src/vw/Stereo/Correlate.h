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

  /// Given a type, these traits classes help to determine a suitable
  /// working type for accumulation operations in the correlator
  template <class T> struct CorrelatorAccumulatorType {};
  template <> struct CorrelatorAccumulatorType<vw::uint8>   { typedef vw::uint16  type; };
  template <> struct CorrelatorAccumulatorType<vw::int8>    { typedef vw::uint16  type; };
  template <> struct CorrelatorAccumulatorType<vw::uint16>  { typedef vw::uint32  type; };
  template <> struct CorrelatorAccumulatorType<vw::int16>   { typedef vw::uint32  type; };
  template <> struct CorrelatorAccumulatorType<vw::uint32>  { typedef vw::uint64  type; };
  template <> struct CorrelatorAccumulatorType<vw::int32>   { typedef vw::uint64  type; };
  template <> struct CorrelatorAccumulatorType<vw::float32> { typedef vw::float32 type; };
  template <> struct CorrelatorAccumulatorType<vw::float64> { typedef vw::float64 type; };

VW_DEFINE_EXCEPTION(CorrelatorErr, vw::Exception);

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

  class AbsDiffCostFunc {
    // These functors allow us to specialize the behavior of the image
    // differencing operation, which is part of measuring the sum of
    // absolute difference (SOAD) between two images.  For 8-bit slog
    // images, we take the xor (^) of the two images.
    static inline float absdiff (const float val1, const float val2) { return fabs(val1-val2); }
    static inline double absdiff (const double val1, const double val2) { return fabs(val1-val2); }
    static inline uint8 absdiff (const uint8 val1, const uint8 val2) { return val1 > val2 ? (val1-val2) : (val2-val1); }
    static inline uint16 absdiff (const uint16 val1, const uint16 val2) { return val1 > val2 ? (val1-val2) : (val2-val1); }
    static inline uint32 absdiff (const uint32 val1, const uint32 val2) { return val1 > val2 ? (val1-val2) : (val2-val1); }
    static inline int8 absdiff (const int8 val1, const int8 val2) { return abs(val1-val2); }
    static inline int16 absdiff (const int16 val1, const int16 val2) { return abs(val1-val2); }
    static inline int32 absdiff (const int32 val1, const int32 val2) { return abs(val1-val2); }
    
  public:
    template <class ChannelT>
    ChannelT operator()(ChannelT const& x, ChannelT const& y) {
      return absdiff(x, y);
    }
  };

  struct SqDiffCostFunc {
    template <class ChannelT>
    ChannelT operator()(ChannelT const& x, ChannelT const& y) {
      return (x - y) * (x - y);
    }
  };

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
  
    typename CorrelatorAccumulatorType<ChannelT>::type ret = 0;
    AbsDiffCostFunc cost_fn;
    for (int rr= 0; rr< kern_height; rr++) {
     for (int cc= 0; cc< kern_width; cc++) {
        ret += cost_fn(new_img0[cc], new_img1[cc]);
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
                            bool verbose = false);

  /// This routine cross checks L2R and R2L, placing the final version
  /// of the disparity map in L2R.
  void cross_corr_consistency_check(ImageView<PixelDisparity<float> > &L2R, 
                                    ImageView<PixelDisparity<float> > &R2L,
                                    double cross_corr_threshold, bool verbose = false);

}} // namespace vw::stereo

#endif // __VW_STEREO_CORRELATOR_H__
