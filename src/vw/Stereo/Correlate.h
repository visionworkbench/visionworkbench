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
  void subpixel_correlation_affine_2d(ImageView<PixelDisparity<float> > &disparity_map,
                                      ImageView<ChannelT> const& left_image,
                                      ImageView<ChannelT> const& right_image,
                                      int kern_width, int kern_height,
                                      bool do_horizontal_subpixel = true,
                                      bool do_vertical_subpixel = true,
                                      bool verbose = false);

  template <class ChannelT> 
  void subpixel_correlation_parabola(ImageView<PixelDisparity<float> > &disparity_map,
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
