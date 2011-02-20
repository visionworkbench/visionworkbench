// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_STEREO_CORRELATE_H__
#define __VW_STEREO_CORRELATE_H__

#include <vw/Stereo/DisparityMap.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/ImageMath.h>
#include <vw/Math/LinearAlgebra.h>
#include <limits.h>

namespace vw {
namespace stereo {

  enum CorrelatorType { ABS_DIFF_CORRELATOR = 0,
                        SQR_DIFF_CORRELATOR = 1,
                        NORM_XCORR_CORRELATOR = 2 };

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

  // Return type for the correlator .. remember SOAD can only return positive values
  template <class T>
  inline typename boost::enable_if_c<std::numeric_limits<typename CorrelatorAccumulatorType<T>::type>::has_quiet_NaN,typename CorrelatorAccumulatorType<T>::type>::type
  CorrelatorFailureValue() {
    typedef typename CorrelatorAccumulatorType<T>::type accum_type;
    return std::numeric_limits<accum_type>::quiet_NaN();
  }

  template <class T>
  inline typename boost::disable_if_c<std::numeric_limits<typename CorrelatorAccumulatorType<T>::type>::has_quiet_NaN,typename CorrelatorAccumulatorType<T>::type>::type
  CorrelatorFailureValue() {
    typedef typename CorrelatorAccumulatorType<T>::type accum_type;
    return std::numeric_limits<accum_type>::max();
  }

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

  // Gaussian blur pre-processing
  //
  // Default gaussian blur standard deviation is 1.5 pixels.
  class BlurStereoPreprocessingFilter {
    float m_blur_sigma;

  public:
    typedef ImageView<float> result_type;

    BlurStereoPreprocessingFilter(float blur_sigma = 1.5) : m_blur_sigma(blur_sigma) {}

    template <class ViewT>
    result_type operator()(ImageViewBase<ViewT> const& view) const {
      return gaussian_filter(channel_cast<float>(view.impl()),m_blur_sigma);
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
  template <class PixelT>
  inline typename CorrelatorAccumulatorType<typename CompoundChannelType<PixelT>::type>::type
  compute_soad(PixelT *img0, PixelT *img1,
               int32 r, int32 c,                   // row and column in img0
               int32 hdisp, int32 vdisp,           // Current disparity offset from (c,r) for img1
               int32 kern_width, int32 kern_height,// Kernel dimensions
               int32 width, int32 height) {        // Image dimensions
    typedef typename CompoundChannelType<PixelT>::type channel_type;
    typedef typename CorrelatorAccumulatorType<channel_type>::type accum_type;

    r -= kern_height/2;
    c -= kern_width/2;
    if (r<0         || c<0       || r+kern_height>=height       || c+kern_width>=width ||
        r+vdisp < 0 || c+hdisp<0 || r+vdisp+kern_height>=height || c+hdisp+kern_width>=width) {
      return CorrelatorFailureValue<PixelT>();
    }

    PixelT *new_img0 = img0;
    PixelT *new_img1 = img1;

    new_img0 += c + r*width;
    new_img1 += (c+hdisp) + (r+vdisp)*width;

    accum_type ret = 0;
    AbsDiffCostFunc cost_fn;
    for (int32 rr= 0; rr< kern_height; rr++) {
     for (int32 cc= 0; cc< kern_width; cc++) {
        ret += cost_fn(new_img0[cc], new_img1[cc]);
      }
      new_img0 += width;
      new_img1 += width;
    }
    return ret;
  }

  /// Compute the sum of the absolute difference between a template
  /// region taken from img1 and the window centered at (c,r) in img0.
  template <class ViewT>
  inline typename CorrelatorAccumulatorType<typename CompoundChannelType<typename ViewT::pixel_type>::type>::type
  compute_soad(ImageViewBase<ViewT> const& img0,
               ImageViewBase<ViewT> const& img1,
               int32 r, int32 c,                   // row and column in img0
               int32 hdisp, int32 vdisp,           // Current disparity offset from (c,r) for img1
               int32 kern_width, int32 kern_height,
               BBox2i const& left_bbox,
               BBox2i const& right_bbox) {// Kernel dimensions
    typedef typename CompoundChannelType<typename ViewT::pixel_type>::type channel_type;
    typedef typename CorrelatorAccumulatorType<channel_type>::type accum_type;

    r -= kern_height/2;
    c -= kern_width/2;
    if ( r < left_bbox.min().y()  || c<left_bbox.min().x() ||
         r + kern_height >= left_bbox.max().y()            ||
         c + kern_width >= left_bbox.max().x() ) {
      return CorrelatorFailureValue<typename ViewT::pixel_type>();
    }

    if ( r + vdisp < right_bbox.min().y() ||
         c + hdisp<right_bbox.min().x()   ||
         r + vdisp+kern_height >= right_bbox.max().y() ||
         c + hdisp+kern_width >= right_bbox.max().x() ) {
      return CorrelatorFailureValue<typename ViewT::pixel_type>();
    }


    accum_type ret = 0;
    AbsDiffCostFunc cost_fn;
    for (int32 rr= 0; rr< kern_height; rr++) {
     for (int32 cc= 0; cc< kern_width; cc++) {
       ret += cost_fn(img0.impl()(c+cc,r+rr), img1.impl()(c+cc+hdisp,r+rr+vdisp));
      }
    }
    return ret;
  }

  /// For a given set of images, compute the optimal disparity (minimum
  /// SOAD) at position left_image(i,j) for the given correlation window
  /// settings.
  ///
  /// The left_image and right_image must have the same dimensions, but
  /// this is only checked here if debugging is enabled.
  template <class ChannelT>
  inline PixelMask<Vector2f>
  compute_disparity(ImageView<ChannelT> &left_image,
                    ImageView<ChannelT> &right_image,
                    int32 i, int32 j,
                    int32 kern_width, int32 kern_height,
                    int32 min_h_disp, int32 max_h_disp,
                    int32 min_v_disp, int32 max_v_disp) {

    typedef typename CorrelatorAccumulatorType<ChannelT>::type accum_type;
    accum_type min_soad = std::numeric_limits<accum_type>::max(); // Impossibly large
    PixelMask<Vector2f> best_disparity;
    invalidate(best_disparity);
    for (int32 ii = min_h_disp; ii <= max_h_disp; ++ii) {
      for (int32 jj = min_v_disp; jj <= max_v_disp; ++jj) {
        accum_type soad = compute_soad(&(left_image(0,0)), &(right_image(0,0)),
                                       j, i, ii, jj,kern_width, kern_height,
                                       left_image.cols(), left_image.rows());
        if (soad != CorrelatorFailureValue<ChannelT>()
            && soad < min_soad) {
          min_soad = soad;
          validate( best_disparity );
          best_disparity[0] = ii;
          best_disparity[1] = jj;
        }
      }
    }
    return best_disparity;
  }

  inline int
  adjust_weight_image(ImageView<float> &weight,
                      ImageView<PixelMask<Vector2f> > const& disparity_map_patch,
                      ImageView<float> const& weight_template) {

    int32 center_pix_x = weight_template.cols()/2;
    int32 center_pix_y = weight_template.rows()/2;
    PixelMask<Vector2f> center_pix =
      disparity_map_patch(center_pix_x, center_pix_y);

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
      vw_throw(LogicErr() << "subpixel_weight: Sum of weight image was zero.  This isn't supposed to happen!");
    else
      weight /= sum;
    return num_good_pix;
  }

 void cross_corr_consistency_check(ImageView<PixelMask<Vector2f> > &L2R,
                                  ImageView<PixelMask<Vector2f> > const& R2L,
                                  double cross_corr_threshold, bool verbose = false);

#include <vw/Stereo/Correlate.tcc>


}} // namespace vw::stereo

#endif // __VW_STEREO_CORRELATOR_H__
