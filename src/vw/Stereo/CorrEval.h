#ifndef __VW_STEREO_CORR_EVAL_H__
#define __VW_STEREO_CORR_EVAL_H__

#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/PixelMask.h>
#include <vw/Core/Exception.h>
#include <cmath>

// Given two masked aligned images and a float disparity from left to
// right image, find at each pixel the normalized cross-correlation
// (NCC) of the patch of image1 around given pixel of given size, and
// corresponding patch of image2 around pixel + disparity(pixel). Use
// bilinear interpolation (unless rounding to int). Invalid pixels and
// pixels out of range are excluded. If the denominator is 0, set NCC
// to 0.

// Formula:
// 
// NCC(pix) = sum(image1(pix + s) * image2(pix + disp(pix) + s) /
//              sqrt ( sum(image1(pix + s)^2) * sum(image2(pix + disp(pix) + s)^2) ).
//
// Here, s varies over [-kx, kx] x [-ky, ky], where kernel_size =
// (2*kx + 1, 2*ky + 1).

// If 'metric' is not 'ncc' but 'stddev', calculate the average of
// stddev of the images in the two patches.

// See also CostFunctions.h for an implementation of NCC optimized to
// be applied over the entire image at once for a given fixed integer
// disp value.

namespace vw { namespace stereo {

// Calculate left and right patches. It assumes everything was setup properly
// and all the checks have been done.
void calc_patches(// Inputs
                  BBox2i const& bbox, Vector2i const& kernel_size,
                  bool round_to_int, PixelMask<Vector2f> const& disp,
                  BBox2i const& left_box, BBox2i const& right_box,
                  ImageView<PixelMask<float>> const& left,
                  ImageView<PixelMask<float>> const& right,
                  int col, int row, // patches are around this col and row
                  // Outputs
                  ImageView<PixelMask<float>> & left_patch,
                  ImageView<PixelMask<float>> & right_patch);
  
// Calc NCC. Return -1 on failure (normally NCC is non-negative).  
double calc_ncc(ImageView<PixelMask<float>> const& left_patch,
                ImageView<PixelMask<float>> const& right_patch);
  
// Calc stddev. Skip invalid pixels. Return -1 on failure (normally
// stddev is non-negative).
double calc_stddev(ImageView<PixelMask<float>> const& patch);
  
class CorrEval: public ImageViewBase<CorrEval> {
  ImageViewRef<PixelMask<float>> m_left, m_right;
  ImageViewRef<PixelMask<Vector2f>> m_disp;
  Vector2i m_kernel_size;
  std::string m_metric;
  int m_sample_rate;
  bool m_round_to_int;
  int m_extra_padding;

public:
  CorrEval(ImageViewRef<PixelMask<float>>    left,
           ImageViewRef<PixelMask<float>>    right,
           ImageViewRef<PixelMask<Vector2f>> disp,
           Vector2i                   const& kernel_size,
           std::string                const& metric,
           int sample_rate, bool round_to_int,
           int prefilter_mode = 0,
           float prefilter_kernel_width = 0.0):
    m_left(left), m_right(right), m_disp(disp),
    m_kernel_size(kernel_size), m_metric(metric),
    m_sample_rate(sample_rate), m_round_to_int(round_to_int),
    m_extra_padding((int)ceil(prefilter_kernel_width) + 5) {
    
    VW_ASSERT((m_left.cols() == m_disp.cols() && m_left.rows() == m_disp.rows()),
              vw::ArgumentErr()
              << "CorrEval: Left image and disparity must have the same dimensions.\n");

    VW_ASSERT(((m_kernel_size[0] > 0) && (m_kernel_size[0] % 2 == 1) &&
               (m_kernel_size[1] > 0) && (m_kernel_size[1] % 2 == 1)),
              vw::ArgumentErr()
              << "CorrEval: The kernel dimensions must be positive and odd.\n");

    VW_ASSERT((m_metric == "ncc" || metric == "stddev"),
              vw::ArgumentErr() << "CorrEval: Invalid metric: " << m_metric << ".\n");
  }
  
  typedef PixelMask<float> pixel_type;
  typedef PixelMask<float> result_type;
  typedef ProceduralPixelAccessor<CorrEval> pixel_accessor;

  inline int32 cols()   const { return m_left.cols(); }
  inline int32 rows()   const { return m_left.rows(); }
  inline int32 planes() const { return 1; }

  inline pixel_accessor origin() const { return pixel_accessor(*this, 0, 0); }

  inline pixel_type operator()( double/*i*/, double/*j*/, int32/*p*/ = 0 ) const {
    vw_throw(NoImplErr() << "CorrEval::operator() is not implemented.");
    return pixel_type();
  }

  typedef CropView<ImageView<pixel_type>> prerasterize_type;
  prerasterize_type prerasterize(BBox2i const& bbox) const;
  
  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
};
  
CorrEval corr_eval(ImageViewRef<PixelMask<float>>    left,
                   ImageViewRef<PixelMask<float>>    right,
                   ImageViewRef<PixelMask<Vector2f>> disp,
                   Vector2i                   const& kernel_size,
                   std::string                const& metric,
                   int sample_rate, bool round_to_int,
                   int prefilter_mode = 0,
                   float prefilter_kernel_width = 0.0) {
  return CorrEval(left, right, disp, kernel_size, metric,
                  sample_rate, round_to_int,
                  prefilter_mode, prefilter_kernel_width);
}
  
}} // end namespace vw::stereo

#endif //__VW_STEREO_CORR_EVAL__
