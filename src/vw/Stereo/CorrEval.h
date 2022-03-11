#ifndef __VW_STEREO_CORR_EVAL_H__
#define __VW_STEREO_CORR_EVAL_H__

#include <vw/Image/ImageView.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/Interpolation.h>
#include <vw/Core/Exception.h>
#include <vw/Stereo/DisparityMap.h>

// Given two masked images and a float disparity from left to right
// image, find at each pixel the normalized cross-correlation (NCC) of
// the patch of image1 around given pixel of given size, and
// corresponding patch of image2 around pixel + disparity(pixel). Use
// bilinear interpolation. Invalid pixels and pixels out of range are
// excluded. If the denominator is 0, set NCC to 0.

// Formula:
// 
// NCC(pix) = sum( image1(pix + s) * image2(pix + disp(pix) + s) /
//              sqrt ( sum(image1(pix + s)^2) * sum(image2(pix + disp(pix) + s)^2) ).
//
// Here, s varies over [-kx, kx] x [-ky, ky], where kernel_size =
// (2*kx + 1, 2*ky + 1).
  
// See also CostFunctions.h for an implementation of NCC optimized to
// be applied over the entire image at once for a given fixed integer
// disp value.

namespace vw { namespace stereo {

class CorrEval: public ImageViewBase<CorrEval> {
  ImageViewRef<PixelMask<float>> m_left, m_right;
  ImageViewRef<PixelMask<Vector2f>> m_disp;
  Vector2i m_kernel_size;
  
public:
  CorrEval(ImageViewRef<PixelMask<float>> left, ImageViewRef<PixelMask<float>> right,
      ImageViewRef<PixelMask<Vector2f>> disp, Vector2i const& kernel_size):
    m_left(left), m_right(right), m_disp(disp), m_kernel_size(kernel_size) {
    
    VW_ASSERT((m_left.cols() == m_disp.cols() && m_left.rows() == m_disp.rows()),
              vw::ArgumentErr()
              << "CorrEval: Left image and disparity must have the same dimensions.");

    VW_ASSERT(((m_kernel_size[0] > 0) && (m_kernel_size[0] % 2 == 1) &&
               (m_kernel_size[1] > 0) && (m_kernel_size[1] % 2 == 1)),
              vw::ArgumentErr()
              << "CorrEval: The kernel dimensions must be positive and odd.");
  }
  
  typedef PixelMask<float> pixel_type;
  typedef PixelMask<float> result_type;
  typedef ProceduralPixelAccessor<CorrEval> pixel_accessor;

  inline int32 cols()   const { return m_left.cols(); }
  inline int32 rows()   const { return m_left.rows(); }
  inline int32 planes() const { return 1; }

  inline pixel_accessor origin() const { return pixel_accessor(*this, 0, 0); }

  inline pixel_type operator()( double/*i*/, double/*j*/, int32/*p*/ = 0 ) const {
    vw_throw(NoImplErr() << "CorrEval::operator()(...) is not implemented.");
    return pixel_type();
  }

  typedef CropView<ImageView<pixel_type>> prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& bbox) const {

    std::cout << "TODO(oalexan1): Must check with various tile sizes to ensure "
              << "no boundary artifacts!" << std::endl;
    
    // Bring the disparity for the given processing region in memory.
    // It was checked before that it has the correct extent.
    ImageViewRef<PixelMask<Vector2f>> disp = crop(m_disp, bbox);

    Vector2i half_kernel = m_kernel_size/2;
    
    // Need to be able to look beyond the current tile in left image
    // to be able to compute the NCC.
    BBox2i left_box = bbox;
    left_box.expand(half_kernel);

    // For the right image it is more complicated. Need to also
    // consider the disparity and interpolation.
    
    // TODO(oalexan1): if finding NCC in a neighborhood, need to
    // adjust below to use disp values in the neighborhood.
    BBox2i right_box;
    for (int col = 0; col < disp.cols(); col++) {
      for (int row = 0; row < disp.rows(); row++) {
        if (!is_valid(disp(col, row))) 
          continue;
        
        Vector2 left_pix  = bbox.min() + Vector2(col, row);
        Vector2 right_pix = left_pix + disp(col, row).child();
        right_box.grow(Vector2(floor(right_pix.x()), floor(right_pix.y())));
        right_box.grow(Vector2(ceil(right_pix.x()), ceil(right_pix.y())));
      }
    }

    right_box.expand(half_kernel); // Take into account the kernel
    right_box.expand(BilinearInterpolation::pixel_buffer); // Due to interpolation
    right_box.expand(2); // because right_box is exclusive in the upper right, and +1 just in case

    // An invalid pixel value used for edge extension
    PixelMask<float> nodata_pix(0); nodata_pix.invalidate();
    ValueEdgeExtension<PixelMask<float>> nodata_ext(nodata_pix); 

    // Crop portions of the inputs and bring them in memory. Extend them if need be
    // with invalid data.
    ImageView<PixelMask<float>> left  = crop(edge_extend(m_left, nodata_ext), left_box);
    ImageView<PixelMask<float>> right = crop(edge_extend(m_right, nodata_ext), right_box);
    
    // Interpolate into the right image
    ImageViewRef<PixelMask<float>> interp_right 
      = interpolate(right, BilinearInterpolation(), nodata_ext);

    // Allocate room for the patches
    ImageView<PixelMask<float>> left_patch(m_kernel_size[0], m_kernel_size[1]);
    ImageView<PixelMask<float>> right_patch(m_kernel_size[0], m_kernel_size[1]);
    
    // Process the tile
    ImageView<result_type> tile(bbox.width(), bbox.height());
    for (int col = 0; col < tile.cols(); col++) {
      for (int row = 0; row < tile.rows(); row++) {

        // Start the tile as invalid
        tile(col, row) = PixelMask<float>(0.0);
        tile(col, row).invalidate();

        PixelMask<Vector2f> d = disp(col, row);
        if (!is_valid(d))
          continue;
          
        // Iterate over the patch
        for (int c = 0; c < m_kernel_size[0]; c++) {
          for (int r = 0; r < m_kernel_size[1]; r++) {

            // Left pixel and right pixels in the full images.
            // The left pix is int, but the right pix is not because disp is float.
            // Make it double for added precision.
            Vector2i left_pix(col + bbox.min().x() + c - half_kernel[0],
                              row + bbox.min().y() + r - half_kernel[1]);
            Vector2  right_pix = Vector2(left_pix) + Vector2(d.child());

            // Compensate for the fact that we will access cropped
            // versions (sometime the cropped version will cut and
            // sometimes will extend the original images).
            left_pix  -= left_box.min();
            right_pix -= right_box.min();

            if (!bounding_box(left).contains(left_pix) || !bounding_box(right).contains(right_pix))
              vw_throw(ArgumentErr() << "Out of bounds in the NCC calculation. "
                       << "This is not expected.");
            
            left_patch(c, r)  = left(left_pix[0], left_pix[1]);           // access int pix
            right_patch(c, r) = interp_right(right_pix[0], right_pix[1]); // interp float pix
          }
        }

        // Do NCC
        double num = 0.0, den1 = 0.0, den2 = 0.0;
        for (int c = 0; c < m_kernel_size[0]; c++) {
          for (int r = 0; r < m_kernel_size[1]; r++) {
            if (!is_valid(left_patch(c, r) || !is_valid(right_patch(c, r)))) 
              continue;

            double a = left_patch(c, r).child();
            double b = right_patch(c, r).child();
            num  += a*b;
            den1 += a*a;
            den2 += b*b;
          }
        }

        if (den1 > 0.0 && den2 > 0.0) {
          tile(col, row).validate();
          tile(col, row).child() = num / sqrt(den1 * den2);
        }
        
      }
    }
    
    return prerasterize_type(tile, -bbox.min().x(), -bbox.min().y(),
                             cols(), rows());
  }

  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
};
  
CorrEval corr_eval(ImageViewRef<PixelMask<float>> left, ImageViewRef<PixelMask<float>> right,
                   ImageViewRef<PixelMask<Vector2f>> disp, Vector2i const& kernel_size){
  return CorrEval(left, right, disp, kernel_size);
}
  
}} // end namespace vw::stereo

#endif //__VW_STEREO_CORR_EVAL__
