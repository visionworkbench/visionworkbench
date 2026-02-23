#include <vw/Image/ImageView.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/Interpolation.h>
#include <vw/Core/Exception.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Stereo/Correlation.h>
#include <vw/Stereo/CorrEval.h>

// See CorrEval.h for documentation.

namespace vw { namespace stereo {

// Calculate left and right patches. It assumes everything was setup properly
// and all the checks have been done.
void calc_patches(// Inputs
                  BBox2i const& bbox, Vector2i const& kernel_size, bool round_to_int,
                  PixelMask<Vector2f> const& disp,
                  BBox2i const& left_box, BBox2i const& right_box,
                  ImageView<PixelMask<float>> const& left,
                  ImageView<PixelMask<float>> const& right,
                  int col, int row, // patches are around this col and row
                  // Outputs
                  ImageView<PixelMask<float>> & left_patch,
                  ImageView<PixelMask<float>> & right_patch) {

  // An invalid pixel value used for edge extension
  PixelMask<float> nodata_pix(0); nodata_pix.invalidate();
  ValueEdgeExtension<PixelMask<float>> nodata_ext(nodata_pix); 
  
  // Interpolate into the right image. Avoid using an ImageViewRef to
  // avoid a per-pixel virtual function overhead. The 'auto' keyword
  // will use the exact type.
  auto interp_right = interpolate(right, BilinearInterpolation(), nodata_ext);
    
  Vector2i half_kernel = kernel_size/2;
  
  // Iterate over the patch
  for (int c = 0; c < kernel_size[0]; c++) {
    for (int r = 0; r < kernel_size[1]; r++) {
      
      // Left pixel and right pixels in the full images.
      // The left pix is int, but the right pix is not because disp is float.
      // Make it double for added precision.
      Vector2i left_pix(col + bbox.min().x() + c - half_kernel[0],
                        row + bbox.min().y() + r - half_kernel[1]);
      Vector2  right_pix = Vector2(left_pix) + Vector2(disp.child());
                                                       
      // Compensate for the fact that we will access cropped image
      // versions (which may either cut or extend the original
      // images).
      left_pix  -= left_box.min();
      right_pix -= right_box.min();

      // Sanity check for left_pix. We do not check right_pix, as maybe filtering
      // messed it up and it went out of bounds. In that case interpolation
      /// will simply return an invalid value.
      if (!bounding_box(left).contains(left_pix))
        vw_throw(ArgumentErr() << "Out of bounds in the NCC calculation. "
                 << "This is not expected.");
      
      left_patch(c, r)  = left(left_pix[0], left_pix[1]);           // access int pix
      
      if (!round_to_int) 
        right_patch(c, r) = interp_right(right_pix[0], right_pix[1]); // interp float pix
      else
        right_patch(c, r) = right(right_pix[0], right_pix[1]); // do not interpolate
    }
  }
}
  
// Calc NCC. Return -1 on failure (normally NCC is non-negative).  
double calc_ncc(ImageView<PixelMask<float>> const& left_patch,
                ImageView<PixelMask<float>> const& right_patch) {

  if (left_patch.cols() != right_patch.cols() || left_patch.rows() != right_patch.rows()) 
    vw_throw(ArgumentErr() << "The left and right patches have different dimensions.");
  
  double num = 0.0, den1 = 0.0, den2 = 0.0;
  for (int c = 0; c < left_patch.cols(); c++) {
    for (int r = 0; r < left_patch.rows(); r++) {
      if (!is_valid(left_patch(c, r) || !is_valid(right_patch(c, r)))) 
        continue;
      
      double a = left_patch(c, r).child();
      double b = right_patch(c, r).child();
      num  += a*b;
      den1 += a*a;
      den2 += b*b;
    }
  }
  
  if (den1 > 0.0 && den2 > 0.0) 
    return  num / sqrt(den1 * den2);
    
  return -1.0;
}

// Calc stddev. Skip invalid pixels. Return -1 on failure (normally
// stddev is non-negative).
double calc_stddev(ImageView<PixelMask<float>> const& patch) {

  // Find the mean
  int num = 0;
  double mean = 0.0;
  for (int c = 0; c < patch.cols(); c++) {
    for (int r = 0; r < patch.rows(); r++) {
      if (!is_valid(patch(c, r))) 
        continue;
      
      num  += 1;
      mean += patch(c, r).child();
    }
  }
  
  if (num == 0) 
    return -1.0;

  mean /= num;

  double sum = 0.0;
  num = 0.0;
  
  for (int c = 0; c < patch.cols(); c++) {
    for (int r = 0; r < patch.rows(); r++) {
      if (!is_valid(patch(c, r))) 
        continue;
      
      num += 1;
      sum += (patch(c, r).child() - mean) * (patch(c, r).child() - mean);
    }
  }
  
  if (num == 0) 
    return -1.0;

  return sqrt(sum / num);
}

CorrEval::prerasterize_type CorrEval::prerasterize(vw::BBox2i const& bbox) const {
  
  // Bring the disparity for the given processing region in memory.
  // It was checked before that it has the correct extent.
  ImageView<PixelMask<Vector2f>> disp = crop(m_disp, bbox);

  if (m_round_to_int) {
    for (int col = 0; col < disp.cols(); col++) {
      for (int row = 0; row < disp.rows(); row++) {
        // Round the disparity (both valid and invalid values, the validity is not affected)
        disp(col, row).child() = round(disp(col, row).child());
      }
    }
    
    // TODO(oalexan1): Consider subdividing regions as done in stereo
    // correlation.  For disparities that have integer values, and
    // which vary little over a large area, that may be more efficient
    // than the per-pixel approach. However, if the range of
    // disparities in a region is more than the kernel size, the
    // per-pixel approach should do better. So there has to be a check
    // and much testing.
    
    // Also, the best_of_search_convolution() logic in Correlation.cc
    // needs some modifications, since there the disparity with lowest
    // cost function is kept, but here we must keep the cost function
    // for the given known disparity regardless of cost function
    // value.
    
    //  ImageView<PixelMask<Vector2i>> int_disp(disp.cols(), disp.rows());
    //  std::vector<stereo::SearchParam> zones; 
    //  subdivide_regions(int_disp, bounding_box(int_disp),
    //   zones, m_kernel_size);
    // Now must iterate over regions.
  }

  Vector2i half_kernel = m_kernel_size/2;

  // Need to be able to look beyond the current tile in the left image
  // to compute the NCC. Add extra padding beyond the NCC kernel to
  // ensure operations like prefiltering have enough context at tile
  // boundaries, avoiding boundary artifacts.
  BBox2i left_box = bbox;
  left_box.expand(half_kernel + Vector2i(m_extra_padding, m_extra_padding));

  // For the right image it is more complicated. Need to also
  // consider the disparity and interpolation.
  // Note: The memory usage can be high for a large disparity.
  
  // TODO(oalexan1): When finding the curvature of NCC will need to further
  // expand the box given the neighborhood we will use then. 
  BBox2i right_box;
  for (int col = 0; col < disp.cols(); col++) {
    for (int row = 0; row < disp.rows(); row++) {
      if (!is_valid(disp(col, row))) 
        continue;

      if (col % m_sample_rate != 0 || row % m_sample_rate != 0) 
        continue;
      
      Vector2 left_pix  = bbox.min() + Vector2(col, row);
      Vector2 right_pix = left_pix + disp(col, row).child();
      right_box.grow(Vector2(floor(right_pix.x()), floor(right_pix.y())));
      right_box.grow(Vector2(ceil(right_pix.x()), ceil(right_pix.y())));
    }
  }
  
  right_box.expand(half_kernel); // Take into account the kernel
  right_box.expand(BilinearInterpolation::pixel_buffer); // Due to interpolation
  right_box.expand(2); // because right_box is exclusive in the upper-right, and +1 just in case
  right_box.expand(m_extra_padding); // Extra padding for prefilter context

  // An invalid pixel value used for edge extension
  PixelMask<float> nodata_pix(0); nodata_pix.invalidate();
  ValueEdgeExtension<PixelMask<float>> nodata_ext(nodata_pix); 
  
  // Crop portions of the inputs and bring them in memory. Extend them if need be
  // with invalid data to not go out of range later. Data validity will be checked.
  ImageView<PixelMask<float>> left  = crop(edge_extend(m_left, nodata_ext), left_box);
  ImageView<PixelMask<float>> right = crop(edge_extend(m_right, nodata_ext), right_box);
  
  // Allocate room for the patches
  ImageView<PixelMask<float>> left_patch(m_kernel_size[0], m_kernel_size[1]);
  ImageView<PixelMask<float>> right_patch(m_kernel_size[0], m_kernel_size[1]);

  // Create the tile with the result
  ImageView<result_type> tile(bbox.width(), bbox.height());
  for (int col = 0; col < tile.cols(); col++) {
    for (int row = 0; row < tile.rows(); row++) {
      
      // Start the tile as invalid
      tile(col, row) = PixelMask<float>(0.0);
      tile(col, row).invalidate();
      
      PixelMask<Vector2f> d = disp(col, row);
      if (!is_valid(d))
        continue;

      if (col % m_sample_rate != 0 || row % m_sample_rate != 0) 
        continue;
      
      calc_patches(// Inputs
                   bbox, m_kernel_size, m_round_to_int, d,
                   left_box, right_box,  
                   left, right,  
                   col, row,  // patches are around this col and row
                   // Outputs
                   left_patch, right_patch);
      
      if (m_metric == "ncc") {
        double ncc = calc_ncc(left_patch, right_patch);
        if (ncc >= 0) {
          tile(col, row).validate();
          tile(col, row).child() = ncc;
        }
      } else if (m_metric == "stddev") {
        double left_stddev = calc_stddev(left_patch);
        double right_stddev = calc_stddev(right_patch);
        if (left_stddev >= 0.0 && right_stddev >= 0.0) {
          tile(col, row).validate();
          tile(col, row).child() = (left_stddev + right_stddev)/2.0;
        }
      }
    }
  }

  return prerasterize_type(tile, -bbox.min().x(), -bbox.min().y(),
                           cols(), rows());
}
  
}} // end namespace vw::stereo

