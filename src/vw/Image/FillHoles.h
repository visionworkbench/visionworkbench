// __BEGIN_LICENSE__
//  Copyright (c) 2006-2026, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__

/// \file FillHoles.h
///
/// Naive hole-filling view for masked images. For better results,
/// use the grassfire-based fill_holes_grass() in InpaintView.h.

#ifndef __VW_IMAGE_FILL_HOLES_H__
#define __VW_IMAGE_FILL_HOLES_H__

#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/PixelAccessors.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/AlgorithmFunctions.h>

#include <algorithm>
#include <cmath>

namespace vw {

// A class to fill holes in images. Given the number
// hole_fill_len, which gives the size of holes, look left, right,
// up, and down, as far as hole_fill_len/2. Must have valid pixels
// both left and right, otherwise both up and down, else do
// nothing. The motivation here is to fill in only pixels
// "surrounded" by valid pixels.

// Use one of the two approaches:
//
// 1. Interpolate using the left and right values. Interpolate using
// the up and down values. Average the results. Fast.
//
// 2. Find the weighted average of the points in the window
// of size hole_fill_len centered at the current point. Slow.

// After holes are filled, do several passes to average the results.

// The image being passed in should be a PixelMask.
template <class ImageT>
class FillHoles: public ImageViewBase<FillHoles<ImageT>> {
  ImageT m_img;
  int m_hole_fill_mode, m_hole_fill_num_smooth_iter;
  int m_hole_fill_half, m_smooth_half;
  ImageView<double> m_dist_kernel, m_gauss_kernel;
  typedef typename ImageT::pixel_type PixelT;

public:

  typedef PixelT pixel_type;
  typedef PixelT result_type;
  typedef ProceduralPixelAccessor<FillHoles> pixel_accessor;

  FillHoles(ImageViewBase<ImageT> const& img, int hole_fill_mode,
            int hole_fill_num_smooth_iter,
            int hole_fill_len):
    m_img(img.impl()), m_hole_fill_mode(hole_fill_mode),
    m_hole_fill_num_smooth_iter(hole_fill_num_smooth_iter),
    m_hole_fill_half((hole_fill_len+1)/2), m_smooth_half(2) {

    if (m_hole_fill_mode == 2) {
      // Use the inverse square distance kernel
      int h = m_hole_fill_half;
      m_dist_kernel.set_size(2*h+1, 2*h+1);
      for (int c = 0; c < m_dist_kernel.cols(); c++) {
        for (int r = 0; r < m_dist_kernel.rows(); r++) {
          double r2 = double(c-h)*(c-h) + double(r-h)*(r-h);
          m_dist_kernel(c, r) = 1.0/std::max(r2, 1.0);
        }
      }
    }

    if (m_hole_fill_num_smooth_iter > 0) {
      // Use a gaussian kernel
      int h = m_smooth_half;
      m_gauss_kernel.set_size(2*h+1, 2*h+1);
      double val = 0.25; // value to reach at kernel edge
      double sigma = -log(val)/double(h*h);
      for (int c = 0; c < m_gauss_kernel.cols(); c++) {
        for (int r = 0; r < m_gauss_kernel.rows(); r++) {
          double r2 = double(c-h)*(c-h) + double(r-h)*(r-h);
          m_gauss_kernel(c, r) = exp(-sigma*r2);
        }
      }
    }
  }

  inline int32 cols  () const { return m_img.cols(); }
  inline int32 rows  () const { return m_img.rows(); }
  inline int32 planes() const { return 1; }

  inline pixel_accessor origin() const { return pixel_accessor(*this); }

  inline result_type operator()(size_t i, size_t j, size_t p=0) const {
    vw_throw(NoImplErr() << "FillHoles: operator() not implemented.\n");
  }

  typedef CropView<ImageView<PixelT>> prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& bbox) const {

    // Crop into an expanded box as to have enough pixels to do
    // averaging with given window at every pixel in the current box.
    int h = m_hole_fill_half; // shorten
    BBox2i biased_box = bbox;
    biased_box.expand(h+1);
    biased_box.crop(bounding_box(m_img));
    ImageView<PixelT> img(crop(m_img, biased_box));
    ImageView<PixelT> filled_img = copy(img);
    int nc = img.cols(), nr = img.rows(); // shorten

    for (int row = 0; row < nr; row++) {
      for (int col = 0; col < nc; col++) {

        if (is_valid(img(col, row))) continue; // skip valid

        // Look left, right, up, down, and find the closest valid pixels
        double r0 = -1, r1 = -1, c0 = -1, c1 = -1; // no good indices yet
        for (int k = row-1; k >= std::max(row-h, 0); k--)
          if (is_valid(img(col, k))) { r0 = k; break; }
        if (r0 >= 0) {
          // Found a point to the left, try to also find one to the right
          for (int k = row+1; k <= std::min(row + h, nr-1); k++)
            if (is_valid(img(col, k))) { r1 = k; break; }
        }
        for (int k = col-1; k >= std::max(col-h, 0); k--)
          if (is_valid(img(k, row))) { c0 = k; break; }
        if (c0 >= 0) {
          // Found a point up, try to also find one down
          for (int k = col+1; k <= std::min(col + h, nc-1); k++)
            if (is_valid(img(k, row))) { c1 = k; break; }
        }

        // skip if no good neighbors
        if ((r0 < 0 || r1 < 0) && (c0 < 0 || c1 < 0)) continue;

        double sum = 0.0;
        PixelT V; V.validate();
        if (m_hole_fill_mode == 1) {
          // Interpolate between left and right, then between top and
          // bottom. Average the results.
          if (r0 >= 0 && r1 >= 0) {
            V += ((r1-row)*img(col, r0) + (row-r0)*img(col, r1))/double(r1-r0);
            sum++;
          }
          if (c0 >= 0 && c1 >= 0) {
            V += ((c1-col)*img(c0, row) + (col-c0)*img(c1, row))/double(c1-c0);
            sum++;
          }
        } else {
          // Weighted average with given kernel
          for (int c = std::max(col-h, 0); c <= std::min(col+h, nc-1); c++) {
            for (int r = std::max(row-h, 0); r <= std::min(row+h, nr-1); r++) {
              if (!is_valid(img(c, r))) continue;
              double wt = m_dist_kernel(c-col+h, r-row+h);
              V   += wt*img(c, r);
              sum += wt;
            }
          }
        }

        if (sum > 0) filled_img(col, row) = V/sum;
      }
    }

    // Smooth the resulting image by repeated convolutions
    // with a small gaussian kernel
    for (int i = 0; i < m_hole_fill_num_smooth_iter; i++) {

      ImageView<PixelT> curr_img = copy(filled_img);

      for (int row = 0; row < nr; row++) {
        for (int col = 0; col < nc; col++) {
          if (is_valid(img(col, row))) continue; // skip valid
          if (!is_valid(filled_img(col, row))) continue; // don't add more valid

          int nh = m_smooth_half;
          double sum = 0.0;
          PixelT V; V.validate();
          for (int c = std::max(col-nh, 0); c <= std::min(col+nh, nc-1); c++) {
            for (int r = std::max(row-nh, 0); r <= std::min(row+nh, nr-1); r++) {
              if (!is_valid(curr_img(c, r))) continue;
              double wt = m_gauss_kernel(c-col+nh, r-row+nh);
              V   += wt*curr_img(c, r);
              sum += wt;
            }
          }
          if (sum > 0) filled_img(col, row) = V/sum;
        }
      }
    } // end smoothing iterations

    return prerasterize_type(filled_img,
                             -biased_box.min().x(), -biased_box.min().y(),
                             cols(), rows());
  }

  template <class ImgT>
  inline void rasterize(ImgT const& img, BBox2i const& bbox) const {
    vw::rasterize(prerasterize(bbox), img, bbox);
  }
};

// Fill holes. This algorithm is very naive and works not so well.
// The grassfire-based in-painting algorithm from ASP works much better.
template <class ImageT> FillHoles<ImageT>
fill_holes(ImageViewBase<ImageT> const& img, int hole_fill_mode,
           int hole_fill_num_smooth_iter, int hole_fill_len) {
  return FillHoles<ImageT>(img.impl(), hole_fill_mode,
                           hole_fill_num_smooth_iter, hole_fill_len);
}

} // namespace vw

#endif // __VW_IMAGE_FILL_HOLES_H__
