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

/// \file ImageSurface.h
///
/// Surface-analysis views: MeanFillTransparent, ComputeNormals,
/// DotProd, and TwoThresholdFill. Split from Algorithms.h
/// for compile-time reduction.

#ifndef __VW_IMAGE_IMAGESURFACE_H__
#define __VW_IMAGE_IMAGESURFACE_H__

#include <vw/Image/AlgorithmFunctions.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Image/PerPixelAccessorViews.h>

namespace vw {

// MeanFillTransparent

// This is a preprocess step that sets the value of transparent
// pixels to the mean of the nearby opaque pixels.

template <class ImageT>
class MeanFillTransparent:
  public ImageViewBase<MeanFillTransparent<ImageT>> {
  ImageT m_image;

  template <class SrcAccessT>
  typename SrcAccessT::pixel_type
  inline accumulate_mean(SrcAccessT const& src) const {
    typedef typename SrcAccessT::pixel_type result_type;
    typedef typename CompoundChannelType<result_type>::type channel_type;
    typedef typename PixelWithoutAlpha<result_type>::type non_a_type;
    typedef typename AccumulatorType<channel_type>::type acc_type;
    typedef typename PixelChannelCast<non_a_type, acc_type>::type
      non_a_acc_type;
    non_a_acc_type sum_value;
    acc_type weight = 0;

    SrcAccessT px = src;
    px.next_col();
    sum_value += non_a_acc_type(non_alpha_channels(*px)) *
                 acc_type(alpha_channel(*px));
    weight += acc_type(alpha_channel(*px));
    px.next_row();
    sum_value += non_a_acc_type(non_alpha_channels(*px)) *
                 acc_type(alpha_channel(*px));
    weight += acc_type(alpha_channel(*px));
    px.prev_col();
    sum_value += non_a_acc_type(non_alpha_channels(*px)) *
                 acc_type(alpha_channel(*px));
    weight += acc_type(alpha_channel(*px));
    px.prev_col();
    sum_value += non_a_acc_type(non_alpha_channels(*px)) *
                 acc_type(alpha_channel(*px));
    weight += acc_type(alpha_channel(*px));
    px.prev_row();
    sum_value += non_a_acc_type(non_alpha_channels(*px)) *
                 acc_type(alpha_channel(*px));
    weight += acc_type(alpha_channel(*px));
    px.prev_row();
    sum_value += non_a_acc_type(non_alpha_channels(*px)) *
                 acc_type(alpha_channel(*px));
    weight += acc_type(alpha_channel(*px));
    px.next_col();
    sum_value += non_a_acc_type(non_alpha_channels(*px)) *
                 acc_type(alpha_channel(*px));
    weight += acc_type(alpha_channel(*px));
    px.next_col();
    sum_value += non_a_acc_type(non_alpha_channels(*px)) *
                 acc_type(alpha_channel(*px));
    weight += acc_type(alpha_channel(*px));

    if (weight <= 0)
      return result_type();

    result_type result(sum_value / weight);
    alpha_channel(result) = ChannelRange<channel_type>::min();
    return result;
  }

public:
  typedef typename ImageT::pixel_type pixel_type;
  typedef pixel_type result_type;
  typedef ProceduralPixelAccessor<MeanFillTransparent> pixel_accessor;

  MeanFillTransparent(ImageT const& image) : m_image(image) {}

  inline int32 cols() const { return m_image.cols(); }
  inline int32 rows() const { return m_image.rows(); }
  inline int32 planes() const { return m_image.planes(); }
  inline pixel_accessor origin() const { return pixel_accessor(*this); }

  inline result_type helper(int32 x, int32 y, int32 p,
                            true_type) const {
    if (is_transparent(m_image(x, y, p))) {
      if (x > 1 && y > 1 && x + 1 < cols() && y + 1 < rows())
        return accumulate_mean(m_image.origin().advance(x, y, p));
      else
        return accumulate_mean(
          edge_extend(m_image, ConstantEdgeExtension()).origin()
            .advance(x, y, p));
    }
    return m_image(x, y, p);
  }

  inline result_type helper(int32 x, int32 y, int32 p,
                            false_type) const {
    return m_image(x, y, p);
  }

  inline result_type operator()(int32 x, int32 y, int32 p = 0) const {
    return helper(x, y, p,
      typename PixelHasAlpha<pixel_type>::type());
  }

  typedef MeanFillTransparent<CropView<ImageView<result_type>>>
    prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& bbox) const {
    BBox2i actual = bbox;
    actual.expand(1);
    ImageView<result_type> src =
      edge_extend(m_image, actual, ConstantEdgeExtension());
    return prerasterize_type(crop(src,
      -actual.min()[0], -actual.min()[1], cols(), rows()));
  }

  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
};

template <class SourceT>
MeanFillTransparent<SourceT>
inline mean_fill_transparent(ImageViewBase<SourceT> const& src) {
  return MeanFillTransparent<SourceT>(src.impl());
}

// ComputeNormals

class ComputeNormalsFunc:
  public ReturnFixedType<PixelMask<Vector3f>> {
  float m_u_scale, m_v_scale;

public:
  ComputeNormalsFunc(float u_scale, float v_scale):
    m_u_scale(u_scale), m_v_scale(v_scale) {}

  BBox2i work_area() const {
    return BBox2i(Vector2i(-1, -1), Vector2i(1, 1));
  }

  // Compute the surface normal with Horn's 3x3 weighted central difference (the
  // same stencil as gdaldem hillshade). It is symmetric, so it avoids the
  // half-pixel shift and the quantization speckle of a one-sided difference.
  //
  // Validity is judged on the RAW data, before any fill: off-image cells (from
  // the caller's ZeroEdgeExtension) and no-data cells are both invalid. A normal
  // is produced only when the window is not degenerate: the center ROW (west,
  // center, east) must have at least 2 of 3 valid, and the center COLUMN (north,
  // center, south) must have at least 2 of 3 valid. That guarantees a real
  // east-west and a real north-south difference, so we never fabricate a flat
  // normal for an under-determined pixel. Note this still shades a pixel whose
  // center is no-data, as long as both opposite cardinal neighbors exist (Horn's
  // gradient does not use the center value), and it shades a pixel with only one
  // of west/east (or north/south) missing. Otherwise the pixel is masked.
  //
  // Once the gate passes, missing cells are filled with a stand-in (the center
  // value, or the mean of the valid cardinal neighbors when the center itself is
  // missing) so the partial differences are well defined. A full-data window
  // reduces to plain Horn, so interior and image-border pixels of a hole-free
  // DEM are unchanged. This is stricter than gdaldem near voids on purpose; the
  // pc_align alignment path calls gdaldem directly and is unaffected.
  template <class PixelAccessorT>
  PixelMask<Vector3f> operator()(PixelAccessorT const& accessor_loc) const {

    // Gather the 3x3 window values and raw validity. z[row+1][col+1], with row
    // growing downward (south) and col growing to the right (east).
    float z[3][3];
    bool  ok[3][3];
    for (int dy = -1; dy <= 1; dy++) {
      for (int dx = -1; dx <= 1; dx++) {
        PixelAccessorT a = accessor_loc;
        a.advance(dx, dy);
        ok[dy + 1][dx + 1] = !is_transparent(*a);
        z [dy + 1][dx + 1] = ok[dy + 1][dx + 1] ? float(*a) : 0.0f;
      }
    }

    // Degeneracy gate on the raw center row and center column, before any fill.
    int row_valid = int(ok[1][0]) + int(ok[1][1]) + int(ok[1][2]); // W, C, E
    int col_valid = int(ok[0][1]) + int(ok[1][1]) + int(ok[2][1]); // N, C, S
    if (row_valid < 2 || col_valid < 2)
      return PixelMask<Vector3f>();

    // Fill value for missing cells: the center if valid, else the mean of the
    // valid cardinal neighbors (the gate guarantees at least two of them exist).
    float fill;
    if (ok[1][1]) {
      fill = z[1][1];
    } else {
      float s = 0.0f; int n = 0;
      if (ok[1][0]) { s += z[1][0]; n++; } // W
      if (ok[1][2]) { s += z[1][2]; n++; } // E
      if (ok[0][1]) { s += z[0][1]; n++; } // N
      if (ok[2][1]) { s += z[2][1]; n++; } // S
      fill = s / float(n);
    }
    for (int r = 0; r < 3; r++)
      for (int c = 0; c < 3; c++)
        if (!ok[r][c]) z[r][c] = fill;

    // Horn's weighted central differences (east - west, south - north). The
    // factor of 8 from the weights is folded into the run components below.
    float dzdx = (z[0][2] + 2*z[1][2] + z[2][2]) - (z[0][0] + 2*z[1][0] + z[2][0]);
    float dzdy = (z[2][0] + 2*z[2][1] + z[2][2]) - (z[0][0] + 2*z[0][1] + z[0][2]);

    Vector3f n1(8.0f * m_u_scale, 0, dzdx);
    Vector3f n2(0, 8.0f * m_v_scale, dzdy);

    return normalize(cross_prod(n1, n2));
  }
};

// Use ZeroEdgeExtension so that reads past the image edge come back as an
// invalid (masked) pixel. ComputeNormalsFunc then treats off-image cells as
// no-data when it applies its degeneracy gate, i.e. validity is judged on the
// raw data before any fill, exactly as gdaldem judges the true image extent.
template <class ViewT>
UnaryPerPixelAccessorView<EdgeExtensionView<ViewT, ZeroEdgeExtension>,
                          ComputeNormalsFunc>
compute_normals(ImageViewBase<ViewT> const& image,
                float u_scale, float v_scale) {
  return UnaryPerPixelAccessorView<
    EdgeExtensionView<ViewT, ZeroEdgeExtension>, ComputeNormalsFunc>(
    edge_extend(image.impl(), ZeroEdgeExtension()),
    ComputeNormalsFunc(u_scale, v_scale));
}

// DotProduct

/// Perform the dot product between each pixel and a constant vector.
class DotProdFunc:
  public ReturnFixedType<PixelMask<PixelGray<float>>> {
  Vector3f m_vec;
public:
  DotProdFunc(Vector3f const& vec) : m_vec(normalize(vec)) {}
  PixelMask<PixelGray<float>>
  operator()(PixelMask<Vector3f> const& pix) const {
    if (is_transparent(pix))
      return PixelMask<PixelGray<float>>();
    else
      return dot_prod(pix.child(), m_vec);
  }
};

template <class ViewT>
UnaryPerPixelView<ViewT, DotProdFunc>
dot_prod(ImageViewBase<ViewT> const& view, Vector3f const& vec) {
  return UnaryPerPixelView<ViewT, DotProdFunc>(view.impl(),
                                               DotProdFunc(vec));
}

// TwoThresholdFill

/// Apply a double threshold to an image.
template <class ImageT>
class TwoThresholdFill:
  public ImageViewBase<TwoThresholdFill<ImageT>> {

  ImageT const& m_image;
  int    m_expand_size;
  double m_low_threshold;
  double m_high_threshold;
  uint8  m_output_false;
  uint8  m_output_true;
public:
  TwoThresholdFill(ImageViewBase<ImageT> const& image, int expand_size,
                   double low_threshold, double high_threshold,
                   uint8 output_false = 0, uint8 output_true = 1):
    m_image(image.impl()), m_expand_size(expand_size),
    m_low_threshold(low_threshold), m_high_threshold(high_threshold),
    m_output_false(output_false), m_output_true(output_true) {}

  typedef uint8      pixel_type;
  typedef pixel_type result_type;
  typedef ProceduralPixelAccessor<TwoThresholdFill> pixel_accessor;

  inline int32 cols() const { return m_image.cols(); }
  inline int32 rows() const { return m_image.rows(); }
  inline int32 planes() const { return 1; }

  inline pixel_accessor origin() const {
    return pixel_accessor(*this, 0, 0);
  }

  inline pixel_type operator()(double i, double j, int32 p = 0) const {
    vw_throw(NoImplErr() << "operator()(...) is not implemented");
    return pixel_type();
  }

  typedef CropView<ImageView<pixel_type>> prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& bbox) const {
    BBox2i big_bbox = bbox;
    big_bbox.expand(m_expand_size);
    big_bbox.crop(bounding_box(m_image));

    ImageView<pixel_type> output_tile(big_bbox.width(),
                                      big_bbox.height());
    ValueEdgeExtension<pixel_type> edge_wrapper(pixel_type(0));

    ImageView<typename ImageT::pixel_type> input_tile =
      crop(m_image, big_bbox);

    // First pass: top-left to bottom-right
    for (int r = 0; r < output_tile.rows(); ++r) {
      for (int c = 0; c < output_tile.cols(); ++c) {
        if ((input_tile(c, r) > m_high_threshold) ||
            ((input_tile(c, r) > m_low_threshold) &&
             ((edge_wrapper(output_tile, c - 1, r - 1) > 0) ||
              (edge_wrapper(output_tile, c,     r - 1) > 0) ||
              (edge_wrapper(output_tile, c + 1, r - 1) > 0) ||
              (edge_wrapper(output_tile, c - 1, r    ) > 0))))
          output_tile(c, r) = m_output_true;
        else
          output_tile(c, r) = m_output_false;
      }
    }

    // Second pass: bottom-right to top-left
    for (int r = output_tile.rows() - 1; r >= 0; --r) {
      for (int c = output_tile.cols() - 1; c >= 0; --c) {
        if (output_tile(c, r) == m_output_true)
          continue;
        if ((input_tile(c, r) > m_low_threshold) &&
            ((edge_wrapper(output_tile, c + 1, r + 1) > 0) ||
             (edge_wrapper(output_tile, c,     r + 1) > 0) ||
             (edge_wrapper(output_tile, c - 1, r + 1) > 0) ||
             (edge_wrapper(output_tile, c + 1, r    ) > 0)))
          output_tile(c, r) = m_output_true;
      }
    }

    return prerasterize_type(output_tile,
      -big_bbox.min().x(), -big_bbox.min().y(), cols(), rows());
  }

  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
};

/// Applies a flood fill from pixels above the high threshold through
/// pixels above the low threshold.
template <class ImageT>
TwoThresholdFill<ImageT>
two_threshold_fill(ImageViewBase<ImageT> const& image, int expand_size,
                   double low_threshold, double high_threshold,
                   uint8 output_false = 0, uint8 output_true = 1) {
  return TwoThresholdFill<ImageT>(image.impl(), expand_size,
    low_threshold, high_threshold, output_false, output_true);
}

} // namespace vw

#endif // __VW_IMAGE_IMAGESURFACE_H__
