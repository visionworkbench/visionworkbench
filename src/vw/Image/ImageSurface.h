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
    return BBox2i(Vector2i(0, 0), Vector2i(1, 1));
  }

  template <class PixelAccessorT>
  PixelMask<Vector3f> operator()(PixelAccessorT const& accessor_loc) const {
    PixelAccessorT acc = accessor_loc;

    if (is_transparent(*acc))
      return PixelMask<Vector3f>();
    float alt1 = *acc;

    acc.advance(1, 0);
    if (is_transparent(*acc))
      return PixelMask<Vector3f>();
    float alt2 = *acc;

    acc.advance(-1, 1);
    if (is_transparent(*acc))
      return PixelMask<Vector3f>();
    float alt3 = *acc;

    Vector3f n1(m_u_scale, 0, alt2 - alt1);
    Vector3f n2(0, m_v_scale, alt3 - alt1);

    return normalize(cross_prod(n1, n2));
  }
};

template <class ViewT>
UnaryPerPixelAccessorView<EdgeExtensionView<ViewT, ConstantEdgeExtension>,
                          ComputeNormalsFunc>
compute_normals(ImageViewBase<ViewT> const& image,
                float u_scale, float v_scale) {
  return UnaryPerPixelAccessorView<
    EdgeExtensionView<ViewT, ConstantEdgeExtension>, ComputeNormalsFunc>(
    edge_extend(image.impl(), ConstantEdgeExtension()),
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
