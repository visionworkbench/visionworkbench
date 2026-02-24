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

/// \file ImageTransform.h
///
/// Spatial transformation views: transpose, rotation by 90-degree
/// increments, vertical flip, and weighted RGB-to-gray conversion.
/// All views are shallow (no data copy). Split from Manipulation.h
/// for compile-time reduction.

#ifndef __VW_IMAGE_IMAGETRANSFORM_H__
#define __VW_IMAGE_IMAGETRANSFORM_H__

#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/PerPixelViews.h>

namespace vw {

// Transpose

template <class ChildT>
class TransposePixelAccessor {
  ChildT m_child;
public:
  typedef typename ChildT::pixel_type pixel_type;
  typedef typename ChildT::result_type result_type;
  typedef typename ChildT::offset_type offset_type;
  TransposePixelAccessor(ChildT const& acc) : m_child(acc) {}

  inline TransposePixelAccessor& next_col() { m_child.next_row(); return *this; }
  inline TransposePixelAccessor& prev_col() { m_child.prev_row(); return *this; }
  inline TransposePixelAccessor& next_row() { m_child.next_col(); return *this; }
  inline TransposePixelAccessor& prev_row() { m_child.prev_col(); return *this; }
  inline TransposePixelAccessor& next_plane() { m_child.next_plane(); return *this; }
  inline TransposePixelAccessor& prev_plane() { m_child.prev_plane(); return *this; }
  inline TransposePixelAccessor& advance(offset_type di, offset_type dj,
                                          ssize_t dp = 0) {
    m_child.advance(dj, di, dp); return *this;
  }

  inline result_type operator*() const { return *m_child; }
};

template <class ImageT>
class TransposeView: public ImageViewBase<TransposeView<ImageT>> {
  ImageT m_child;
public:
  typedef typename ImageT::pixel_type pixel_type;
  typedef typename ImageT::result_type result_type;
  typedef TransposePixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

  TransposeView(ImageT const& image) : m_child(image) {}

  inline int32 cols() const { return m_child.rows(); }
  inline int32 rows() const { return m_child.cols(); }
  inline int32 planes() const { return m_child.planes(); }

  inline pixel_accessor origin() const { return m_child.origin(); }

  inline result_type operator()(int32 i, int32 j, int32 p = 0) const {
    return m_child(j, i, p);
  }

  template <class ViewT>
  TransposeView const& operator=(ImageViewBase<ViewT> const& view) const {
    view.impl().rasterize(*this,
      BBox2i(0, 0, view.impl().cols(), view.impl().rows()));
    return *this;
  }

  template <class ViewT>
  TransposeView& operator=(ImageViewBase<ViewT> const& view) {
    *const_cast<const TransposeView*>(this) = view.impl();
    return *this;
  }

  ImageT const& child() const { return m_child; }

  typedef TransposeView<typename ImageT::prerasterize_type> prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& bbox) const {
    BBox2i child_bbox(bbox.min().y(), bbox.min().x(),
                      bbox.height(), bbox.width());
    return prerasterize_type(m_child.prerasterize(child_bbox));
  }
  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
};

template <class ImageT>
struct IsMultiplyAccessible<TransposeView<ImageT>>:
  public IsMultiplyAccessible<ImageT> {};

/// Transpose an image.
template <class ImageT>
TransposeView<ImageT> transpose(ImageViewBase<ImageT> const& v) {
  return TransposeView<ImageT>(v.impl());
}

// Rotate180

template <class ChildT>
class Rotate180PixelAccessor {
  ChildT m_child;
public:
  typedef typename ChildT::pixel_type  pixel_type;
  typedef typename ChildT::result_type result_type;
  typedef typename ChildT::offset_type offset_type;
  Rotate180PixelAccessor(ChildT const& image) : m_child(image) {}

  inline Rotate180PixelAccessor& next_col() { m_child.prev_col(); return *this; }
  inline Rotate180PixelAccessor& prev_col() { m_child.next_col(); return *this; }
  inline Rotate180PixelAccessor& next_row() { m_child.prev_row(); return *this; }
  inline Rotate180PixelAccessor& prev_row() { m_child.next_row(); return *this; }
  inline Rotate180PixelAccessor& next_plane() { m_child.next_plane(); return *this; }
  inline Rotate180PixelAccessor& prev_plane() { m_child.prev_plane(); return *this; }
  inline Rotate180PixelAccessor& advance(offset_type di, offset_type dj,
                                          ssize_t dp = 0) {
    m_child.advance(-di, -dj, dp); return *this;
  }

  inline result_type operator*() const { return *m_child; }
};

template <class ImageT>
class Rotate180View: public ImageViewBase<Rotate180View<ImageT>> {
  ImageT m_child;
public:
  typedef typename ImageT::pixel_type pixel_type;
  typedef typename ImageT::result_type result_type;
  typedef Rotate180PixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

  Rotate180View(ImageT const& image) : m_child(image) {}

  inline int32 cols() const { return m_child.cols(); }
  inline int32 rows() const { return m_child.rows(); }
  inline int32 planes() const { return m_child.planes(); }

  inline pixel_accessor origin() const {
    return m_child.origin().advance(cols() - 1, rows() - 1);
  }

  inline result_type operator()(int32 i, int32 j, int32 p = 0) const {
    return m_child(cols() - 1 - i, rows() - 1 - j, p);
  }

  template <class ViewT>
  Rotate180View const& operator=(ImageViewBase<ViewT> const& view) const {
    view.impl().rasterize(*this,
      BBox2i(0, 0, view.impl().cols(), view.impl().rows()));
    return *this;
  }

  template <class ViewT>
  Rotate180View& operator=(ImageViewBase<ViewT> const& view) {
    *const_cast<const Rotate180View*>(this) = view.impl();
    return *this;
  }

  ImageT const& child() const { return m_child; }

  typedef Rotate180View<typename ImageT::prerasterize_type> prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& bbox) const {
    BBox2i child_bbox(cols() - bbox.max().x(), rows() - bbox.max().y(),
                      bbox.width(), bbox.height());
    return prerasterize_type(m_child.prerasterize(child_bbox));
  }
  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
};

template <class ImageT>
struct IsMultiplyAccessible<Rotate180View<ImageT>>:
  public IsMultiplyAccessible<ImageT> {};

/// Rotate an image 180 degrees.
template <class ImageT>
Rotate180View<ImageT> rotate_180(ImageViewBase<ImageT> const& v) {
  return Rotate180View<ImageT>(v.impl());
}

// Rotate90CW

template <class ChildT>
class Rotate90CWPixelAccessor {
  ChildT m_child;
public:
  typedef typename ChildT::pixel_type pixel_type;
  typedef typename ChildT::result_type result_type;
  typedef typename ChildT::offset_type offset_type;

  Rotate90CWPixelAccessor(ChildT const& acc) : m_child(acc) {}
  inline Rotate90CWPixelAccessor& next_col() { m_child.prev_row(); return *this; }
  inline Rotate90CWPixelAccessor& prev_col() { m_child.next_row(); return *this; }
  inline Rotate90CWPixelAccessor& next_row() { m_child.next_col(); return *this; }
  inline Rotate90CWPixelAccessor& prev_row() { m_child.prev_col(); return *this; }
  inline Rotate90CWPixelAccessor& next_plane() { m_child.next_plane(); return *this; }
  inline Rotate90CWPixelAccessor& prev_plane() { m_child.prev_plane(); return *this; }
  inline Rotate90CWPixelAccessor& advance(offset_type di, offset_type dj,
                                           ssize_t dp = 0) {
    m_child.advance(dj, -di, dp); return *this;
  }

  inline result_type operator*() const { return *m_child; }
};

template <class ImageT>
class Rotate90CWView: public ImageViewBase<Rotate90CWView<ImageT>> {
  ImageT m_child;
public:
  typedef typename ImageT::pixel_type pixel_type;
  typedef typename ImageT::result_type result_type;
  typedef Rotate90CWPixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

  Rotate90CWView(ImageT const& image) : m_child(image) {}

  inline int32 cols() const { return m_child.rows(); }
  inline int32 rows() const { return m_child.cols(); }
  inline int32 planes() const { return m_child.planes(); }

  inline pixel_accessor origin() const {
    return m_child.origin().advance(0, cols() - 1);
  }

  inline result_type operator()(int32 i, int32 j, int32 p = 0) const {
    return m_child(j, cols() - 1 - i, p);
  }

  template <class ViewT>
  Rotate90CWView const& operator=(ImageViewBase<ViewT> const& view) const {
    view.impl().rasterize(*this,
      BBox2i(0, 0, view.impl().cols(), view.impl().rows()));
    return *this;
  }

  template <class ViewT>
  Rotate90CWView& operator=(ImageViewBase<ViewT> const& view) {
    *const_cast<const Rotate90CWView*>(this) = view.impl();
    return *this;
  }

  ImageT const& child() const { return m_child; }

  typedef Rotate90CWView<typename ImageT::prerasterize_type> prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& bbox) const {
    BBox2i child_bbox(bbox.min().y(), cols() - bbox.max().x(),
                      bbox.height(), bbox.width());
    return prerasterize_type(m_child.prerasterize(child_bbox));
  }
  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
};

template <class ImageT>
struct IsMultiplyAccessible<Rotate90CWView<ImageT>>:
  public IsMultiplyAccessible<ImageT> {};

/// Rotate an image 90 degrees clockwise.
template <class ImageT>
Rotate90CWView<ImageT> rotate_90_cw(ImageViewBase<ImageT> const& v) {
  return Rotate90CWView<ImageT>(v.impl());
}

// Rotate90CCW

template <class ChildT>
class Rotate90CCWPixelAccessor {
  ChildT m_child;
public:
  typedef typename ChildT::pixel_type pixel_type;
  typedef typename ChildT::result_type result_type;
  typedef typename ChildT::offset_type offset_type;
  Rotate90CCWPixelAccessor(ChildT const& acc) : m_child(acc) {}

  inline Rotate90CCWPixelAccessor& next_col() { m_child.next_row(); return *this; }
  inline Rotate90CCWPixelAccessor& prev_col() { m_child.prev_row(); return *this; }
  inline Rotate90CCWPixelAccessor& next_row() { m_child.prev_col(); return *this; }
  inline Rotate90CCWPixelAccessor& prev_row() { m_child.next_col(); return *this; }
  inline Rotate90CCWPixelAccessor& next_plane() { m_child.next_plane(); return *this; }
  inline Rotate90CCWPixelAccessor& prev_plane() { m_child.prev_plane(); return *this; }
  inline Rotate90CCWPixelAccessor& advance(offset_type di, offset_type dj,
                                            ssize_t dp = 0) {
    m_child.advance(-dj, di, dp); return *this;
  }

  inline result_type operator*() const { return *m_child; }
};

template <class ImageT>
class Rotate90CCWView: public ImageViewBase<Rotate90CCWView<ImageT>> {
  ImageT m_child;
public:
  typedef typename ImageT::pixel_type pixel_type;
  typedef typename ImageT::result_type result_type;
  typedef Rotate90CCWPixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

  Rotate90CCWView(ImageT const& image) : m_child(image) {}

  inline int32 cols() const { return m_child.rows(); }
  inline int32 rows() const { return m_child.cols(); }
  inline int32 planes() const { return m_child.planes(); }

  inline pixel_accessor origin() const {
    return m_child.origin().advance(rows() - 1, 0);
  }

  inline result_type operator()(int32 i, int32 j, int32 p = 0) const {
    return m_child(rows() - 1 - j, i, p);
  }

  template <class ViewT>
  Rotate90CCWView const& operator=(ImageViewBase<ViewT> const& view) const {
    view.impl().rasterize(*this,
      BBox2i(0, 0, view.impl().cols(), view.impl().rows()));
    return *this;
  }

  template <class ViewT>
  Rotate90CCWView& operator=(ImageViewBase<ViewT> const& view) {
    *const_cast<const Rotate90CCWView*>(this) = view.impl();
    return *this;
  }

  ImageT const& child() const { return m_child; }

  typedef Rotate90CCWView<typename ImageT::prerasterize_type> prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& bbox) const {
    BBox2i child_bbox(rows() - bbox.max().y(), bbox.min().x(),
                      bbox.height(), bbox.width());
    return prerasterize_type(m_child.prerasterize(child_bbox));
  }
  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
};

template <class ImageT>
struct IsMultiplyAccessible<Rotate90CCWView<ImageT>>:
  public IsMultiplyAccessible<ImageT> {};

/// Rotate an image 90 degrees counter-clockwise.
template <class ImageT>
Rotate90CCWView<ImageT> rotate_90_ccw(ImageViewBase<ImageT> const& v) {
  return Rotate90CCWView<ImageT>(v.impl());
}

// FlipVertical

template <class ChildT>
class FlipVerticalPixelAccessor {
  ChildT m_child;
public:
  typedef typename ChildT::pixel_type pixel_type;
  typedef typename ChildT::result_type result_type;
  typedef typename ChildT::offset_type offset_type;
  FlipVerticalPixelAccessor(ChildT const& acc) : m_child(acc) {}

  inline FlipVerticalPixelAccessor& next_col() { m_child.next_col(); return *this; }
  inline FlipVerticalPixelAccessor& prev_col() { m_child.prev_col(); return *this; }
  inline FlipVerticalPixelAccessor& next_row() { m_child.prev_row(); return *this; }
  inline FlipVerticalPixelAccessor& prev_row() { m_child.next_row(); return *this; }
  inline FlipVerticalPixelAccessor& next_plane() { m_child.next_plane(); return *this; }
  inline FlipVerticalPixelAccessor& prev_plane() { m_child.prev_plane(); return *this; }
  inline FlipVerticalPixelAccessor& advance(offset_type di, offset_type dj,
                                             ssize_t dp = 0) {
    m_child.advance(di, -dj, dp); return *this;
  }

  inline result_type operator*() const { return *m_child; }
};

template <class ImageT>
class FlipVerticalView: public ImageViewBase<FlipVerticalView<ImageT>> {
  ImageT m_child;
public:
  typedef typename ImageT::pixel_type pixel_type;
  typedef typename ImageT::result_type result_type;
  typedef FlipVerticalPixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

  FlipVerticalView(ImageT const& image) : m_child(image) {}

  inline int32 cols() const { return m_child.cols(); }
  inline int32 rows() const { return m_child.rows(); }
  inline int32 planes() const { return m_child.planes(); }

  inline pixel_accessor origin() const {
    return m_child.origin().advance(0, rows() - 1);
  }

  inline result_type operator()(int32 i, int32 j, int32 p = 0) const {
    return m_child(i, rows() - 1 - j, p);
  }

  ImageT const& child() const { return m_child; }

  template <class ViewT>
  FlipVerticalView const& operator=(ImageViewBase<ViewT> const& view) const {
    view.impl().rasterize(*this,
      BBox2i(0, 0, view.impl().cols(), view.impl().rows()));
    return *this;
  }

  template <class ViewT>
  FlipVerticalView& operator=(ImageViewBase<ViewT> const& view) {
    *const_cast<const FlipVerticalView*>(this) = view.impl();
    return *this;
  }

  typedef FlipVerticalView<typename ImageT::prerasterize_type> prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& bbox) const {
    BBox2i child_bbox(bbox.min().x(), rows() - bbox.max().y(),
                      bbox.width(), bbox.height());
    return prerasterize_type(m_child.prerasterize(child_bbox));
  }
  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
};

template <class ImageT>
struct IsMultiplyAccessible<FlipVerticalView<ImageT>>:
  public IsMultiplyAccessible<ImageT> {};

/// Flip an image vertically.
template <class ImageT>
FlipVerticalView<ImageT> flip_vertical(ImageViewBase<ImageT> const& v) {
  return FlipVerticalView<ImageT>(v.impl());
}

// WeightedRGBToGray

/// A weighted rgb-to-gray pixel conversion functor.
class WeightedRGBToGrayFunctor {
  double m_rw, m_gw, m_bw;
public:
  template <class ArgsT> struct result {};
  template <class FuncT, class ChannelT>
  struct result<FuncT(PixelRGB<ChannelT>)> {
    typedef PixelGray<ChannelT> type;
  };
  template <class FuncT, class ChannelT>
  struct result<FuncT(PixelRGBA<ChannelT>)> {
    typedef PixelGrayA<ChannelT> type;
  };
  WeightedRGBToGrayFunctor(double rw, double gw, double bw):
    m_rw(rw), m_gw(gw), m_bw(bw) {}
  template <class ChannelT>
  inline PixelGrayA<ChannelT> operator()(PixelRGBA<ChannelT> const& rgb) const {
    return weighted_rgb_to_gray(rgb, m_rw, m_gw, m_bw);
  }
  template <class ChannelT>
  inline PixelGray<ChannelT> operator()(PixelRGB<ChannelT> const& rgb) const {
    return weighted_rgb_to_gray(rgb, m_rw, m_gw, m_bw);
  }
};

/// Weighted conversion from PixelRGBA to PixelGrayA using user-specified weights.
template <class ImageT>
inline UnaryPerPixelView<ImageT, WeightedRGBToGrayFunctor>
weighted_rgb_to_gray(ImageViewBase<ImageT> const& image,
                     double rw, double gw, double bw) {
  return UnaryPerPixelView<ImageT, WeightedRGBToGrayFunctor>(
    image.impl(), WeightedRGBToGrayFunctor(rw, gw, bw));
}

/// Weighted conversion from PixelRGBA to PixelGrayA using the default weights.
template <class ImageT>
inline UnaryPerPixelView<ImageT, WeightedRGBToGrayFunctor>
weighted_rgb_to_gray(ImageViewBase<ImageT> const& image) {
  WeightedRGBToGrayFunctor func(VW_RGB_TO_GRAY_R_WEIGHT,
                                VW_RGB_TO_GRAY_G_WEIGHT,
                                VW_RGB_TO_GRAY_B_WEIGHT);
  return UnaryPerPixelView<ImageT, WeightedRGBToGrayFunctor>(
    image.impl(), func);
}

} // namespace vw

#endif // __VW_IMAGE_IMAGETRANSFORM_H__
