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

/// \file ImageChannels.h
///
/// Channel and plane selection views, pixel/channel cast views.
/// Includes select_col/row/plane/channel, channels_to_planes,
/// pixel_cast, pixel_cast_rescale, channel_cast, and
/// channel_cast_round_and_clamp. Split from Manipulation.h
/// for compile-time reduction.

#ifndef __VW_IMAGE_IMAGECHANNELS_H__
#define __VW_IMAGE_IMAGECHANNELS_H__

#include <vw/Image/ImageView.h>
#include <vw/Image/PerPixelViews.h>

#include <boost/mpl/if.hpp>

namespace vw {

// select_col()

template <class ImageT>
class SelectColView: public ImageViewBase<SelectColView<ImageT>> {
  ImageT m_child;
  int32 m_col;
public:
  typedef typename ImageT::pixel_type pixel_type;
  typedef typename ImageT::result_type result_type;
  typedef typename ImageT::pixel_accessor pixel_accessor;

  SelectColView(ImageT const& image, int32 col):
    m_child(image), m_col(col) {}

  int32 cols() const { return 1; }
  inline int32 rows() const { return m_child.rows(); }
  inline int32 planes() const { return m_child.planes(); }

  inline pixel_accessor origin() const {
    return m_child.origin().advance(m_col, 0, 0);
  }
  inline result_type operator()(int32 /*i*/, int32 j, int32 p = 0) const {
    return m_child(m_col, j, p);
  }

  SelectColView const& operator=(SelectColView const& view) const {
    view.rasterize(*this, BBox2i(0, 0, view.cols(), view.rows()));
    return *this;
  }

  template <class ViewT>
  SelectColView const& operator=(ImageViewBase<ViewT> const& view) const {
    view.impl().rasterize(*this,
      BBox2i(0, 0, view.impl().cols(), view.impl().rows()));
    return *this;
  }

  typedef SelectColView<typename ImageT::prerasterize_type> prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& bbox) const {
    return prerasterize_type(m_child.prerasterize(bbox), m_col);
  }
  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
};

template <class ImageT>
struct IsMultiplyAccessible<SelectColView<ImageT>>:
  public IsMultiplyAccessible<ImageT> {};

/// Extracts a single column of an image.
template <class ImageT>
SelectColView<ImageT> select_col(ImageViewBase<ImageT> const& v, int32 col) {
  return SelectColView<ImageT>(v.impl(), col);
}

// select_row()

template <class ImageT>
class SelectRowView: public ImageViewBase<SelectRowView<ImageT>> {
  ImageT m_child;
  int32 m_row;
public:
  typedef typename ImageT::pixel_type pixel_type;
  typedef typename ImageT::result_type result_type;
  typedef typename ImageT::pixel_accessor pixel_accessor;

  SelectRowView(ImageT const& image, int32 row):
    m_child(image), m_row(row) {}

  inline int32 cols() const { return m_child.cols(); }
  inline int32 rows() const { return 1; }
  inline int32 planes() const { return m_child.planes(); }

  inline pixel_accessor origin() const {
    return m_child.origin().advance(0, m_row, 0);
  }

  inline result_type operator()(int32 i, int32 /*j*/, int32 p = 0) const {
    return m_child(i, m_row, p);
  }

  SelectRowView const& operator=(SelectRowView const& view) const {
    view.rasterize(*this, BBox2i(0, 0, view.cols(), view.rows()));
    return *this;
  }

  template <class ViewT>
  SelectRowView const& operator=(ImageViewBase<ViewT> const& view) const {
    view.impl().rasterize(*this,
      BBox2i(0, 0, view.impl().cols(), view.impl().rows()));
    return *this;
  }

  typedef SelectRowView<typename ImageT::prerasterize_type> prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& bbox) const {
    return prerasterize_type(m_child.prerasterize(bbox), m_row);
  }
  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
};

template <class ImageT>
struct IsMultiplyAccessible<SelectRowView<ImageT>>:
  public IsMultiplyAccessible<ImageT> {};

/// Extracts a single row of an image.
template <class ImageT>
SelectRowView<ImageT> select_row(ImageViewBase<ImageT> const& v, int32 row) {
  return SelectRowView<ImageT>(v.impl(), row);
}

// select_plane()

template <class ImageT>
class SelectPlaneView: public ImageViewBase<SelectPlaneView<ImageT>> {
  ImageT m_child;
  int32 m_plane;
public:
  typedef typename ImageT::pixel_type pixel_type;
  typedef typename ImageT::result_type result_type;
  typedef typename ImageT::pixel_accessor pixel_accessor;

  SelectPlaneView(ImageT const& image, int32 plane):
    m_child(image), m_plane(plane) {}

  inline int32 cols() const { return m_child.cols(); }
  inline int32 rows() const { return m_child.rows(); }
  inline int32 planes() const { return 1; }

  inline pixel_accessor origin() const {
    return m_child.origin().advance(0, 0, m_plane);
  }
  inline result_type operator()(int32 i, int32 j, int32 p = 0) const {
    return m_child(i, j, m_plane + p);
  }

  template <class ViewT>
  SelectPlaneView const& operator=(ImageViewBase<ViewT> const& view) const {
    view.impl().rasterize(*this,
      BBox2i(0, 0, view.impl().cols(), view.impl().rows()));
    return *this;
  }

  template <class ViewT>
  SelectPlaneView& operator=(ImageViewBase<ViewT> const& view) {
    *const_cast<const SelectPlaneView*>(this) = view.impl();
    return *this;
  }

  typedef SelectPlaneView<typename ImageT::prerasterize_type> prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& bbox) const {
    return prerasterize_type(m_child.prerasterize(bbox), m_plane);
  }
  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
};

template <class ImageT>
struct IsMultiplyAccessible<SelectPlaneView<ImageT>>:
  public IsMultiplyAccessible<ImageT> {};

/// Extracts a single plane of a multi-plane image.
template <class ImageT>
SelectPlaneView<ImageT> select_plane(ImageViewBase<ImageT> const& v,
                                     int32 plane) {
  return SelectPlaneView<ImageT>(v.impl(), plane);
}

// select_channel()

/// A channel selecting functor, used by select_channel().
template <class ImageT>
struct SelectChannelFunctor {
  int32 m_channel;
public:
  SelectChannelFunctor(int32 channel) : m_channel(channel) {}

  typedef typename CompoundChannelType<typename ImageT::pixel_type>::type
    channel_type;
  typedef typename CopyCVR<typename ImageT::result_type, channel_type>::type
    result_type;

  result_type operator()(typename ImageT::result_type pixel) const {
    return compound_select_channel<result_type>(pixel, m_channel);
  }
};

/// Extracts a single channel of a multi-channel image.
template <class ImageT>
UnaryPerPixelView<ImageT, SelectChannelFunctor<ImageT>>
inline select_channel(ImageViewBase<ImageT>& image, int32 channel) {
  return UnaryPerPixelView<ImageT, SelectChannelFunctor<ImageT>>(
    image.impl(), SelectChannelFunctor<ImageT>(channel));
}

/// Extracts a single channel of a multi-channel image (const overload).
template <class ImageT>
UnaryPerPixelView<ImageT, SelectChannelFunctor<const ImageT>>
inline select_channel(ImageViewBase<ImageT> const& image, int32 channel) {
  return UnaryPerPixelView<ImageT, SelectChannelFunctor<const ImageT>>(
    image.impl(), SelectChannelFunctor<const ImageT>(channel));
}

/// A convenience function to select the alpha channel of an image.
template <class ImageT>
UnaryPerPixelView<ImageT, SelectChannelFunctor<ImageT>>
inline select_alpha_channel(ImageViewBase<ImageT>& image) {
  if (!PixelHasAlpha<typename ImageT::pixel_type>::value)
    vw_throw(ArgumentErr()
             << "Image has no alpha channel in call to select_alpha_channel()");
  return select_channel(image,
    PixelNumChannels<typename ImageT::pixel_type>::value - 1);
}

/// A convenience function to select the alpha channel of an image
/// (const overload).
template <class ImageT>
UnaryPerPixelView<ImageT, SelectChannelFunctor<const ImageT>>
inline select_alpha_channel(ImageViewBase<ImageT> const& image) {
  if (!PixelHasAlpha<typename ImageT::pixel_type>::value)
    vw_throw(ArgumentErr()
             << "Image has no alpha channel in call to select_alpha_channel()");
  return select_channel(image,
    PixelNumChannels<typename ImageT::pixel_type>::value - 1);
}

// channels_to_planes()

/// A channels-to-planes pixel accessor adaptor.
template <class ChildT>
class ChannelsToPlanesAccessor {
  ChildT m_child;
  int32 m_channel;
public:
  typedef typename CompoundChannelType<typename ChildT::pixel_type>::type
    pixel_type;
  typedef typename CopyCVR<typename ChildT::result_type, pixel_type>::type
    result_type;
  typedef typename ChildT::offset_type offset_type;

  ChannelsToPlanesAccessor(ChildT const& acc) : m_child(acc), m_channel(0) {}
  inline ChannelsToPlanesAccessor& next_col() { m_child.next_col(); return *this; }
  inline ChannelsToPlanesAccessor& prev_col() { m_child.prev_col(); return *this; }
  inline ChannelsToPlanesAccessor& next_row() { m_child.next_row(); return *this; }
  inline ChannelsToPlanesAccessor& prev_row() { m_child.prev_row(); return *this; }
  inline ChannelsToPlanesAccessor& next_plane() { ++m_channel; return *this; }
  inline ChannelsToPlanesAccessor& prev_plane() { --m_channel; return *this; }
  inline ChannelsToPlanesAccessor& advance(offset_type di, offset_type dj,
                                            ssize_t dp = 0) {
    m_child.advance(di, dj); m_channel += dp; return *this;
  }

  inline result_type operator*() const {
    return compound_select_channel<result_type>(*m_child, m_channel);
  }
};

/// A view that turns a one plane, multi-channel view into a
/// multi-plane, one channel view.
template <class ImageT>
class ChannelsToPlanesView:
  public ImageViewBase<ChannelsToPlanesView<ImageT>> {
  ImageT m_child;
public:
  typedef typename CompoundChannelType<typename ImageT::pixel_type>::type
    pixel_type;
  typedef typename CopyCVR<typename ImageT::result_type, pixel_type>::type
    result_type;

  typedef typename boost::mpl::if_<
    IsCompound<typename ImageT::pixel_type>,
    ChannelsToPlanesAccessor<typename ImageT::pixel_accessor>,
    typename ImageT::pixel_accessor>::type pixel_accessor;

  ChannelsToPlanesView(ImageT const& image) : m_child(image) {
    VW_ASSERT(m_child.planes() == 1,
              ArgumentErr()
              << "ChannelsToPlanesView: The image must be single plane.");
  }

  inline int32 cols() const { return m_child.cols(); }
  inline int32 rows() const { return m_child.rows(); }
  inline int32 planes() const { return m_child.channels(); }

  inline pixel_accessor origin() const {
    return pixel_accessor(m_child.origin());
  }

  inline result_type operator()(int32 i, int32 j, int32 p = 0) const {
    return compound_select_channel<result_type>(m_child(i, j), p);
  }

  template <class ViewT>
  ChannelsToPlanesView const& operator=(ImageViewBase<ViewT> const& view) const {
    view.impl().rasterize(*this,
      BBox2i(0, 0, view.impl().cols(), view.impl().rows()));
    return *this;
  }

  template <class ViewT>
  ChannelsToPlanesView& operator=(ImageViewBase<ViewT> const& view) {
    *const_cast<const ChannelsToPlanesView*>(this) = view.impl();
    return *this;
  }

  ImageT const& child() const { return m_child; }

  typedef ChannelsToPlanesView<typename ImageT::prerasterize_type>
    prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& bbox) const {
    return prerasterize_type(m_child.prerasterize(bbox));
  }
  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
};

template <class ImageT>
struct IsMultiplyAccessible<ChannelsToPlanesView<ImageT>>:
  public IsMultiplyAccessible<ImageT> {};

/// Adapts a multi-channel image view so the channels are treated as planes.
template <class ImageT>
ChannelsToPlanesView<ImageT>
channels_to_planes(ImageViewBase<ImageT> const& v) {
  return ChannelsToPlanesView<ImageT>(v.impl());
}

// pixel_cast()

/// A pixel casting functor, used by pixel_cast().
template <class PixelT>
struct PixelCastFunctor: ReturnFixedType<PixelT> {
  template <class ArgT>
  inline PixelT operator()(ArgT pixel) const {
    return pixel_cast<PixelT>(pixel);
  }
};

/// Create a new image view by statically casting the pixels to a new type.
template <class PixelT, class ImageT>
inline UnaryPerPixelView<ImageT, PixelCastFunctor<PixelT>>
pixel_cast(ImageViewBase<ImageT> const& image) {
  return UnaryPerPixelView<ImageT, PixelCastFunctor<PixelT>>(image.impl());
}

// pixel_cast_rescale()

/// A pixel casting functor, used by pixel_cast_rescale().
template <class PixelT>
struct PixelCastRescaleFunctor: ReturnFixedType<PixelT> {
  template <class ArgT>
  inline PixelT operator()(ArgT pixel) const {
    return pixel_cast_rescale<PixelT>(pixel);
  }
};

/// Create a new image view by casting and rescaling the pixels to a new type.
template <class PixelT, class ImageT>
inline UnaryPerPixelView<ImageT, PixelCastRescaleFunctor<PixelT>>
pixel_cast_rescale(ImageViewBase<ImageT> const& image) {
  return UnaryPerPixelView<ImageT, PixelCastRescaleFunctor<PixelT>>(
    image.impl());
}

// channel_cast()

/// A pixel channel casting functor, used by channel_cast().
template <class ChannelT>
struct PixelChannelCastFunctor:
  UnaryReturnBinaryTemplateBind2nd<CompoundChannelCast, ChannelT> {
  template <class ArgT>
  inline typename CompoundChannelCast<ArgT, ChannelT>::type
  operator()(ArgT const& pixel) const {
    return channel_cast<ChannelT>(pixel);
  }
};

/// Create a new image view by statically casting the channels of the
/// pixels to a new type.
template <class ChannelT, class ImageT>
inline UnaryPerPixelView<ImageT, PixelChannelCastFunctor<ChannelT>>
channel_cast(ImageViewBase<ImageT> const& image) {
  return UnaryPerPixelView<ImageT, PixelChannelCastFunctor<ChannelT>>(
    image.impl());
}

// channel_cast_round_and_clamp()

/// A pixel channel casting, rounding and clamping functor.
/// Uses the function with the same name from vw/Image/PixelTypeInfo.h.
template <class ChannelT>
struct PixelChannelCastRoundClampFunctor:
  UnaryReturnBinaryTemplateBind2nd<CompoundChannelCast, ChannelT> {
  template <class ArgT>
  inline typename CompoundChannelCast<ArgT, ChannelT>::type
  operator()(ArgT const& pixel) const {
    return channel_cast_round_and_clamp<ChannelT>(pixel);
  }
};

/// Create a new image view by casting, rounding, and clamping the
/// channels of the pixels to a new type.
template <class ChannelT, class ImageT>
inline UnaryPerPixelView<ImageT, PixelChannelCastRoundClampFunctor<ChannelT>>
channel_cast_round_and_clamp(ImageViewBase<ImageT> const& image) {
  return UnaryPerPixelView<ImageT,
    PixelChannelCastRoundClampFunctor<ChannelT>>(image.impl());
}

// channel_cast_rescale()

/// A pixel channel casting and rescaling functor.
template <class ChannelT>
struct PixelChannelCastRescaleFunctor:
  UnaryReturnBinaryTemplateBind2nd<CompoundChannelCast, ChannelT> {
  template <class ArgT>
  inline typename CompoundChannelCast<ArgT, ChannelT>::type
  operator()(ArgT const& pixel) const {
    return channel_cast_rescale<ChannelT>(pixel);
  }
};

/// Create a new image view by casting and rescaling the channels of
/// the pixels to a new type.
template <class ChannelT, class ImageT>
inline UnaryPerPixelView<ImageT, PixelChannelCastRescaleFunctor<ChannelT>>
channel_cast_rescale(ImageViewBase<ImageT> const& image) {
  return UnaryPerPixelView<ImageT,
    PixelChannelCastRescaleFunctor<ChannelT>>(image.impl());
}

} // namespace vw

#endif // __VW_IMAGE_IMAGECHANNELS_H__
