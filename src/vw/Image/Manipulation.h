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

/// \file Manipulation.h
///
/// Core image manipulation: copy, crop, and subsample views.
/// Spatial transforms (transpose, rotate, flip) are in ImageTransform.h.
/// Channel/plane selection and casting are in ImageChannels.h.

#ifndef __VW_IMAGE_MANIPULATION_H__
#define __VW_IMAGE_MANIPULATION_H__

#include <boost/mpl/logical.hpp>

#include <vw/Image/ImageView.h>
#include <vw/Image/AlgorithmFunctions.h>

// Backward compatibility: include the split-out headers so existing
// code that includes Manipulation.h keeps working.
#include <vw/Image/ImageTransform.h>
#include <vw/Image/ImageChannels.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/PerPixelViews.h>

namespace vw {

// copy()

template <class ImageT>
class CopyView: public ImageViewBase<CopyView<ImageT>> {
  ImageView<typename ImageT::pixel_type> m_child;
public:
  typedef typename ImageView<typename ImageT::pixel_type>::pixel_type
    pixel_type;
  typedef pixel_type const& result_type;
  typedef typename ImageView<typename ImageT::pixel_type>::pixel_accessor
    pixel_accessor;

  CopyView(ImageT const& image):
    m_child(image.cols(), image.rows(), image.planes()) {
    image.rasterize(m_child, BBox2i(0, 0, image.cols(), image.rows()));
  }

  inline int32 cols() const { return m_child.cols(); }
  inline int32 rows() const { return m_child.rows(); }
  inline int32 planes() const { return m_child.planes(); }

  inline pixel_accessor origin() const { return m_child.origin(); }
  inline result_type operator()(int32 i, int32 j, int32 p = 0) const {
    return m_child(i, j, p);
  }

  typedef CopyView prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& /*bbox*/) const {
    return *this;
  }
  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
    m_child.rasterize(dest, bbox);
  }
};

template <class ImageT>
struct IsMultiplyAccessible<CopyView<ImageT>>: public true_type {};

/// Make a (deep) copy of an image.
template <class ImageT>
CopyView<ImageT> copy(ImageViewBase<ImageT> const& v) {
  return CopyView<ImageT>(v.impl());
}

// crop()

template <class ImageT>
class CropView: public ImageViewBase<CropView<ImageT>> {
  typedef typename boost::mpl::if_<IsFloatingPointIndexable<ImageT>,
    double, int32>::type offset_type;

  ImageT m_child;
  offset_type m_ci, m_cj;
  int32 m_di, m_dj;

public:
  typedef typename ImageT::pixel_type pixel_type;
  typedef typename ImageT::result_type result_type;
  typedef typename ImageT::pixel_accessor pixel_accessor;

  CropView(ImageT const& image, offset_type const upper_left_i,
           offset_type const upper_left_j,
           int32 const width, int32 const height):
    m_child(image), m_ci(upper_left_i), m_cj(upper_left_j),
    m_di(width), m_dj(height) {}

  template<class RealT>
  CropView(ImageT const& image, BBox<RealT, 2> const& bbox):
    m_child(image),
    m_ci((offset_type)(bbox.min()[0])),
    m_cj((offset_type)(bbox.min()[1])),
    m_di(int32(.5 + (bbox.width()))),
    m_dj(int32(.5 + (bbox.height()))) {}

  inline int32 cols() const { return m_di; }
  inline int32 rows() const { return m_dj; }
  inline int32 planes() const { return m_child.planes(); }

  inline pixel_accessor origin() const {
    return m_child.origin().advance(m_ci, m_cj);
  }

  inline result_type operator()(offset_type i, offset_type j,
                                int32 p = 0) const {
    return m_child(m_ci + i, m_cj + j, p);
  }

  CropView const& operator=(CropView const& view) const {
    view.rasterize(*this,
      BBox2i(0, 0, view.impl().cols(), view.impl().rows()));
    return *this;
  }

  template <class ViewT>
  CropView const& operator=(ImageViewBase<ViewT> const& view) const {
    view.impl().rasterize(*this,
      BBox2i(0, 0, view.impl().cols(), view.impl().rows()));
    return *this;
  }

  ImageT const& child() const { return m_child; }

  typedef CropView<typename ImageT::prerasterize_type> prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& bbox) const {
    return prerasterize_type(
      m_child.prerasterize(bbox + Vector2i(m_ci, m_cj)),
      m_ci, m_cj, m_di, m_dj);
  }
  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
};

template <class ImageT>
struct IsFloatingPointIndexable<CropView<ImageT>>:
  public IsFloatingPointIndexable<ImageT> {};

template <class ImageT>
struct IsMultiplyAccessible<CropView<ImageT>>:
  public IsMultiplyAccessible<ImageT> {};

/// Crop an image.
template <class ImageT>
inline CropView<ImageT> crop(ImageViewBase<ImageT> const& v,
                              int32 upper_left_x, int32 upper_left_y,
                              int32 width, int32 height) {
  return CropView<ImageT>(v.impl(), upper_left_x, upper_left_y,
                           width, height);
}

/// Crop an image.
template <class ImageT, class BBoxRealT>
inline CropView<ImageT> crop(ImageViewBase<ImageT> const& v,
                              BBox<BBoxRealT, 2> const& bbox) {
  return CropView<ImageT>(v.impl(), bbox);
}

// subsample()

template <class ChildT>
class SubsamplePixelAccessor {
  ChildT m_child;
  int32 m_xdelta, m_ydelta;
public:
  typedef typename ChildT::pixel_type pixel_type;
  typedef typename ChildT::result_type result_type;
  typedef typename ChildT::offset_type offset_type;
  SubsamplePixelAccessor(ChildT const& acc, int32 xdelta, int32 ydelta):
    m_child(acc), m_xdelta(xdelta), m_ydelta(ydelta) {}

  inline SubsamplePixelAccessor& next_col() {
    m_child.advance(m_xdelta, 0); return *this;
  }
  inline SubsamplePixelAccessor& prev_col() {
    m_child.advance(-m_xdelta, 0); return *this;
  }
  inline SubsamplePixelAccessor& next_row() {
    m_child.advance(0, m_ydelta); return *this;
  }
  inline SubsamplePixelAccessor& prev_row() {
    m_child.advance(0, -m_ydelta); return *this;
  }
  inline SubsamplePixelAccessor& next_plane() {
    m_child.next_plane(); return *this;
  }
  inline SubsamplePixelAccessor& prev_plane() {
    m_child.prev_plane(); return *this;
  }
  inline SubsamplePixelAccessor& advance(offset_type di, offset_type dj,
                                          ssize_t dp = 0) {
    m_child.advance((offset_type)m_xdelta * di,
                    (offset_type)m_ydelta * dj, dp);
    return *this;
  }

  inline result_type operator*() const { return *m_child; }
};

template <class ImageT>
class SubsampleView: public ImageViewBase<SubsampleView<ImageT>> {
  ImageT m_child;
  int32 m_xdelta, m_ydelta;
public:
  typedef typename ImageT::pixel_type pixel_type;
  typedef typename ImageT::result_type result_type;
  typedef SubsamplePixelAccessor<typename ImageT::pixel_accessor>
    pixel_accessor;

  SubsampleView(ImageT const& image, int32 subsampling_factor):
    m_child(image), m_xdelta(subsampling_factor),
    m_ydelta(subsampling_factor) {
    VW_ASSERT(m_xdelta > 0 && m_ydelta > 0,
              ArgumentErr()
              << "SubsampleView: Arguments must be greater than zero.");
  }
  SubsampleView(ImageT const& image, int32 xfactor, int32 yfactor):
    m_child(image), m_xdelta(xfactor), m_ydelta(yfactor) {
    VW_ASSERT(m_xdelta > 0 && m_ydelta > 0,
              ArgumentErr()
              << "SubsampleView: Arguments must be greater than zero.");
  }

  inline int32 cols() const {
    return 1 + (m_child.cols() - 1) / m_xdelta;
  }
  inline int32 rows() const {
    return 1 + (m_child.rows() - 1) / m_ydelta;
  }
  inline int32 planes() const { return m_child.planes(); }

  inline pixel_accessor origin() const {
    return pixel_accessor(m_child.origin(), m_xdelta, m_ydelta);
  }
  inline result_type operator()(int32 i, int32 j, int32 p = 0) const {
    return m_child(m_xdelta * i, m_ydelta * j, p);
  }

  ImageT const& child() const { return m_child; }

  typedef CropView<ImageView<pixel_type>> prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& bbox) const {
    ImageView<pixel_type> buffer(bbox.width(), bbox.height());

    Vector2i sub_region_size(bbox.width() / m_xdelta,
                             bbox.height() / m_ydelta);
    if (sub_region_size.x() < 2)
      sub_region_size.x() = 2;
    if (sub_region_size.y() < 2)
      sub_region_size.y() = 2;

    typedef std::vector<BBox2i> ContainerT;
    ContainerT bboxes = subdivide_bbox(bbox, sub_region_size.x(),
                                       sub_region_size.y());

    typedef SubsampleView<typename ImageT::prerasterize_type>
      input_pre_type;

    for (ContainerT::const_iterator b = bboxes.begin();
         b != bboxes.end(); ++b) {
      vw::rasterize(
        input_pre_type(
          m_child.prerasterize(
            BBox2i(m_xdelta * (*b).min().x(),
                   m_ydelta * (*b).min().y(),
                   m_xdelta * ((*b).width() - 1) + 1,
                   m_xdelta * ((*b).height() - 1) + 1)),
          m_xdelta, m_ydelta),
        crop(buffer, *b - bbox.min()), *b);
    }

    return crop(buffer, -bbox.min().x(), -bbox.min().y(),
                cols(), rows());
  }
  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
};

template <class ImageT>
struct IsMultiplyAccessible<SubsampleView<ImageT>>:
  public IsMultiplyAccessible<ImageT> {};

/// Subsample an image by an integer factor.
template <class ImageT>
inline SubsampleView<ImageT> subsample(ImageT const& v,
                                        int32 subsampling_factor) {
  return SubsampleView<ImageT>(v, subsampling_factor);
}

/// Subsample an image by integer factors in x and y.
template <class ImageT>
inline SubsampleView<ImageT> subsample(ImageT const& v,
                                        int32 xfactor, int32 yfactor) {
  return SubsampleView<ImageT>(v, xfactor, yfactor);
}

} // namespace vw

#endif // __VW_IMAGE_MANIPULATION_H__
