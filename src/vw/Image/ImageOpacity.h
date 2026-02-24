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

/// \file ImageOpacity.h
///
/// Image-level opacity/transparency queries and nonzero bounding box.
/// Split from AlgorithmFunctions.h for compile-time reduction.

#ifndef __VW_IMAGE_IMAGEOPACITY_H__
#define __VW_IMAGE_IMAGEOPACITY_H__

#include <vw/Image/ImageView.h>

namespace vw {

// is_opaque()

template <class ImageT>
bool is_opaque_helper(ImageT const& image, true_type) {
  for (int32 y = 0; y < image.rows(); ++y)
    for (int32 x = 0; x < image.cols(); ++x)
      if (!(is_opaque(image(x, y))))
        return false;
  return true;
}

template <class ImageT>
bool is_opaque_helper(ImageT const& /*image*/, false_type) {
  return true;
}

/// Returns true if the given image is entirely opaque, or false if
/// it is at least partially transparent.
template <class ImageT>
bool is_opaque(ImageViewBase<ImageT> const& image) {
  return is_opaque_helper(image.impl(),
    typename PixelHasAlpha<typename ImageT::pixel_type>::type());
}

// is_transparent()

template <class ImageT>
bool is_transparent_helper(ImageT const& image, true_type) {
  for (int32 y = 0; y < image.rows(); ++y)
    for (int32 x = 0; x < image.cols(); ++x)
      if (!is_transparent(image(x, y)))
        return false;
  return true;
}

template <class ImageT>
bool is_transparent_helper(ImageT const& /*image*/, false_type) {
  return false;
}

/// Returns true if the given image is entirely transparent, or false if
/// it is opaque or only partially transparent.
template <class ImageT>
bool is_transparent(ImageViewBase<ImageT> const& image) {
  return is_transparent_helper(image.impl(),
    typename PixelHasAlpha<typename ImageT::pixel_type>::type());
}

// nonzero_data_bounding_box()

template <class ViewT>
BBox2i nonzero_data_bounding_box(ImageViewBase<ViewT> const& image_) {
  const typename ViewT::pixel_type zero = typename ViewT::pixel_type();
  ViewT const& image = static_cast<ViewT const&>(image_);
  int32 x = 0, y = 0, cols = 0, rows = 0;
  int32 i, j, icols = image.cols(), irows = image.rows();
  for (j = 0; j < irows; ++j) {
    for (i = 0; i < icols; ++i)
      if (image(i, j) != zero)
        break;
    if (i != icols)
      break;
  }
  if (j != irows) {
    y = j;
    for (j = irows - 1; j; --j) {
      for (i = 0; i < icols; ++i)
        if (image(i, j) != zero)
          break;
      if (i != icols)
        break;
    }
    rows = j - y + 1;
    for (i = 0; i < icols; ++i) {
      for (j = y; j < y + rows; ++j)
        if (image(i, j) != zero)
          break;
      if (j != y + rows)
        break;
    }
    x = i;
    for (i = icols - 1; i; --i) {
      for (j = y; j < y + rows; ++j)
        if (image(i, j) != zero)
          break;
      if (j != y + rows)
        break;
    }
    cols = i - x + 1;
  }
  return BBox2i(x, y, cols, rows);
}

} // namespace vw

#endif // __VW_IMAGE_IMAGEOPACITY_H__
