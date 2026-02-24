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

/// \file AlgorithmFunctions.h
///
/// Core non-lazy algorithms: fill, bounding_box, subdivide_bbox.
/// Opacity queries (is_opaque, is_transparent, nonzero_data_bounding_box)
/// are in ImageOpacity.h.

#ifndef __VW_IMAGE_ALGORITHM_FUNCTIONS_H__
#define __VW_IMAGE_ALGORITHM_FUNCTIONS_H__

#include <vw/Image/ImageView.h>
#include <vw/Image/MaskViews.h>

namespace vw {

// fill()

/// Fill an image with a constant pixel value
template <class ImageT, class ValueT>
void fill(ImageViewBase<ImageT> const& image, ValueT value_) {
  int32 planes = image.impl().planes();
  int32 rows = image.impl().rows();
  int32 cols = image.impl().cols();
  typename ImageT::pixel_type value(value_);
  typename ImageT::pixel_accessor plane = image.impl().origin();
  for (int32 p = planes; p; --p) {
    typename ImageT::pixel_accessor row = plane;
    for (int32 r = rows; r; --r) {
      typename ImageT::pixel_accessor col = row;
      for (int32 c = cols; c; --c) {
        *col = value;
        col.next_col();
      }
      row.next_row();
    }
    plane.next_plane();
  }
}

// bounding_box()

template <class ViewT>
BBox2i bounding_box(ImageViewBase<ViewT> const& image_) {
  return BBox2i(0, 0, image_.impl().cols(), image_.impl().rows());
}

// subdivide_bbox()

/// A utility routine that, given an image, returns a vector of
/// bounding boxes for sub-regions of the image of the specified
/// size. Note that bounding boxes along the right and bottom edges
/// of the image will not have the specified dimension unless the
/// image width and height are perfectly divisible by the bounding
/// box width and height, respectively.
/// - If include_partials is set to false, only full size boxes will be
///   included in the output.
/// - If full_size is true, the boxes at the right and bottom edges will
///   be grown inward to have full size (they will overlap with earlier
///   boxes). This assumes include_partials = true.
/// - Output tiles are in raster order, top left to bottom right.
std::vector<BBox2i>
subdivide_bbox(BBox2i const& object, int32 block_width, int32 block_height,
               bool include_partials = true, bool full_size = false);

template <class T>
inline std::vector<BBox2i>
subdivide_bbox(ImageViewBase<T> const& view,
               int32 block_width, int32 block_height,
               bool include_partials = true) {
  return subdivide_bbox(bounding_box(view.impl()),
                        block_width, block_height, include_partials);
}

} // namespace vw

#endif // __VW_IMAGE_ALGORITHM_FUNCTIONS_H__
