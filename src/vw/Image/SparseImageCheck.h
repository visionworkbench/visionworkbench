// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
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


/// \file SparseImageCheck.h
///
///
#ifndef __VW_IMAGE_SPARSE_IMAGE_CHECK_H__
#define __VW_IMAGE_SPARSE_IMAGE_CHECK_H__

#include <vw/Math/BBox.h>

namespace vw {

  /// Checks to see if a given bounding box intersects the wrapped image.
  /// - This is the default version for most image views, but some classes
  ///   should specialize the class to work differently when it wraps them.
  template <class SrcViewT>
  class SparseImageCheck {
    BBox2i m_src_bbox;
    
  public:
    /// Constructor wraps another image object.
    SparseImageCheck(SrcViewT const& source) : m_src_bbox(0,0,source.cols(),source.rows()) {}
    
    /// Return true if the given bounding box intersects the image data.
    bool operator() (BBox2i const& bbox) { return bbox.intersects(m_src_bbox); }
  };

  // A helper function to make it easy to test for sparisty.
  template <class ImageT>
    bool sparse_check( ImageT const& image, BBox2i const& bbox ) {
    return SparseImageCheck<ImageT>(image)(bbox);
  }

} // namespace vw

#endif // __VW_IMAGE_SPARSE_IMAGE_CHECK_H__
