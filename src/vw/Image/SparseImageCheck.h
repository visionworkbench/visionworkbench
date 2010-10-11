// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file SparseImageCheck.h
///
///
#ifndef __VW_IMAGE_SPARSE_IMAGE_CHECK_H__
#define __VW_IMAGE_SPARSE_IMAGE_CHECK_H__

#include <vw/Math/BBox.h>

namespace vw {

  // This is the default version for most image views.  It simply
  // checks for intersection with the entire bounding box of the
  // image.
  template <class SrcViewT>
  class SparseImageCheck {
    BBox2i m_src_bbox;
  public:
    SparseImageCheck(SrcViewT const& source) : m_src_bbox(0,0,source.cols(),source.rows()) {}
    bool operator() (BBox2i const& bbox) { return bbox.intersects(m_src_bbox); }
  };

  // A helper function to make it easy to test for sparisty.
  template <class ImageT>
    bool sparse_check( ImageT const& image, BBox2i const& bbox ) {
    return SparseImageCheck<ImageT>(image)(bbox);
  }

} // namespace vw

#endif // __VW_IMAGE_SPARSE_IMAGE_CHECK_H__
