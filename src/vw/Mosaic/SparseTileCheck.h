// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__

/// \file SparseTileCheck.h
/// 
/// 
#ifndef __VW_MOSAIC_SPARSE_TILE_CHECK_H__
#define __VW_MOSAIC_SPARSE_TILE_CHECK_H__

#include <vw/Math/BBox.h>

namespace vw {
namespace mosaic {

  // This is the default version for most image views.  It simply
  // checks for intersection with the entire bounding box of the
  // image.
  template <class SrcViewT>
  class SparseTileCheck {
    BBox2i m_src_bbox;
  public:
    SparseTileCheck(SrcViewT const& source) : m_src_bbox(0,0,source.cols(),source.rows()) {}
    bool operator() (BBox2i const& bbox) { return bbox.intersects(m_src_bbox); }
  };
  
}} // namespace vw::mosaic

#endif // __VW_MOSAIC_SPARSE_TILE_CHECK_H__
