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

/// \file BlockCacheImageView.h
///
/// An image view that rasterizes sub blocks of the image and keeps
/// them around in the cache as long as possible.  This is often
/// useful if you are converting between scanline rasterization and
/// block rasterization and vice versa.  For example, if you are
/// processing a viwe that is most efficiently rasterized using a
/// block based approach, but you want to write the results using a
/// scanline based file driver, you can introduce this class as an
/// intermediary.  It will cache the blocks that span the width of the
/// image (as long as you cache is large enough to hold them all) and
/// service the scanline based rasterize requested by the file driver
/// from the cached data.
///
#ifndef __VW_IMAGE_BLOCKCACHEIMAGEVIEW_H__
#define __VW_IMAGE_BLOCKCACHEIMAGEVIEW_H__

#include <string>
#include <map>

#include <vw/Core/Cache.h>
#include <vw/Image/ImageResourceView.h>
#include <vw/Image/ViewImageResource.h>
#include <vw/Image/ImageViewBase.h>

namespace vw {

  /// A view of an image on disk.
  template <class PixelT>
  class BlockCacheView : public ImageResourceView<PixelT>
  {
    typedef ImageResourceView<PixelT> base_type;

  public:
    /// Constructs a BlockCacheView of the given input view
    template <class ViewT>
    BlockCacheView(ImageViewBase<ViewT> const& view, Vector2i block_size, bool cache=true )
      : base_type( new ViewImageResource( view.impl(), block_size ), cache ) {}

    /// Constructs a BlockCacheView of the given input view using
    /// the specified cache area.
    template <class ViewT>
    BlockCacheView(ImageViewBase<ViewT> const& view, Vector2i block_size, Cache& cache )
      : base_type( new ViewImageResource( view.impl(), block_size ), cache ) {}

    /// Constructs a BlockCacheView of the given resource.  Takes
    /// ownership of the resource object (i.e. deletes it when it's
    /// done using it).
    BlockCacheView( ViewImageResource *resource, bool cache=true )
      : base_type( resource, cache ) {}

    /// Constructs a BlockCacheView of the given resource using the
    /// specified cache area.  Takes ownership of the resource object
    /// (i.e. deletes it when it's done using it).
    BlockCacheView( ViewImageResource *resource, Cache& cache )
      : base_type( resource, cache ) {}

    virtual ~BlockCacheView() {}

    /// Returns the cache block size currently in use.
    Vector2i block_size() const { return this->resource()->native_block_size(); }
  };

} // namespace vw

#endif // __VW_IMAGE_BLOCKCACHEVIEW_H__
