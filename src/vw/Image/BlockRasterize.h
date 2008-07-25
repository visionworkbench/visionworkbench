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

/// \file BlockRasterizeView.h
/// 
/// Defines a wrapper view that causes an image view to be rasterized
/// one block at a time, rather than all at once.  In some situations
/// this can dramatically improve performance by reducing memory
/// utilization.
///
#ifndef __VW_IMAGE_BLOCKRASTERIZE_H__
#define __VW_IMAGE_BLOCKRASTERIZE_H__

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/BlockProcessor.h>

namespace vw {

  /// A wrapper view that rasterizes its child in blocks.
  template <class ImageT>
  class BlockRasterizeView : public ImageViewBase<BlockRasterizeView<ImageT> > {
    ImageT m_child;
    int32 m_block_cols, m_block_rows;
    int m_num_threads;
  public:
    typedef typename ImageT::pixel_type pixel_type;
    typedef typename ImageT::result_type result_type;
    typedef typename ImageT::pixel_accessor pixel_accessor;

    BlockRasterizeView( ImageT const& image,
                        int32 block_cols = 0,
                        int32 block_rows = 0,
                        int num_threads = 0 )
      : m_child(image),
        m_block_cols(block_cols),
        m_block_rows(block_rows),
        m_num_threads(num_threads?num_threads:(Thread::default_num_threads()))
    {
      if( m_block_cols <= 0 ) m_block_cols = image.cols();
      if( m_block_rows <= 0 ) m_block_rows = (std::max)( image.rows()/m_num_threads, 1 );
      // XXX Should the default block size be different for very wide images?
    }

    inline int32 cols() const { return m_child.cols(); }
    inline int32 rows() const { return m_child.rows(); }
    inline int32 planes() const { return m_child.planes(); }

    inline pixel_accessor origin() const { return m_child.origin(); }
    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_child(i,j,p); }

    ImageT& child() { return m_child; }
    ImageT const& child() const { return m_child; }

    /// \cond INTERNAL
    template <class DestT>
    class RasterizeFunctor {
      ImageT const& m_src;
      DestT const& m_dest;
      Vector2i m_offset;
    public:
      RasterizeFunctor( ImageT const& src, DestT const& dest, Vector2i const& offset )
        : m_src(src), m_dest(dest), m_offset(offset) {}
      void operator()( BBox2i const& bbox ) const {
        m_src.rasterize( crop( m_dest, bbox+m_offset ), bbox );
      }
    };

    typedef BlockRasterizeView prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i bbox ) const { return *this; }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const {
      RasterizeFunctor<DestT> rasterizer(m_child,dest,bbox.min());
      BlockProcessor<RasterizeFunctor<DestT> > process(rasterizer, m_block_cols, m_block_rows);
      process(bbox);
    }
    /// \endcond
  };
  
  /// \cond INTERNAL
  template <class ImageT>
  struct IsMultiplyAccessible<BlockRasterizeView<ImageT> > : IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// A free function that creates a BlockRasterizeView wrapping the 
  /// given view, causing it to rasterize in blocks rather than all 
  /// at once.
  template <class ImageT>
  inline BlockRasterizeView<ImageT> block_rasterize( ImageViewBase<ImageT> const& image,
                                                     int32 block_cols = 0,
                                                     int32 block_rows = 0,
                                                     int num_threads = 0 ) {
    return BlockRasterizeView<ImageT>( image.impl(), block_cols, block_rows, num_threads );
  }

} // namespace vw

#endif // __VW_IMAGE_BLOCKRASTERIZE_H__
