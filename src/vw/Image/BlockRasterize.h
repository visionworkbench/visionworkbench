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

namespace vw {

  /// A wrapper view that rasterizes its child in blocks.
  template <class ImageT>
  class BlockRasterizeView : public ImageViewBase<BlockRasterizeView<ImageT> > {
    ImageT m_child;
    int m_blocksize;
  public:
    typedef typename ImageT::pixel_type pixel_type;
    typedef typename ImageT::result_type result_type;
    typedef typename ImageT::pixel_accessor pixel_accessor;

    BlockRasterizeView( ImageT const& image, int blocksize = 256 ) : m_child(image), m_blocksize(blocksize) {}

    inline int blocksize() const { return m_blocksize; }
    inline void set_blocksize( int blocksize ) { m_blocksize = blocksize; }

    inline unsigned cols() const { return m_child.cols(); }
    inline unsigned rows() const { return m_child.rows(); }
    inline unsigned planes() const { return m_child.planes(); }

    inline pixel_accessor origin() const { return m_child.origin(); }
    inline result_type operator()( int i, int j, int p=0 ) const { return m_child(i,j,p); }

    ImageT& child() { return m_child; }
    ImageT const& child() const { return m_child; }

    /// \cond INTERNAL
    typedef BlockRasterizeView prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i bbox ) const { return *this; }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const {
      for( int y=0; y<bbox.height(); y+=m_blocksize ) {
        int height = (std::min)( m_blocksize, bbox.height()-y );
        for( int x=0; x<bbox.width(); x+=m_blocksize ) {
          int width = (std::min)( m_blocksize, bbox.width()-x );
          BBox2i block_bbox( x+bbox.min().x(), y+bbox.min().y(), width, height );
          m_child.prerasterize( block_bbox ).rasterize( crop( dest, x, y, width, height ), block_bbox );
        }
      }
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
  inline BlockRasterizeView<ImageT> block_rasterize( ImageT const& image, int blocksize = 256 ) {
    return BlockRasterizeView<ImageT>( image, blocksize );
  }

} // namespace vw

#endif // __VW_IMAGE_BLOCKRASTERIZE_H__
