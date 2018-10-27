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


/// \file BlockRasterizeView.h
///
/// Defines a wrapper view that causes an image view to be rasterized
/// one block at a time, rather than all at once, and supports multi-
/// threaded rasterization as well as block-level cacheing.  Even when
/// you don't need those features, in some situations rasterizing one
/// block at a time can dramatically improve performance by reducing
/// memory utilization.
///
#ifndef __VW_IMAGE_BLOCKRASTERIZE_H__
#define __VW_IMAGE_BLOCKRASTERIZE_H__

#include <vw/Core/Cache.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/PixelAccessors.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/BlockProcessor.h>

namespace vw {

  /// A wrapper view that rasterizes its child in blocks.
  template <class ImageT>
  class BlockRasterizeView : public ImageViewBase<BlockRasterizeView<ImageT> > {
  public:
    typedef typename ImageT::pixel_type pixel_type;
    typedef typename ImageT::pixel_type result_type;
    typedef ProceduralPixelAccessor<BlockRasterizeView> pixel_accessor;

    BlockRasterizeView( ImageT const& image, Vector2i const block_size,
                        int num_threads = 0, Cache *cache = NULL )
      : m_child           ( new ImageT(image) ),
        m_block_size      ( block_size ),
        m_num_threads     ( num_threads ),
        m_cache_ptr       ( cache )
    {
      if( m_block_size.x() <= 0 || m_block_size.y() <= 0 )
        m_block_size = image_block::get_default_block_size<pixel_type>(image.rows(), image.cols(), image.planes());

      if (m_cache_ptr) // Manager is not needed if not using a cache.
        m_block_manager.initialize(m_cache_ptr, m_block_size, m_child);
    }

    inline int32 cols  () const { return m_child->cols();   }
    inline int32 rows  () const { return m_child->rows();   }
    inline int32 planes() const { return m_child->planes(); }
    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0, 0 ); }

    inline result_type operator()( int32 x, int32 y, int32 p = 0 ) const {
#if VW_DEBUG_LEVEL > 1
      VW_OUT(VerboseDebugMessage, "image") << "ImageResourceView rasterizing pixel (" 
                                           << x << "," << y << "," << p << ")" << std::endl;
#endif
      if ( m_cache_ptr ) {
        // Note that requesting a value from a handle forces that data to be generated.
        // Early-out optimization for single-block resources
        if( m_block_manager.only_one_block() ) {
          const Cache::Handle<image_block::BlockGenerator<ImageT> >& handle = m_block_manager.quick_single_block();
          result_type result = handle->operator()( x, y, p );
          handle.release();
          return result;
        }

        Vector2i block_index = m_block_manager.get_block_index(Vector2i(x,y));
        const Cache::Handle<image_block::BlockGenerator<ImageT> >& handle
            = m_block_manager.block(block_index);
        Vector2i start_pixel = m_block_manager.get_block_start_pixel(block_index);
        result_type result = handle->operator()( x - start_pixel.x(),
                                                 y - start_pixel.y(), p );
        handle.release(); // Let the Cache know we are finished with this block of data.
        return result;
      }
      else // No cache, access the image directly.  This could be brutally slow!
        return (*m_child)(x,y,p);
    }

    ImageT      & child()       { return *m_child; }
    ImageT const& child() const { return *m_child; }

    typedef CropView<ImageView<pixel_type> > prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      // Init output data
      ImageView<pixel_type> buf( bbox.width(), bbox.height(), planes() );
      // Fill in the output data from this view
      rasterize( buf, bbox );
      // "Fake" the bbox image so it looks like a full size image.
      return CropView<ImageView<pixel_type> >( buf, BBox2i(-bbox.min().x(),-bbox.min().y(),cols(),rows()) );
    }

    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
      // Create functor to rasterize this image into the destination image
      RasterizeFunctor<DestT> rasterizer( *this, dest, bbox.min() );
      // Set up block processor to call the functor in parallel blocks.
      image_block::BlockProcessor<RasterizeFunctor<DestT> > process( rasterizer, m_block_size, m_num_threads );
      // Tell the block processor to do all the work.
      process(bbox);
    }

  private:
    // These function objects are spawned to rasterize the child image.
    // One functor is created per child thread, and they are called
    // in succession with bounding boxes that are each contained within one block.
    template <class DestT>
    class RasterizeFunctor {
      BlockRasterizeView const& m_view;
      DestT const& m_dest;
      Vector2i     m_offset;
    public:
      RasterizeFunctor( BlockRasterizeView const& view, DestT const& dest, Vector2i const& offset )
        : m_view(view), m_dest(dest), m_offset(offset) {}

      /// Rasterize part of m_view into m_dest.
      void operator()( BBox2i const& bbox ) const {
#if VW_DEBUG_LEVEL > 1
        VW_OUT(VerboseDebugMessage, "image") << "BlockRasterizeView::RasterizeFunctor( " << bbox << " )" << std::endl;
#endif
        if( m_view.m_cache_ptr ) {
          // Ask the cache managing object to get the image tile, we might already have it.
          Vector2i block_index = m_view.m_block_manager.get_block_index(bbox);

          const Cache::Handle<image_block::BlockGenerator<ImageT> >& handle
            = m_view.m_block_manager.block(block_index);
          handle->rasterize( crop( m_dest, bbox-m_offset ),
                             bbox - m_view.m_block_manager.get_block_start_pixel(block_index) );
          handle.release();
        }
        else // No cache, generate the image tile from scratch.
          m_view.child().rasterize( crop( m_dest, bbox-m_offset ), bbox );
      }
    }; // End class RasterizeFunctor

    // Allows RasterizeFunctor to access cache-related members.
    template <class DestT> friend class RasterizeFunctor;

    // We store this by shared pointer so it doesn't move when we copy
    // the BlockRasterizeView, since the BlockGenerators point to it.
    boost::shared_ptr<ImageT> m_child;
    Vector2i m_block_size;
    int32    m_num_threads;
    Cache   *m_cache_ptr;

    /// This object keeps track of the BlockGenerator for each image tile (if using a cache)
    image_block::BlockGeneratorManager<ImageT> m_block_manager;
  };

  /// Create a BlockRasterizeView with no caching.
  template <class ImageT>
  inline BlockRasterizeView<ImageT> block_rasterize( ImageViewBase<ImageT> const& image,
                                                     Vector2i const& block_size, int num_threads = 0 ) {
    return BlockRasterizeView<ImageT>( image.impl(), block_size, num_threads );
  }

  /// Create a BlockRasterizeView using the vw system Cache object.
  template <class ImageT>
  inline BlockRasterizeView<ImageT> block_cache( ImageViewBase<ImageT> const& image,
                                                 Vector2i const& block_size, int num_threads = 0 ) {
    return BlockRasterizeView<ImageT>( image.impl(), block_size, num_threads, &vw_system_cache() );
  }

  /// Create a BlockRasterizeView using the provided Cache object.
  template <class ImageT>
  inline BlockRasterizeView<ImageT> block_cache( ImageViewBase<ImageT> const& image,
                                                 Vector2i const& block_size, int num_threads, Cache& cache ) {
    return BlockRasterizeView<ImageT>( image.impl(), block_size, num_threads, &cache );
  }

} // namespace vw

#endif // __VW_IMAGE_BLOCKRASTERIZE_H__
