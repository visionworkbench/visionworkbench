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


/// \file BlockImageOperator.h
///
/// Execute a function on an entire image one block at a time.
/// This is similar to BlockRasterizeView except that it does 
/// not provide a view interface.  The results must be handled
/// by the functor object since this class provides no output.
///
/// FuncT must implement the following function:
///
/// template <class T>
/// void operator()(ImageView<T> const& image, BBox2i const& bbox) {
/// --> Note that bbox is the region that image was cropped from!
///


#ifndef __VW_IMAGE_BLOCKIMAGEOPERATOR_H__
#define __VW_IMAGE_BLOCKIMAGEOPERATOR_H__

#include <vw/Core/Cache.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/PixelAccessors.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/BlockProcessor.h>

namespace vw {

  /// A wrapper view that processes the input image in blocks.
  template <class ImageT, class FuncT>
  class BlockImageOperator {
  public:
    typedef typename ImageT::pixel_type pixel_type;

    BlockImageOperator( ImageT const& image, FuncT & functor,
                        Vector2i const& block_size,
                        int num_threads = 0, Cache *cache = NULL )
      : m_child           ( new ImageT(image) ),
        m_functor         ( functor ),
        m_block_size      ( block_size ),
        m_num_threads     ( num_threads ),
        m_cache_ptr       ( cache )
    {
      if( m_block_size.x() <= 0 || m_block_size.y() <= 0 )
        m_block_size = image_block::get_default_block_size<pixel_type>(image.rows(), image.cols(), image.planes());

      if (m_cache_ptr) // Manager is not needed if not using a cache.
        m_block_manager.initialize(m_cache_ptr, m_block_size, m_child);
    }


    /// Apply the functor to the specified section of the image.
    /// - This function can be run on a bbox of any size, including the whole image.
    void execute_func(BBox2i const& bbox) {

      // Create functor which will call the input functor on individual image blocks.
      RunFunctor functor_wrapper(*this);

      // Set up block processor to call the functor in parallel blocks.
      image_block::BlockProcessor<RunFunctor> processor( functor_wrapper, m_block_size, m_num_threads );

      // Tell the block processor to do all the work.
      processor(bbox);
    }

  private:

    // These functions are here to be convenient for other private functions.
    ImageT      & child()       { return *m_child; }
    ImageT const& child() const { return *m_child; }

    // These function objects are spawned to rasterize the child image.
    // One functor is created per child thread, and they are called
    // in succession with bounding boxes that are each contained within one block.
    class RunFunctor {
      BlockImageOperator const& m_view;
    public:
      RunFunctor( BlockImageOperator const& view)
        : m_view(view) {}

      /// Rasterize part of m_view into a buffer, then call m_functor on it.
      void operator()( BBox2i const& bbox ) const {

        // Set up the buffer where the image data will be copied to.
        ImageView<BlockImageOperator::pixel_type> temp_image( bbox.width(), bbox.height(), 
                                                              m_view.child().planes() );

        if( m_view.m_cache_ptr ) {
          // Ask the cache managing object to get the image tile, we might already have it.
          Vector2i block_index = m_view.m_block_manager.get_block_index(bbox);

          const Cache::Handle<image_block::BlockGenerator<ImageT> >& handle
            = m_view.m_block_manager.block(block_index);
          handle->rasterize( temp_image, bbox - m_view.m_block_manager.get_block_start_pixel(block_index));
          handle.release();
        }
        else // No cache, generate the image tile from scratch.
          m_view.child().rasterize(temp_image, bbox);

        // Now that we have rasterized the required portion of the image, call our functor on it.
        m_view.m_functor(temp_image, bbox);
      }
    }; // End class RunFunctor

    // Allows RasterizeFunctor to access cache-related members.
    friend class RunFunctor;

    // We store this by shared pointer so it doesn't move when we copy
    // the BlockImageOperator, since the BlockGenerators point to it.
    boost::shared_ptr<ImageT> m_child;
    FuncT    &m_functor;
    Vector2i  m_block_size;
    int32     m_num_threads;
    Cache    *m_cache_ptr;

    /// This object keeps track of the BlockGenerator for each image tile (if using a cache)
    image_block::BlockGeneratorManager<ImageT> m_block_manager;

  }; // End class BlockImageOperator


  /// Create a BlockRasterizeView with no caching.
  /// - If each tile will only be accessed once then there is no need for a cache!
  template <class ImageT, class FuncT>
  inline void block_op( ImageViewBase<ImageT> const& image,
                        FuncT & functor,
                        Vector2i const& block_size, int num_threads = 0 ) {
    BlockImageOperator<ImageT, FuncT> block_op(image.impl(), functor, block_size, num_threads);
    block_op.execute_func(bounding_box(image));
  }

  /// Create a BlockImageOperator using the vw system Cache object.
  /// - Use this if there will be repeat tile access.
  template <class ImageT, class FuncT>
  inline void block_op_cache( ImageViewBase<ImageT> const& image,
                              FuncT & functor,
                              Vector2i const& block_size, int num_threads = 0 ) {

    BlockImageOperator<ImageT, FuncT> block_op(image.impl(), functor, block_size,
                                               num_threads, &vw_system_cache());
    block_op.execute_func(bounding_box(image));
  }

} // namespace vw

#endif // __VW_IMAGE_BLOCKIMAGEOPERATOR_H__
