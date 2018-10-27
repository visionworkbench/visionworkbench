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


/// \file BlockProcessor.h
///
/// This class provides a bit of infrastructure to make it easier to
/// implement algorithms in a block-oriented way.  You specify an
/// arbitrary function, a grid spacing, and a maximum number of
/// processing threads.  You can then call the block processor,
/// passing it an arbitrarily large bounding box.  It will chop that
/// bounding box up into blocks and call the callback function on
/// each block, spawning as many child threads as you request.
///
/// Strictly speaking, this doesn't need to be in the Image module.
/// However, it was designed for large image processing, it depends
/// on Math but does not really belong there, and there's nothing
/// between there and here in the dependency chain.  So, here it is.
///
#ifndef __VW_IMAGE_BLOCKPROCESSOR_H__
#define __VW_IMAGE_BLOCKPROCESSOR_H__

#include <vw/Core/Settings.h>
#include <vw/Core/Thread.h>
#include <vw/Math/BBox.h>

namespace vw {

/// These things require careful use and are put in a namespace to keep 
///  them from being accidentally used.
namespace image_block {
  
  /// Class to call m_func(BBox2i) in parallel.  It is up to the FuncT
  ///  type to handle what that should do.
  template <class FuncT>
  class BlockProcessor {
    FuncT    m_func;
    Vector2i m_block_size;
    uint32   m_num_threads;
  public:

    /// Create a BlockProcessor object with the specified parameters.
    /// - The function will get executed in "units" of block_size simultaneously
    ///   by the specified number of threads.
    /// - The func object must have an operator(BBox2i) function that does whatever
    ///   it is you want done.
    BlockProcessor( FuncT const& func, Vector2i const& block_size, uint32 threads = 0 )
      : m_func(func), m_block_size(block_size),
        m_num_threads(threads?threads:(vw_settings().default_num_threads())) {}

    /// We will construct and call one BlockThread per thread.
    class BlockThread {
    public:
      // All the BlockThread objects share a reference to a shared Info object,
      // which stores information about what block should be processed next.
      class Info {
      public:
        Info( FuncT const& func, BBox2i const& total_bbox, Vector2i const& block_size )
          : m_func(func), m_total_bbox(total_bbox),
            m_block_bbox(round_down(total_bbox.min().x(),block_size.x()),round_down(total_bbox.min().y(),block_size.y()),block_size.x(),block_size.y()),
            m_block_size(block_size) {}

        // Return the next block bbox to process.
        BBox2i bbox() const {
          BBox2i block_bbox = m_block_bbox;
          block_bbox.crop( m_total_bbox );
          return block_bbox;
        }

        // Return the processing function.
        FuncT const& func() const {
          return m_func;
        }

        // Are we finished?
        bool complete() const {
          return ( m_block_bbox.min().y() >= m_total_bbox.max().y() );
        }

        // Returns the info mutex, for locking.
        Mutex& mutex() {
          return m_mutex;
        }

        // Advance the block_bbox to point to the next block to process.
        void advance() {
          m_block_bbox.min().x() += m_block_size.x();
          if( m_block_bbox.min().x() >= m_total_bbox.max().x() ) {
            m_block_bbox.min().x() = round_down(m_total_bbox.min().x(),m_block_size.x());
            m_block_bbox.min().y() += m_block_size.y();
            m_block_bbox.max().y() = m_block_bbox.min().y() + m_block_size.y();
          }
          m_block_bbox.max().x() = m_block_bbox.min().x() + m_block_size.x();
        }

      private:
        // This hideous nonsense rounds an integer value *down* to the nearest
        // multple of the given modulus.  It's this hideous partly because
        // it avoids modular arithematic on negative numbers, which is technically
        // implementation-defined in all but the most recent C/C++ standards.
        static int32 round_down(int32 val, int32 mod) {
          return val + ((val>=0) ? (-(val%mod)) : (((-val-1)%mod)-mod+1));
        }

        FuncT const& m_func;
        BBox2i   m_total_bbox, m_block_bbox;
        Vector2i m_block_size;
        Mutex    m_mutex;
      }; // End class Info

      BlockThread( Info &info ) : info(info) {}

      void operator()() {
        while( true ) {
          BBox2i bbox;
          {
            // Grab the next bbox to process, and update it with the
            // subsequent bbox for the next thread to grab.
            Mutex::Lock lock(info.mutex());
            if( info.complete() )
              return;
            bbox = info.bbox();
            info.advance();
          }
          info.func()( bbox );
        }
      }

    private:
      Info &info;
    }; // End class BlockThread

    /// Break bbox into sections of block_size, then call
    ///  func(sub_bbox) for each of them.
    inline void operator()( BBox2i bbox ) const {
      typename BlockThread::Info info( m_func, bbox, m_block_size );

      // Avoid threads altogether in the single-threaded case.
      // Annoyingly, this still creates an unnecessary Mutex.
      if( m_num_threads == 1 ) {
        BlockThread bt( info );
        return bt();
      }

      std::vector<boost::shared_ptr<BlockThread> > generators;
      std::vector<boost::shared_ptr<Thread     > > threads;

      for( uint32 i=0; i<m_num_threads; ++i ) {
        boost::shared_ptr<BlockThread> generator( new BlockThread( info ) );
        generators.push_back( generator );
        boost::shared_ptr<Thread> thread( new Thread( generator ) );
        threads.push_back( thread );
      }

      for( uint32 i=0; i<m_num_threads; ++i ) {
        threads[i]->join();
      }
    }

  }; // End class BlockProcessor

//-----------------------------------------------------------------------------------
// Things below here are other helpful block processing utilities.


  /// Get the default block size to use for block image processing.
  template<typename PixelT>
  Vector2i get_default_block_size(int rows, int cols, int planes=1) {

    const int32 default_blocksize = 2*1024*1024; // 2 megabytes
    // XXX Should the default block configuration be different for
    // very wide images?  Either way we will guess wrong some of
    // the time, so advanced users will have to know what they're
    // doing in any case.
    int32 block_rows = default_blocksize / (planes*cols*int32(sizeof(PixelT)));
    if( block_rows < 1 ) block_rows = 1;
    else if( block_rows > rows ) block_rows = rows;
    return Vector2i( cols, block_rows );
  }


  /// These objects rasterize a full block of image data to be stored in the cache.
  /// - Set up with a source image and an ROI.  When generate() is called, and
  ///    ImageView object is created containing that ROI from the source image.
  /// - This is the class that is used by the Cache class to generate nearly all images
  ///   loaded by vision workbench!
  template <class ImageT>
  class BlockGenerator {
    boost::shared_ptr<ImageT> m_child; ///< Source image view.
    BBox2i m_bbox; ///< ROI of the source image.
  public:
    typedef ImageView<typename ImageT::pixel_type> value_type; ///< A plain image view

    /// Set this object up to represent a region of an image view.
    BlockGenerator( boost::shared_ptr<ImageT> const& child, BBox2i const& bbox )
      : m_child( child ), m_bbox( bbox ) {}

    /// Return the size in bytes that the rasterized object occupies.
    size_t size() const {
      return m_bbox.width() * m_bbox.height() * m_child->planes() * sizeof(typename ImageT::pixel_type);
    }

    /// Rasterize this object into memory from whatever its source is.
    boost::shared_ptr<value_type > generate() const {
      boost::shared_ptr<value_type > ptr( new value_type( m_bbox.width(), m_bbox.height(), m_child->planes() ) );
      m_child->rasterize( *ptr, m_bbox );
      return ptr;
    }
  }; // End class BlockGenerator

  /// Manages a table of BlockGenerator objects spanning an entire image.
  template <class ImageT>
  class BlockGeneratorManager {

    Cache   *m_cache_ptr;
    Vector2i m_block_size;
    int      m_table_width, m_table_height;
    size_t   m_block_table_size;
    boost::shared_array<Cache::Handle<BlockGenerator<ImageT> > > m_block_table;

  public:

    BlockGeneratorManager()
     : m_cache_ptr(0), m_block_size(0,0), m_table_width(0), m_table_height(0), m_block_table_size(0) {}

    /// Fill up m_block_table with a set of BlockGenerator objects.
    void initialize(Cache *cache, Vector2i const block_size, boost::shared_ptr<ImageT> image) {
      m_cache_ptr  = cache;
      m_block_size = block_size;

      // Error checking
      if( m_block_size.x() <= 0 || m_block_size.y() <= 0 ) {
        vw_throw( ArgumentErr() << "BlockGeneratorManager: Illegal block size: " << m_block_size);
      }
      if( !m_cache_ptr ) {
        vw_throw( ArgumentErr() << "BlockGeneratorManager: No cache provided!");
      }

      // Compute the table layout
      m_table_width  = (image->cols()-1) / m_block_size.x() + 1;
      m_table_height = (image->rows()-1) / m_block_size.y() + 1;
      m_block_table_size = m_table_height * m_table_width;
      m_block_table.reset( new Cache::Handle<BlockGenerator<ImageT> >[ m_block_table_size ] );
      BBox2i view_bbox(0,0,image->cols(),image->rows());

      // Iterate through the block positions and insert a generator object for each block
      // into m_block_table.
      for( int32 iy=0; iy<m_table_height; ++iy ) {
        for( int32 ix=0; ix<m_table_width; ++ix ) {
          BBox2i bbox( ix*m_block_size.x(), iy*m_block_size.y(), m_block_size.x(), m_block_size.y() );
          bbox.crop( view_bbox );
          block(ix,iy) = m_cache_ptr->insert( BlockGenerator<ImageT>( image, bbox ) );
        }
      } // End loop through the blocks

    } // End initialize()

    /// Return the block index for a given input pixel.
    Vector2i get_block_index( Vector2i pixel ) const {
      return Vector2i(pixel.x()/m_block_size.x(), pixel.y()/m_block_size.y());
    }

    /// Return the block index for a given input bounding box (which must fall entirely within one block!)
    Vector2i get_block_index( BBox2i bbox ) const {
      Vector2i block_index = get_block_index(bbox.min());
#if VW_DEBUG_LEVEL > 1
      Vector2i block_index_max = get_block_index(bbox.max()-Vector2i(1,1));
      if( block_index != block_index_max ) {
        vw_throw(LogicErr() << "BlockGeneratorManager: bbox " << bbox << " spans more than one cache block!");
      }
#endif
      return block_index;
    }

    /// Return the top left pixel coordinate of a given block.
    Vector2i get_block_start_pixel( Vector2i block_index ) const {
      return Vector2i(block_index.x()*m_block_size.x(),
                      block_index.y()*m_block_size.y());
    }

    /// Throw if the block index is out of bounds.
    void check_block_index( Vector2i block_index ) const {
      int ix = block_index.x();
      int iy = block_index.y();
      if( ix<0 || ix>=m_table_width || iy<0 || iy>=m_table_height )
        vw_throw( ArgumentErr() << "BlockGeneratorManager: Block indices out of bounds, (" << ix
                  << "," << iy << ") of (" << m_table_width << "," << m_table_height << ")" );
    }

    /// Fetch the block generator for the requested block (const ref)
    const Cache::Handle<BlockGenerator<ImageT> >& block( Vector2i block_index ) const {
      int ix = block_index.x();
      int iy = block_index.y();
      check_block_index(block_index);
      return m_block_table[block_index.x() + block_index.y()*m_table_width];
    }

    /// Fetch the block generator for the requested block
    Cache::Handle<BlockGenerator<ImageT> >& block( Vector2i block_index )  {
      int ix = block_index.x();
      int iy = block_index.y();
      check_block_index(block_index);
      return m_block_table[block_index.x() + block_index.y()*m_table_width];
    }

    /// Overload
    Cache::Handle<BlockGenerator<ImageT> >& block( int ix, int iy )  {
      check_block_index(Vector2i(ix, iy));
      return m_block_table[ix + iy*m_table_width];
    }

    /// Return true if there is only a single block
    bool only_one_block() const { return (m_block_table_size==1); }

    /// Shortcut function if you know there is only one block.
    const Cache::Handle<BlockGenerator<ImageT> >& quick_single_block()  const {
      return m_block_table[0];
    }
    Cache::Handle<BlockGenerator<ImageT> >& quick_single_block()  {
      return m_block_table[0];
    }

  }; // End class BlockGeneratorManager


} // namespace image_block
} // namespace vw

#endif // __VW_IMAGE_BLOCKPROCESSOR_H__
