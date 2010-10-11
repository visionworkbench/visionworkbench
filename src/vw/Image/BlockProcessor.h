// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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

  template <class FuncT>
  class BlockProcessor {
    FuncT m_func;
    Vector2i m_block_size;
    int m_num_threads;
  public:
    BlockProcessor( FuncT const& func, Vector2i const& block_size, int threads = 0 )
      : m_func(func), m_block_size(block_size),
        m_num_threads(threads?threads:(vw_settings().default_num_threads())) {}

    // We will construct and call one BlockThread per thread.
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
        BBox2i m_total_bbox, m_block_bbox;
        Vector2i m_block_size;
        Mutex m_mutex;
      };

      BlockThread( Info &info ) : info(info) {}

      void operator()() {
        while( true ) {
          BBox2i bbox;
          {
            // Grab the next bbox to process, and update it with the
            // subsequent bbox for the next thread to grab.
            Mutex::Lock lock(info.mutex());
            if( info.complete() ) return;
            bbox = info.bbox();
            info.advance();
          }
          info.func()( bbox );
        }
      }

    private:
      Info &info;
    };

    inline void operator()( BBox2i bbox ) const {
      typename BlockThread::Info info( m_func, bbox, m_block_size );

      // Avoid threads altogether in the single-threaded case.
      // Annoyingly, this still creates an unnecessary Mutex.
      if( m_num_threads == 1 ) {
        BlockThread bt( info );
        return bt();
      }

      std::vector<boost::shared_ptr<BlockThread> > generators;
      std::vector<boost::shared_ptr<Thread> > threads;

      for( int i=0; i<m_num_threads; ++i ) {
        boost::shared_ptr<BlockThread> generator( new BlockThread( info ) );
        generators.push_back( generator );
        boost::shared_ptr<Thread> thread( new Thread( generator ) );
        threads.push_back( thread );
      }

      for( int i=0; i<m_num_threads; ++i ) {
        threads[i]->join();
      }
    }

  };

} // namespace vw

#endif // __VW_IMAGE_BLOCKPROCESSOR_H__
