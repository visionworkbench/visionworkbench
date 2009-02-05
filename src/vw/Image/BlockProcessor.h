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
    int32 m_block_cols, m_block_rows;
    int m_num_threads;
  public:
    BlockProcessor( FuncT const& func, int32 block_cols, int32 block_rows, int threads = 0 )
      : m_func(func), m_block_cols(block_cols), m_block_rows(block_rows),
        m_num_threads(threads?threads:(Settings::system_settings().default_num_threads())) {}

    class BlockThread {
    public:
      struct Info {
        FuncT const& m_func;
        BBox2i total_bbox, block_bbox;
        int32 block_cols, block_rows;
        Mutex mutex;
        Info( FuncT const& func, BBox2i total_bbox, int32 block_cols, int32 block_rows )
          : m_func(func), total_bbox(total_bbox),
            block_bbox(0,0,(std::min)(block_cols,total_bbox.width()),(std::min)(block_rows,total_bbox.height())),
            block_cols(block_cols), block_rows(block_rows) {}
      };
      
      BlockThread( Info &info ) : info(info) {}

      void operator()() {
        while( true ) {
          BBox2i bbox;
          {
            Mutex::Lock lock(info.mutex);
            if( info.block_bbox.min().y() >= info.total_bbox.height() ) return;
            bbox = info.block_bbox;
            info.block_bbox.min().x() += info.block_cols;
            if( info.block_bbox.min().x() >= info.total_bbox.width() ) {
              info.block_bbox.min().x() = 0;
              info.block_bbox.min().y() += info.block_rows;
              info.block_bbox.max().y() = (std::min)( info.block_rows + info.block_bbox.min().y(), info.total_bbox.height() );
            }
            info.block_bbox.max().x() = (std::min)( info.block_cols + info.block_bbox.min().x(), info.total_bbox.width() );
          }
          info.m_func( bbox );
        }
      }

    private:
      Info &info;
    };

    inline void operator()( BBox2i bbox ) const {
      typename BlockThread::Info info( m_func, bbox, m_block_cols, m_block_rows );

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
