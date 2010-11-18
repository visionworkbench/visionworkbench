// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file ImageIO.h
///
/// Functions for reading and writing image views to and from
/// image resources.
///
#ifndef __VW_IMAGE_IMAGEIO_H__
#define __VW_IMAGE_IMAGEIO_H__

#include <vw/Core/ProgressCallback.h>
#include <vw/Core/ThreadPool.h>
#include <vw/Image/ImageResource.h>
#include <vw/Image/ImageView.h>

namespace vw {

  // *******************************************************************
  // Image view reading and writing functions.
  // *******************************************************************

  template <class PixelT>
  inline void read_image( ImageView<PixelT>& dst, SrcImageResource const& src, BBox2i const& bbox ) {
    int32 planes = 1;
    if( ! IsCompound<PixelT>::value ) {
      // The image has a fundamental pixel type
      if( src.planes()>1 && src.channels()>1 )
        vw_throw( ArgumentErr() << "Cannot read a multi-plane multi-channel image resource into a single-channel view." );
      planes = (std::max)( src.planes(), src.channels() );
    }
    dst.set_size( bbox.width(), bbox.height(), planes );
    src.read( dst.buffer(), bbox );
  }

  template <class PixelT>
  inline void read_image( ImageView<PixelT>& dst, SrcImageResource const& src ) {
    read_image( dst, src, BBox2i(0,0,src.cols(),src.rows()) );
  }

  template <class PixelT>
  inline void read_image( ImageView<PixelT> const& dst, SrcImageResource const& src, BBox2i const& bbox ) {
    src.read( dst.buffer(), bbox );
  }

  template <class PixelT>
  inline void read_image( ImageView<PixelT> const& dst, SrcImageResource const& src ) {
    read_image( dst, src, BBox2i(0,0,src.cols(),src.rows()) );
  }

  template <class ImageT>
  inline void read_image( ImageViewBase<ImageT> const& dst, SrcImageResource const& src, BBox2i const& bbox ) {
    ImageView<typename ImageT::pixel_type> intermediate;
    read_image( intermediate, src, bbox );
    dst.impl() = intermediate;
  }

  template <class ImageT>
  inline void read_image( ImageViewBase<ImageT> const& dst, SrcImageResource const& src ) {
    read_image( dst, src, BBox2i(0,0,src.cols(),src.rows()) );
  }

  template <class PixelT>
  inline void write_image( DstImageResource &dst, ImageView<PixelT> const& src, BBox2i const& bbox ) {
    dst.write( src.buffer(), bbox );
  }

  template <class PixelT>
  inline void write_image( DstImageResource &dst, ImageView<PixelT> const& src ) {
    write_image( dst, src, BBox2i(0,0,src.cols(),src.rows()) );
  }

  template <class ImageT>
  inline void write_image( DstImageResource &dst, ImageViewBase<ImageT> const& src, BBox2i const& bbox ) {
    ImageView<typename ImageT::pixel_type> intermediate = src;
    write_image( dst, intermediate, bbox );
  }



  // -----------------------------------------------------------------------
  //
  // Doing multi-threaded block rasterization and writing to disk
  // correctly is a little tricky, because blocks almost always must
  // be written _in order_ into the file.  If one of the rasterizing
  // threads for an early block falls behind (which is not all that
  // unlikely), it can cause a backup on the write queue.  Depending
  // on the size of the blocks being rasterized, this can lead to
  // large allocations of memory as rasterized tiles accumulate and
  // sit waiting to be written to disk.
  //
  // To fix this, we need this semaphore to help us meet the following
  // condition:
  //
  // We rasterize _at most_ N blocks at a time, and it will never
  // get more than N blocks ahead of the last block that was written.
  //
  // Of course, one slow rasterization thread can hold up the entire
  // process, but this is the price we pay for guranteed ordering when
  // writing tiles.
  class CountingSemaphore {
    Condition m_block_condition;
    Mutex m_mutex;
    int m_max, m_last_job_index;

  public:
    CountingSemaphore() : m_max( vw_settings().default_num_threads() ),
                          m_last_job_index(-1) {}
    CountingSemaphore( int max ) : m_max(max), m_last_job_index(-1) {}

    // Call to wait for a turn until the number of threads in a area
    // decrements.
    void wait( int job_index ) {
      Mutex::Lock lock(m_mutex);
      while ( job_index > m_last_job_index + m_max ) {
        m_block_condition.wait(lock);
      }
    }

    // Please call when ever a process finishes it's turn
    void notify() {
      {
        Mutex::Lock lock(m_mutex);
        m_last_job_index++;
      }
      m_block_condition.notify_all();
    }
  };

  // This task generator manages the rasterizing and writing of images to disk.
  //
  // Only one thread can be writing to the ImageResource at any given
  // time, however several threads can be rasterizing simultaneously.
  //
  class ThreadedBlockWriter : private boost::noncopyable {

    boost::shared_ptr<FifoWorkQueue> m_rasterize_work_queue;
    boost::shared_ptr<OrderedWorkQueue> m_write_work_queue;
    CountingSemaphore m_write_queue_limit;

    // ----------------------------- TASK TYPES (2) --------------------------

    template <class PixelT>
    class WriteBlockTask : public Task {
      DstImageResource& m_resource;
      ImageView<PixelT> m_image_block;
      BBox2i m_bbox;
      int m_idx;
      CountingSemaphore& m_write_finish_event;

    public:
      WriteBlockTask(DstImageResource& resource, ImageView<PixelT> const& image_block,
                     BBox2i bbox, int idx, CountingSemaphore& write_finish_event) :
      m_resource(resource), m_image_block(image_block), m_bbox(bbox), m_idx(idx), m_write_finish_event(write_finish_event) {}

      virtual ~WriteBlockTask() {}
      virtual void operator() () {
        vw_out(DebugMessage, "image") << "Writing block " << m_idx << " at " << m_bbox << "\n";
        m_resource.write( m_image_block.buffer(), m_bbox );
        m_write_finish_event.notify();
      }
    };

    // -----------------------------

    template <class ViewT>
    class RasterizeBlockTask : public Task {
      ThreadedBlockWriter &m_parent;
      DstImageResource& m_resource;
      ViewT const& m_image;
      BBox2i m_bbox;
      int m_index;
      int m_total_num_blocks;
      SubProgressCallback m_progress_callback;
      CountingSemaphore& m_write_finish_event;

    public:
      RasterizeBlockTask(ThreadedBlockWriter &parent, DstImageResource& resource,
                         ImageViewBase<ViewT> const& image, BBox2i const& bbox,
                         int index, int total_num_blocks,
                         CountingSemaphore& write_finish_event,
                         const ProgressCallback &progress_callback = ProgressCallback::dummy_instance()) :
      m_parent(parent), m_resource(resource), m_image(image.impl()), m_bbox(bbox), m_index(index),
        m_progress_callback(progress_callback,0.0,1.0/float(total_num_blocks)), m_write_finish_event(write_finish_event) {}

      virtual ~RasterizeBlockTask() {}
      virtual void operator()() {

        m_write_finish_event.wait(m_index);

        vw_out(DebugMessage, "image") << "Rasterizing block " << m_index << " at " << m_bbox << "\n";
        // Rasterize the block
        ImageView<typename ViewT::pixel_type> image_block( crop(m_image, m_bbox) );

        // Report progress
        m_progress_callback.report_incremental_progress(1.0);

        // With rasterization complete, we queue up a request to write this block to disk.
        boost::shared_ptr<Task> write_task ( new WriteBlockTask<typename ViewT::pixel_type>( m_resource, image_block, m_bbox, m_index, m_write_finish_event ) );

        m_parent.add_write_task(write_task, m_index);
      }
    };

    // -----------------------------

    void add_write_task(boost::shared_ptr<Task> task, int index) { m_write_work_queue->add_task(task, index); }
    void add_rasterize_task(boost::shared_ptr<Task> task) { m_rasterize_work_queue->add_task(task); }

  public:
    ThreadedBlockWriter() : m_write_queue_limit(vw_settings().write_pool_size()) {
      m_rasterize_work_queue = boost::shared_ptr<FifoWorkQueue>( new FifoWorkQueue() );
      m_write_work_queue = boost::shared_ptr<OrderedWorkQueue>( new OrderedWorkQueue(1) );
    }

    // Add a block to be rasterized.  You can optionally supply an
    // index, which will indicate the order in which this block should
    // be written to disk.
    template <class ViewT>
    void add_block(DstImageResource& resource, ImageViewBase<ViewT> const& image, BBox2i const& bbox, int index, int total_num_blocks,
                   const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() ) {
      boost::shared_ptr<Task> task( new RasterizeBlockTask<ViewT>(*this, resource, image, bbox, index, total_num_blocks, m_write_queue_limit, progress_callback) );
      this->add_rasterize_task(task);
    }

    void process_blocks() {
      m_rasterize_work_queue->join_all();
      m_write_work_queue->join_all();
    }
  };


  /// Write an image view to a resource.
  template <class ImageT>
  void block_write_image( DstImageResource& resource, ImageViewBase<ImageT> const& image,
                          const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() ) {

    VW_ASSERT( image.impl().cols() != 0 && image.impl().rows() != 0 && image.impl().planes() != 0,
               ArgumentErr() << "write_image: cannot write an empty image to a resource" );

    // Set the progress meter to zero.
    progress_callback.report_progress(0);
    if (progress_callback.abort_requested())
      vw_throw( Aborted() << "Aborted by ProgressCallback" );

    // Set up the threaded block writer object, which will manage
    // rasterizing and writing images to disk one block (and one
    // thread) at a time.
    ThreadedBlockWriter block_writer;

    const size_t rows = image.impl().rows();
    const size_t cols = image.impl().cols();

    // Write the image to disk in blocks.  We may need to revisit
    // the order in which these blocks are rasterized, but for now
    // it rasterizes blocks from left to right, then top to bottom.
    Vector2i block_size(cols, rows);
    if (resource.has_block_write())
      block_size = resource.block_write_size();

    int total_num_blocks = ((rows-1)/block_size.y()+1) * ((cols-1)/block_size.x()+1);
    vw_out(DebugMessage,"image") << "ThreadedBlockWriter: writing " << total_num_blocks << " blocks.\n";

    for (size_t j = 0; j < rows; j+= block_size.y()) {
      for (size_t i = 0; i < cols; i+= block_size.x()) {

        // Rasterize and save this image block
        BBox2i current_bbox(Vector2i(i,j),
                            Vector2i(std::min<size_t>(i+block_size.x(),cols),
                                     std::min<size_t>(j+block_size.y(),rows)));

        // Add a task to rasterize this image block.  A seperate task
        // to write the results to disk is generated automatically
        // when rasterization is complete.
        int col_blocks = int( ceil(float(cols)/float(block_size.x())) );
        int i_block_index = int(i/block_size.x());
        int j_block_index = int(j/block_size.y());
        int index = j_block_index*col_blocks+i_block_index;

        vw_out(VerboseDebugMessage,"image") << "ThreadedBlockWriter: Adding block " << index+1 << "/"<< total_num_blocks << " : " << current_bbox << "\n";
        block_writer.add_block(resource, image, current_bbox, index, total_num_blocks, progress_callback );
      }
    }

    // Start the threaded block writer and wait for all tasks to finish.
    block_writer.process_blocks();
    progress_callback.report_finished();
  }


  template <class ImageT>
  void write_image( DstImageResource& resource, ImageViewBase<ImageT> const& image,
                    const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() ) {

    VW_ASSERT( image.impl().cols() != 0 && image.impl().rows() != 0 && image.impl().planes() != 0,
               ArgumentErr() << "write_image: cannot write an empty image to a resource" );

    // Initialize the progress callback
    progress_callback.report_progress(0);

    const int32 rows = boost::numeric_cast<int32>(image.impl().rows());
    const int32 cols = boost::numeric_cast<int32>(image.impl().cols());

    // Write the image to disk in blocks.  We may need to revisit
    // the order in which these blocks are rasterized, but for now
    // it rasterizes blocks from left to right, then top to bottom.
    Vector2i block_size(cols, rows);
    if (resource.has_block_write())
      block_size = resource.block_write_size();

    size_t total_num_blocks = ((rows-1)/block_size.y()+1) * ((cols-1)/block_size.x()+1);
    for (int32 j = 0; j < rows; j+= block_size.y()) {
      for (int32 i = 0; i < cols; i+= block_size.x()) {

        vw_out(DebugMessage, "fileio") << "ImageIO writing block at [" << i << " " << j << "]/["
                                       << rows << " " << cols
                                       << "]    size = " << block_size.x() << " x " <<  block_size.y() << "\n";

        // Update the progress callback.
        if (progress_callback.abort_requested())
          vw_throw( Aborted() << "Aborted by ProgressCallback" );

        float processed_row_blocks = float(j/block_size.y()*((cols-1)/block_size.x()+1));
        float processed_col_blocks = float(i/block_size.x());
        progress_callback.report_progress((processed_row_blocks + processed_col_blocks) / static_cast<float>(total_num_blocks));

        // Rasterize and save this image block
        BBox2i current_bbox(Vector2i(i,j),
                            Vector2i(std::min<int32>(i+block_size.x(),cols),
                                     std::min<int32>(j+block_size.y(),rows)));

        // Rasterize the current image block into a region of memory
        // and send it off to the resource.
        ImageView<typename ImageT::pixel_type> image_block( crop(image.impl(), current_bbox) );
        ImageBuffer buf = image_block.buffer();
        resource.write( buf, current_bbox );

      }
    }
    progress_callback.report_finished();
  }

} // namespace vw

#endif // __VW_IMAGE_IMAGEIO_H__
