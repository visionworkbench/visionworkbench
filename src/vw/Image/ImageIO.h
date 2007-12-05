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

/// \file ImageView.h
/// 
/// Functions for reading and writing image views to and from 
/// image resources.
///
#ifndef __VW_IMAGE_IMAGEIO_H__
#define __VW_IMAGE_IMAGEIO_H__

#include <vw/Core/Exception.h>
#include <vw/Core/ProgressCallback.h>
#include <vw/Core/ThreadPool.h>

#include <vw/Image/ImageResource.h>
#include <vw/Image/ImageView.h>

namespace vw {

  // *******************************************************************
  // Image view reading and writing functions.
  // *******************************************************************

  template <class PixelT>
  inline void read_image( ImageView<PixelT>& dst, ImageResource const& src, BBox2i const& bbox ) {
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
  inline void read_image( ImageView<PixelT>& dst, ImageResource const& src ) {
    read_image( dst, src, BBox2i(0,0,src.cols(),src.rows()) );
  }

  template <class PixelT>
  inline void read_image( ImageView<PixelT> const& dst, ImageResource const& src, BBox2i const& bbox ) {
    src.read( dst.buffer(), bbox );
  }

  template <class PixelT>
  inline void read_image( ImageView<PixelT> const& dst, ImageResource const& src ) {
    read_image( dst, src, BBox2i(0,0,src.cols(),src.rows()) );
  }

  template <class ImageT>
  inline void read_image( ImageViewBase<ImageT> const& dst, ImageResource const& src, BBox2i const& bbox ) {
    ImageView<typename ImageT::pixel_type> intermediate;
    read_image( intermediate, src, bbox );
    dst.impl() = intermediate;
  }

  template <class ImageT>
  inline void read_image( ImageViewBase<ImageT> const& dst, ImageResource const& src ) {
    read_image( dst, src, BBox2i(0,0,src.cols(),src.rows()) );
  }

  template <class PixelT>
  inline void write_image( ImageResource &dst, ImageView<PixelT> const& src, BBox2i const& bbox ) {
    dst.write( src.buffer(), bbox );
  }

  template <class PixelT>
  inline void write_image( ImageResource &dst, ImageView<PixelT> const& src ) {
    write_image( dst, src, BBox2i(0,0,dst.cols(),dst.rows()) );
  }

  template <class ImageT>
  inline void write_image( ImageResource &dst, ImageViewBase<ImageT> const& src, BBox2i const& bbox ) {
    ImageView<typename ImageT::pixel_type> intermediate = src;
    write_image( dst, intermediate, bbox );
  }

  // This task generator manages the rasterizing and writing of images to disk.
  //
  // Only one thread can be writing to the ImageResource at any given
  // time, however several threads can be rasterizing simultaneously.
  // 
  class ThreadedBlockWriter {

    boost::shared_ptr<FifoWorkQueue> m_rasterize_work_queue;
    boost::shared_ptr<OrderedWorkQueue> m_write_work_queue;

    // ----------------------------- TASK TYPES (2) --------------------------------

    template <class PixelT>
    class WriteBlockTask : public Task {
      ImageResource& m_resource;
      ImageView<PixelT> m_image_block;
      BBox2i m_bbox; 

    public:
      WriteBlockTask(ImageResource& resource, ImageView<PixelT> const& image_block, BBox2i bbox) : 
        m_resource(resource), m_image_block(image_block), m_bbox(bbox) {}

      virtual ~WriteBlockTask() {}
      virtual void operator() () { m_resource.write( m_image_block.buffer(), m_bbox ); }
    };

    // ------

    template <class ViewT>
    class RasterizeBlockTask : public Task {
      ThreadedBlockWriter &m_parent;
      ImageResource& m_resource;
      ViewT const& m_image;
      BBox2i m_bbox;
      int m_index;
      
    public:
      RasterizeBlockTask(ThreadedBlockWriter &parent, ImageResource& resource, ImageViewBase<ViewT> const& image, BBox2i const& bbox, int index) :
        m_parent(parent), m_resource(resource), m_image(image.impl()), m_bbox(bbox), m_index(index) {}

      virtual ~RasterizeBlockTask() {}
      virtual void operator()() { 
        ImageView<typename ViewT::pixel_type> image_block( crop(m_image, m_bbox) );

        // Once rasterization is complete, we queue up a request to write this block to disk.
        boost::shared_ptr<Task> write_task ( new WriteBlockTask<typename ViewT::pixel_type>( m_resource, image_block, m_bbox) );
        m_parent.add_write_task(write_task, m_index);
      }
    };

    // ----------------------------- ------------------------ --------------------------------

    void add_write_task(boost::shared_ptr<Task> task, int index) { m_write_work_queue->add_task(task, index); }
    void add_rasterize_task(boost::shared_ptr<Task> task) { m_rasterize_work_queue->add_task(task); }
    
  public:
    ThreadedBlockWriter() {
      m_rasterize_work_queue = boost::shared_ptr<FifoWorkQueue>( new FifoWorkQueue() );
      m_write_work_queue = boost::shared_ptr<OrderedWorkQueue>( new OrderedWorkQueue(1) );
    }


    // Add a block to be rasterized.  You can optionally supply an
    // index, which will indicate the order in which this block should
    // be written to disk.
    template <class ViewT>
    void add_block(ImageResource& resource, ImageViewBase<ViewT> const& image, BBox2i const& bbox, int index = -1) { 
      boost::shared_ptr<Task> task( new RasterizeBlockTask<ViewT>(*this, resource, image, bbox, index) );
      this->add_rasterize_task(task);
    }
                        
    void process_blocks( const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() ) {
      
      // I need to think about how to properly measure progress in the
      // multithreaded system.
      progress_callback.report_progress(0);
      if (progress_callback.abort_requested()) 
        vw_throw( Aborted() << "Aborted by ProgressCallback" );
      //       progress_callback.report_progress((processed_row_blocks + processed_col_blocks)/total_num_blocks);
      progress_callback.report_finished();
      
      m_rasterize_work_queue->join_all();
      m_write_work_queue->join_all();
    }
  };


  /// Write an image view to a resource.
  template <class ImageT>
  void block_write_image( ImageResource& resource, ImageViewBase<ImageT> const& image, 
                          const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() ) {

    VW_ASSERT( image.impl().cols() != 0 && image.impl().rows() != 0 && image.impl().planes() != 0,
               ArgumentErr() << "write_image: cannot write an empty image to a resource" );

    // Set up the threaded block writer object, which will manage
    // rasterizing and writing images to disk one block (and one
    // thread) at a time.
    ThreadedBlockWriter block_writer;
    
    // Write the image to disk in blocks.  We may need to revisit
    // the order in which these blocks are rasterized, but for now
    // it rasterizes blocks from left to right, then top to bottom.
    Vector2i block_size = resource.native_block_size();
    for (int32 j = 0; j < (int32)resource.rows(); j+= block_size.y()) {
      for (int32 i = 0; i < (int32)resource.cols(); i+= block_size.x()) {
        
        // Rasterize and save this image block
        BBox2i current_bbox(Vector2i(i,j),
                            Vector2i(std::min(i+block_size.x(),(int32)(resource.cols())),
                                     std::min(j+block_size.y(),(int32)(resource.rows()))));
        
        // Add a task to rasterize this image block.  A seperate task
        // to write the results to disk is generated automatically
        // when rasterization is complete.
        int col_blocks = int( ceil(float(resource.cols())/float(block_size.x())) );
        int i_block_index = int(i/block_size.x());
        int j_block_index = int(j/block_size.y());
        int index = j_block_index*col_blocks+i_block_index;

        // For debugging:
        //    int total_num_blocks = ((resource.rows()-1)/block_size.y()+1) * ((resource.cols()-1)/block_size.x()+1);
        //   vw_out(0) << "Adding block " << index+1 << "/"<< total_num_blocks << " : " << current_bbox << "\n";

        block_writer.add_block(resource, image, current_bbox, index);
      }
    }
    
    // Start the threaded block writer and wait for all tasks to finish.
    block_writer.process_blocks(progress_callback);
  }
  

  template <class ImageT>
  void write_image( ImageResource& resource, ImageViewBase<ImageT> const& image, 
                    const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() ) {
    
    VW_ASSERT( image.impl().cols() != 0 && image.impl().rows() != 0 && image.impl().planes() != 0,
               ArgumentErr() << "write_image: cannot write an empty image to a resource" );

    // Initialize the progress callback
    progress_callback.report_progress(0);

    // Write the image to disk in blocks.  We may need to revisit
    // the order in which these blocks are rasterized, but for now
    // it rasterizes blocks from left to right, then top to bottom.
    Vector2i block_size = resource.native_block_size();
    int total_num_blocks = ((resource.rows()-1)/block_size.y()+1) * ((resource.cols()-1)/block_size.x()+1);
    for (int32 j = 0; j < (int32)resource.rows(); j+= block_size.y()) {
      for (int32 i = 0; i < (int32)resource.cols(); i+= block_size.x()) {
        
        // Update the progress callback.
        if (progress_callback.abort_requested()) 
          vw_throw( Aborted() << "Aborted by ProgressCallback" );

        float processed_row_blocks = float(j/block_size.y()*((resource.cols()-1)/block_size.x()+1));
        float processed_col_blocks = float(i/block_size.x());
        progress_callback.report_progress((processed_row_blocks + processed_col_blocks)/total_num_blocks);

        // Rasterize and save this image block
        BBox2i current_bbox(Vector2i(i,j),
                            Vector2i(std::min(i+block_size.x(),(int32)(resource.cols())),
                                     std::min(j+block_size.y(),(int32)(resource.rows()))));

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
