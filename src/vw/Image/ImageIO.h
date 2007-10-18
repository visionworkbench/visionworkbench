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

  /// Write an image view to a resource.
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
