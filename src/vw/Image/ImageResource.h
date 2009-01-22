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

/// \file ImageResource.h
/// 
/// Defines the abstract base image resource type.
///
#ifndef __VW_IMAGE_IMAGERESOURCE_H__
#define __VW_IMAGE_IMAGERESOURCE_H__

#include <vw/Math/Vector.h>
#include <vw/Math/BBox.h>

#include <vw/Image/PixelTypeInfo.h>

namespace vw {

  // Forward declaration
  struct ImageBuffer;


  /// Copies image pixel data from the source buffer to the destination 
  /// buffer, converting the pixel format and channel type as required.
  void convert( ImageBuffer const& dst, ImageBuffer const& src, bool rescale=false );


  /// Describes the format of an image, i.e. its dimensions, pixel
  /// structure, and channel type.
  struct ImageFormat {
    int32 cols, rows, planes;
    PixelFormatEnum pixel_format;
    ChannelTypeEnum channel_type;
    
    ImageFormat()
      : cols(0), rows(0), planes(0),
        pixel_format(VW_PIXEL_UNKNOWN),
        channel_type(VW_CHANNEL_UNKNOWN)
    {}

  };


  /// Base class from which specific image resources derive.
  class ImageResource { 
  public:

    virtual ~ImageResource() {};

    /// Returns the number of columns in an image resource.
    virtual int32 cols() const = 0;

    /// Returns the number of rows in an image resource.
    virtual int32 rows() const = 0;

    /// Returns the number of planes in an image resource.
    virtual int32 planes() const = 0;

    /// Returns the number of channels in a image resource.
    int32 channels() const { return num_channels( pixel_format() ); }

    /// Returns the native pixel format of the resource.
    virtual PixelFormatEnum pixel_format() const = 0;

    /// Returns the native channel type of the resource.
    virtual ChannelTypeEnum channel_type() const = 0;

    /// Read the image resource at the given location into the given buffer.
    virtual void read( ImageBuffer const& buf, BBox2i const& bbox ) const = 0;

    /// Write the given buffer to the image resource at the given location.
    virtual void write( ImageBuffer const& buf, BBox2i const& bbox ) = 0;

    /// Returns the optimal block size/alignment for partial reads or writes.
    virtual Vector2i native_block_size() const { return Vector2i(cols(),rows()); }

    /// Force any changes to be written to the resource.
    virtual void flush() {}

  };


  /// Represents a generic image buffer in memory, with dimensions and
  /// pixel format specified at run time.  This class does not
  /// allocate any memory, but rather provides a common format for
  /// describing an existing in-memory buffer of pixels.  The primary
  /// purpose of this class is to provide some common ground for
  /// converting between image formats using the convert() function.
  /// To allocate a fresh buffer for an image, see ImageView.
  struct ImageBuffer {
    void* data;
    ImageFormat format;
    ptrdiff_t cstride, rstride, pstride;
    bool unpremultiplied;
    
    /// Default constructor; constructs an undefined buffer
    ImageBuffer()
      : data(0), format(),
        cstride(0), rstride(0), pstride(0),
        unpremultiplied(false)
    {}

    virtual ~ImageBuffer() {}

    /// Returns the number of columns in the bufffer.
    inline int32 cols() const { return format.cols; }

    /// Returns the number of rows in the bufffer.
    inline int32 rows() const { return format.rows; }

    /// Returns the number of planes in the bufffer.
    inline int32 planes() const { return format.planes; }

    /// Returns the native pixel format of the bufffer.
    inline PixelFormatEnum pixel_format() const { return format.pixel_format; }

    /// Returns the native channel type of the bufffer.
    inline ChannelTypeEnum channel_type() const { return format.channel_type; }

    /// Returns a cropped version of this bufffer.
    inline ImageBuffer cropped( BBox2i const& bbox ) const {
      ImageBuffer self = *this;
      self.data = (uint8*)self.data + cstride*bbox.min().x() + rstride*bbox.min().y();
      self.format.cols = bbox.width();
      self.format.rows = bbox.height();
      return self;
    }

    /// Read the image resource at the given location into the given buffer.
    inline void read( ImageBuffer const& buf, BBox2i const& bbox ) const {
      convert( buf, cropped(bbox) );
    }

    /// Write the given buffer to the image resource at the given location.
    inline void write( ImageBuffer const& buf, BBox2i const& bbox ) {
      convert( cropped(bbox), buf );
    }

    /// Return a pointer to the pixel at (u,v,p)
    inline void* operator()( int32 i, int32 j, int32 p = 0 ) const {
      return ((uint8*)data) + (i*cstride + j*rstride + p*pstride);
    }
    
  };

} // namespace vw

#endif // __VW_IMAGE_IMAGERESOURCE_H__
