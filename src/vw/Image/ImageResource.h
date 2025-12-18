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


/// \file ImageResource.h
///
/// Defines the abstract base image resource types.  These types are for
///  abstracting low level image writing to the level of writing entire
///  ROIs of an image to a buffer in memory or on disk.  The
///  FileIO/DiskImageView.h child of this class is the parent for all
///  disk-based implementations.
///
#ifndef __VW_IMAGE_IMAGERESOURCE_H__
#define __VW_IMAGE_IMAGERESOURCE_H__

#include <vw/Math/Vector.h>
#include <vw/Math/BBox.h>

#include <vw/Image/PixelTypeInfo.h>

namespace vw {

  // Forward declaration
  struct ImageBuffer;

  /// Describes the format of an image, i.e. its dimensions, pixel
  /// structure, and channel type.
  struct ImageFormat {
    uint32          cols, 
                    rows, 
                    planes;
    PixelFormatEnum pixel_format;
    ChannelTypeEnum channel_type;
    bool            premultiplied;

    ImageFormat()
      : cols(0), rows(0), planes(0),
        pixel_format(VW_PIXEL_UNKNOWN),
        channel_type(VW_CHANNEL_UNKNOWN),
        premultiplied(true)
    {}

    // Does this represent a fully-specified data format?
    bool complete() const {
      return   cols != 0
          &&   rows != 0
          && planes != 0
          && num_channels_nothrow(pixel_format) > 0
          && channel_size_nothrow(channel_type) > 0;
    }

    inline bool simple_convert(const ImageFormat& b) const {
      return simple_conversion(channel_type, b.channel_type)
          && simple_conversion(pixel_format, b.pixel_format)
          && premultiplied == b.premultiplied;
    }

    inline bool same_size(const ImageFormat& b) const {
      return cols == b.cols && rows == b.rows && planes == b.planes;
    }

    // These are only valid once you've populated this. No checking is performed.
    size_t cstride  () const {return channel_size(channel_type) * num_channels(pixel_format);}
    size_t rstride  () const {return cstride() * cols;  }
    size_t pstride  () const {return rstride() * rows;  }
    size_t byte_size() const {return pstride() * planes;}
  }; // End class ImageFormat

  /// Dumps an ImageFormat to a std::ostream
  inline std::ostream& operator<<( std::ostream& os, ImageFormat const& f ) {

    os << "ImageFormat rows: " << f.rows << "  cols: " << f.cols
       << "  planes: " << f.planes << "  pixel_format: " << f.pixel_format
       << "  channel_type: " << f.channel_type
       << "  premultiplied: " << f.premultiplied;
    return os;
  }

  /// Copies image pixel data from the source buffer to the destination
  /// buffer, converting the pixel format and channel type as required.
  void convert( ImageBuffer const& dst, ImageBuffer const& src, bool rescale=false );
  
  /// Throws an exception if src cannot be converted to dst using the convert() function.
  /// - Using this function allows us to throw a legible error message instead of gibberish.
  void can_convert(ImageFormat const& dst, ImageFormat const& src);


  /// A read-only image resource
  class SrcImageResource {
    public:
      virtual ~SrcImageResource() {}

      virtual int32 cols  () const {return format().cols;  }
      virtual int32 rows  () const {return format().rows;  }
      virtual int32 planes() const {return format().planes;}

      /// Returns the number of channels in a image resource.
      int32 channels() const { return num_channels( pixel_format() ); }

      /// Returns the native pixel format of the resource.
      virtual PixelFormatEnum pixel_format() const {return format().pixel_format;}

      /// Returns the native channel type of the resource.
      virtual ChannelTypeEnum channel_type() const {return format().channel_type;}

      // Returns the image format as a single object
      virtual ImageFormat format() const = 0;

      /// Read the image resource at the given location into the given buffer.
      virtual void read( ImageBuffer const& buf, BBox2i const& bbox ) const = 0;

      /// Does this resource support block reads?
      // If you override this to true, you must implement the other block_read functions
      virtual bool has_block_read() const = 0;

      /// Returns the preferred block size/alignment for partial reads.
      virtual Vector2i block_read_size() const { return Vector2i(cols(),rows()); }

      // Does this resource have a nodata value?
      // If you override this to true, you must implement the other nodata_read functions
      virtual bool has_nodata_read() const = 0;

      /// Fetch this ImageResource's nodata value
      virtual double nodata_read() const {
        vw_throw(NoImplErr() << "This ImageResource does not support nodata_read().");
        return 0.0;
      }

      /// Return a pointer to the data in the same format as format(). This
      /// might cause a copy, depending on implementation. The shared_ptr will
      /// handle cleanup.
      virtual boost::shared_array<const uint8> native_ptr() const;
      virtual size_t native_size() const;
  };

  /// A write-only image resource
  class DstImageResource {
    public:
      virtual ~DstImageResource() {}

      /// Write the given buffer to the image resource at the given location.
      virtual void write( ImageBuffer const& buf, BBox2i const& bbox ) = 0;

      // Does this resource support block writes?
      // If you override this to true, you must implement the other block_write functions
      virtual bool has_block_write() const = 0;

      /// Gets the preferred block size/alignment for partial writes.
      virtual Vector2i block_write_size() const {
        vw_throw(NoImplErr() << "This ImageResource does not support block writes");
        return Vector2i();
      }

      /// Sets the preferred block size/alignment for partial writes.
      virtual void set_block_write_size(const Vector2i& /*v*/) {
        vw_throw(NoImplErr() << "This ImageResource does not support block writes");
      }

      // Does this resource have an output nodata value?
      // If you override this to true, you must implement the other nodata_write functions
      virtual bool has_nodata_write() const = 0;

      /// Set a nodata value that will be stored in the underlying stream
      virtual void set_nodata_write( double /*value*/ ) {
        vw_throw(NoImplErr() << "This ImageResource does not support set_nodata_write().");
      }

      /// Force any changes to be written to the resource.
      virtual void flush() = 0;
  };

  // A read-write image resource
  class ImageResource : public SrcImageResource, public DstImageResource {};

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
    ssize_t cstride, rstride, pstride;

    /// Default constructor; constructs an undefined buffer
    ImageBuffer()
      : data(0), format(),
        cstride(0), rstride(0), pstride(0)
    {}

    /// Populates stride information from format
    explicit ImageBuffer(ImageFormat format, void *data)
      : data(data), format(format),
        cstride(channel_size(format.channel_type) * num_channels(format.pixel_format)),
        rstride(cstride * format.cols), pstride(rstride * format.rows)
    {}

    virtual ~ImageBuffer() {}

    inline int32 cols  () const { return format.cols;   } ///< Returns the number of columns in the buffer.
    inline int32 rows  () const { return format.rows;   } ///< Returns the number of rows    in the buffer.
    inline int32 planes() const { return format.planes; } ///< Returns the number of planes  in the buffer.

    /// Returns the native pixel format of the bufffer.
    inline PixelFormatEnum pixel_format() const { return format.pixel_format; }

    /// Returns the native channel type of the bufffer.
    inline ChannelTypeEnum channel_type() const { return format.channel_type; }

    /// Returns the size (in bytes) of the data described by this buffer
    inline size_t byte_size() const {
      return planes() * pstride;
    }

    /// Returns a cropped version of this bufffer.
    inline ImageBuffer cropped( BBox2i const& bbox ) const {
      ImageBuffer self = *this;
      self.data = (uint8*)self.data + cstride*bbox.min().x() + rstride*bbox.min().y();
      self.format.cols = bbox.width();
      self.format.rows = bbox.height();
      return self;
    }

    /// Read the image resource at the given location into the given buffer.
    /// - Though the ImageBuffer object is const, the contents of the buffer will change!
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
