// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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

    /// Returns the preferred block size/alignment for partial reads or writes.
    virtual Vector2i block_size() const { return Vector2i(cols(),rows()); }

    /// Set the preferred block size/alignment for partial reads or writes.
    virtual void set_block_size( Vector2i const& ) {
      vw_throw(NoImplErr() << "This ImageResource does not support set_block_size().");
    };

    /// Query whether this ImageResource has a nodata value
    virtual bool has_nodata_value() const { return false; }

    /// Fetch this ImageResource's nodata value
    virtual double nodata_value() const { return 0.0; }

    /// Set the preferred block size/alignment for partial reads or writes.
    virtual void set_nodata_value( double /*value*/ ) {
      vw_throw(NoImplErr() << "This ImageResource does not support set_nodata_value().");
    };

    /// Force any changes to be written to the resource.
    virtual void flush() {}

  };


  /// Acts as a proxy for an arbitrary ImageResource subclass, which is
  /// holds by shared pointer.  This makes it easy to operate on arbitrary
  /// image resources without having to worry about memory management and
  /// object lifetime yourself.
  class ImageResourceRef : public ImageResource {
    boost::shared_ptr<ImageResource> m_resource;
  public:
    // This constructor takes ownership over the resource it's given.
    template <class ResourceT>
    ImageResourceRef( ResourceT *resource ) {
      boost::shared_ptr<ResourceT> resource_ptr( resource );
      m_resource = resource_ptr;
    }

    ImageResourceRef( boost::shared_ptr<ImageResource> resource )
      : m_resource( resource ) {}

    ~ImageResourceRef() {}

    int32 cols() const { return m_resource->cols(); }
    int32 rows() const { return m_resource->rows(); }
    int32 planes() const { return m_resource->planes(); }

    PixelFormatEnum pixel_format() const { return m_resource->pixel_format(); }
    ChannelTypeEnum channel_type() const { return m_resource->channel_type(); }

    void read( ImageBuffer const& buf, BBox2i const& bbox ) const { m_resource->read(buf,bbox); }
    void write( ImageBuffer const& buf, BBox2i const& bbox ) { m_resource->write(buf, bbox); }

    Vector2i block_size() const { return m_resource->block_size(); }
    void set_block_size( Vector2i const& size ) { m_resource->set_block_size(size); }

    void flush() { m_resource->flush(); }
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

    /// Populates stride information from format
    explicit ImageBuffer(ImageFormat format, void *data, bool unpremultiplied = false)
      : data(data), format(format),
        cstride(channel_size(format.channel_type) * num_channels(format.pixel_format)),
        rstride(cstride * format.cols), pstride(rstride * format.rows),
        unpremultiplied(unpremultiplied)
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
