// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// ImageResource to read/write to std streams.

#ifndef __VW_IMAGE_IMAGERESOURCESTREAM_H__
#define __VW_IMAGE_IMAGERESOURCESTREAM_H__

#include <vw/Image/ImageResource.h>
#include <boost/noncopyable.hpp>

namespace vw {

class SrcImageResourceStream : public SrcImageResource, private boost::noncopyable
{
  public:
    typedef std::istream stream_type;

    // Caller is responsible for object lifetime.
    SrcImageResourceStream(stream_type* stream);
    SrcImageResourceStream(stream_type* stream, ImageFormat fmt);
    // Caller gives control of lifetime to this class
    SrcImageResourceStream(boost::shared_ptr<stream_type> stream);
    SrcImageResourceStream(boost::shared_ptr<stream_type> stream, ImageFormat fmt);

    // Will not throw
    virtual const ImageFormat& read_format() const;
    virtual void set_read_format(const ImageFormat& fmt);

    // If m_fmt.complete() is false, these will throw.
    virtual int32 cols() const;
    virtual int32 rows() const;
    virtual int32 planes() const;
    virtual PixelFormatEnum pixel_format() const;
    virtual ChannelTypeEnum channel_type() const;

    // Read the stream and write to the given buffer.
    // If m_fmt is complete, data will be vw::converted to dest. Otherwise, raw
    // bytes will be copied.
    virtual void read( ImageBuffer const& dest, BBox2i const& bbox ) const;

    // Will reset read pointer to start to allow another read.
    // Will throw for unseekable streams.
    virtual void reset();

    virtual bool has_block_read()  const {return false;}
    virtual bool has_nodata_read() const {return false;}

  protected:
    boost::shared_ptr<stream_type> m_stream;
    ImageFormat m_fmt;
};

class DstImageResourceStream : public DstImageResource, private boost::noncopyable
{
  public:
    typedef std::ostream stream_type;

    // Caller is responsible for object lifetime.
    DstImageResourceStream(stream_type* stream);
    DstImageResourceStream(stream_type* stream, ImageFormat fmt);
    // Caller gives control of lifetime to this class
    DstImageResourceStream(boost::shared_ptr<stream_type> stream);
    DstImageResourceStream(boost::shared_ptr<stream_type> stream, ImageFormat fmt);

    virtual const ImageFormat& write_format() const;
    virtual void set_write_format(const ImageFormat& fmt);

    // Read the given buffer and write to the stream.
    // If m_fmt is complete, data will be vw::converted to dest. Otherwise, raw
    // bytes will be copied.
    virtual void write( ImageBuffer const& buf, BBox2i const& bbox );

    virtual void flush();

    // Will reset write pointer to start to allow an overwrite.
    // Will throw for unseekable streams.
    virtual void reset();

    virtual bool has_block_write()  const {return false;}
    virtual bool has_nodata_write() const {return false;}

  protected:
    boost::shared_ptr<stream_type> m_stream;
    ImageFormat m_fmt;
};

class ImageResourceStream : public ImageResource, public SrcImageResourceStream, public DstImageResourceStream
{
  public:
    typedef std::iostream stream_type;

    // Caller is responsible for object lifetime.
    ImageResourceStream(stream_type* stream)
      : SrcImageResourceStream(stream), DstImageResourceStream(stream) {};
    ImageResourceStream(stream_type* stream, ImageFormat fmt)
      : SrcImageResourceStream(stream, fmt), DstImageResourceStream(stream, fmt) {};
    // Caller gives control of lifetime to this class
    ImageResourceStream(boost::shared_ptr<stream_type> stream)
      : SrcImageResourceStream(stream), DstImageResourceStream(stream) {};
    ImageResourceStream(boost::shared_ptr<stream_type> stream, ImageFormat fmt)
      : SrcImageResourceStream(stream, fmt), DstImageResourceStream(stream, fmt) {};

    // If m_fmt.complete() is false, these will throw.
    virtual int32 cols() const                   {return SrcImageResourceStream::cols();}
    virtual int32 rows() const                   {return SrcImageResourceStream::rows();}
    virtual int32 planes() const                 {return SrcImageResourceStream::planes();}
    virtual PixelFormatEnum pixel_format() const {return SrcImageResourceStream::pixel_format();}
    virtual ChannelTypeEnum channel_type() const {return SrcImageResourceStream::channel_type();}

    // Read the stream and write to the given buffer.
    // If m_fmt is complete, data will be vw::converted to dest. Otherwise, raw
    // bytes will be copied.
    virtual void read( ImageBuffer const& dest, BBox2i const& bbox ) const {
      SrcImageResourceStream::read(dest, bbox);
    }

    // Read the given buffer and write to the stream.
    // If m_fmt is complete, data will be vw::converted to dest. Otherwise, raw
    // bytes will be copied.
    virtual void write( ImageBuffer const& buf, BBox2i const& bbox ) {
      DstImageResourceStream::write(buf, bbox);
    }

    // Will throw for unseekable streams.
    virtual void reset() {
      SrcImageResourceStream::reset();
      DstImageResourceStream::reset();
    }

    virtual void flush() {
      DstImageResourceStream::flush();
    }

    virtual bool has_block_write() const  {return DstImageResourceStream::has_block_write();}
    virtual bool has_nodata_write() const {return DstImageResourceStream::has_nodata_write();}
    virtual bool has_block_read() const   {return SrcImageResourceStream::has_block_read();}
    virtual bool has_nodata_read() const  {return SrcImageResourceStream::has_nodata_read();}
};

}

#endif
