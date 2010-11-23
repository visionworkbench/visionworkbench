// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Image/ImageResourceStream.h>
#include <vw/Core/Debugging.h>

using namespace vw;

namespace {
  void noop_istream(SrcImageResourceStream::stream_type*)   {}
  void noop_ostream(DstImageResourceStream::stream_type*)   {}

  void perform_read(std::istream* stream, char* buf, size_t size) {
    stream->read(buf, size);
    if (stream->fail()) {
      stream->clear();
      vw::vw_throw(vw::IOErr() <<  "Failed to read enough data from stream");
    }
  }

  void perform_write(std::ostream* stream, const char* buf, size_t size) {
    stream->write(buf, size);
    if (stream->bad()) {
      stream->clear();
      vw::vw_throw(vw::IOErr() <<  "Failed to write all data to stream");
    }
  }
}

SrcImageResourceStream::SrcImageResourceStream(stream_type* stream)
  : m_stream(stream, noop_istream) {}

SrcImageResourceStream::SrcImageResourceStream(stream_type* stream, ImageFormat fmt)
  : m_stream(stream, noop_istream), m_fmt(fmt) {
  VW_ASSERT(m_fmt.complete(), ArgumentErr() << "ImageFormat must fully describe the image data");
}
SrcImageResourceStream::SrcImageResourceStream(boost::shared_ptr<stream_type> stream)
  : m_stream(stream) {}
SrcImageResourceStream::SrcImageResourceStream(boost::shared_ptr<stream_type> stream, ImageFormat fmt)
  : m_stream(stream), m_fmt(fmt) {
  VW_ASSERT(m_fmt.complete(), ArgumentErr() << "ImageFormat must fully describe the image data");
}

int32 SrcImageResourceStream::cols() const {
  VW_ASSERT(m_fmt.complete(), LogicErr() << "Function only callable on complete image format " << VW_CURRENT_FUNCTION);
  return m_fmt.cols;
}

int32 SrcImageResourceStream::rows() const {
  VW_ASSERT(m_fmt.complete(), LogicErr() << "Function only callable on complete image format " << VW_CURRENT_FUNCTION);
  return m_fmt.rows;
}

int32 SrcImageResourceStream::planes() const {
  VW_ASSERT(m_fmt.complete(), LogicErr() << "Function only callable on complete image format " << VW_CURRENT_FUNCTION);
  return m_fmt.planes;
}

PixelFormatEnum SrcImageResourceStream::pixel_format() const {
  VW_ASSERT(m_fmt.complete(), LogicErr() << "Function only callable on complete image format " << VW_CURRENT_FUNCTION);
  return m_fmt.pixel_format;
}

ChannelTypeEnum SrcImageResourceStream::channel_type() const {
  VW_ASSERT(m_fmt.complete(), LogicErr() << "Function only callable on complete image format " << VW_CURRENT_FUNCTION);
  return m_fmt.channel_type;
}

const ImageFormat& SrcImageResourceStream::read_format() const {
  return m_fmt;
}

void SrcImageResourceStream::set_read_format(const ImageFormat& fmt) {
  VW_ASSERT(m_fmt.complete(), ArgumentErr() << "ImageFormat must fully describe the image data");
  m_fmt = fmt;
}

// Read the image resource at the given location into the given buffer.
void SrcImageResourceStream::read( ImageBuffer const& dst_buf, BBox2i const& bbox ) const {
  VW_ASSERT(dst_buf.format.cols == uint32(bbox.width()) && dst_buf.format.rows == uint32(bbox.height()),
      LogicErr() << VW_CURRENT_FUNCTION << ": Only complete reads are supported" );

  VW_ASSERT(!m_stream->fail(), IOErr() << "Can't read from stream (the bad or fail flag is already up)");

  // If we're doing unformatted reads, or conversion is simple
  if (!m_fmt.complete() || (m_fmt.complete() && m_fmt.same_size(dst_buf.format) && m_fmt.simple_convert(dst_buf.format))) {
    perform_read(m_stream.get(), reinterpret_cast<char*>(dst_buf.data), dst_buf.byte_size());
    return;
  }

  // Do the complex conversion
  ImageBuffer tmp_src(m_fmt, NULL);
  boost::scoped_array<char> data(new char[tmp_src.byte_size()]);
  tmp_src.data = data.get();

  perform_read(m_stream.get(), reinterpret_cast<char*>(tmp_src.data), tmp_src.byte_size());

  convert(dst_buf, tmp_src, false);
}

void SrcImageResourceStream::reset() {
  VW_ASSERT(!m_stream->fail(), IOErr() << "Can't seek in stream (the bad or fail flag is already up)");
  m_stream->seekg(0);
  VW_ASSERT(!m_stream->fail(), IOErr() << "Failed to seek. Is this input stream seekable?");
}

DstImageResourceStream::DstImageResourceStream(stream_type* stream)
  : m_stream(stream, noop_ostream) {}
DstImageResourceStream::DstImageResourceStream(boost::shared_ptr<stream_type> stream)
  : m_stream(stream) {}
DstImageResourceStream::DstImageResourceStream(stream_type* stream, ImageFormat fmt)
  : m_stream(stream, noop_ostream), m_fmt(fmt) {
  VW_ASSERT(m_fmt.complete(), ArgumentErr() << "ImageFormat must fully describe the image data");
}
DstImageResourceStream::DstImageResourceStream(boost::shared_ptr<stream_type> stream, ImageFormat fmt)
  : m_stream(stream), m_fmt(fmt) {
  VW_ASSERT(m_fmt.complete(), ArgumentErr() << "ImageFormat must fully describe the image data");
}

const ImageFormat& DstImageResourceStream::write_format() const {
  return m_fmt;
}

void DstImageResourceStream::set_write_format(const ImageFormat& fmt) {
  VW_ASSERT(m_fmt.complete(), ArgumentErr() << "ImageFormat must fully describe the image data");
  m_fmt = fmt;
}


// Write the given buffer to the image resource at the given location.
void DstImageResourceStream::write( ImageBuffer const& src_buf, BBox2i const& bbox ) {
  VW_ASSERT(src_buf.format.cols == uint32(bbox.width()) && src_buf.format.rows == uint32(bbox.height()),
      LogicErr() << VW_CURRENT_FUNCTION << ": Only complete writes are supported" );

  VW_ASSERT(!m_stream->bad(), IOErr() << "Can't write to stream (the bad flag is already up)");

  // If we're doing unformatted writes, or conversion is simple
  if (!m_fmt.complete() || (m_fmt.complete() && m_fmt.same_size(src_buf.format) && m_fmt.simple_convert(src_buf.format))) {
    perform_write(m_stream.get(), reinterpret_cast<const char*>(src_buf.data), src_buf.byte_size());
    return;
  }

  // Do the complex conversion
  ImageBuffer tmp_dst(m_fmt, NULL);
  boost::scoped_array<char> data(new char[tmp_dst.byte_size()]);
  tmp_dst.data = data.get();

  convert(tmp_dst, src_buf, false);
  perform_write(m_stream.get(), reinterpret_cast<const char*>(tmp_dst.data), tmp_dst.byte_size());
}

void DstImageResourceStream::flush() {
  m_stream->flush();
}

void DstImageResourceStream::reset() {
  VW_ASSERT(!m_stream->bad(), IOErr() << "Can't seek in stream (the bad flag is already up)");
  m_stream->seekp(0);
  VW_ASSERT(!m_stream->bad(), IOErr() << "Failed to seek. Is this output stream seekable?");
}
