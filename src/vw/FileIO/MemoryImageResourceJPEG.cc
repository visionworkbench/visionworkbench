// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/FileIO/MemoryImageResourceJPEG.h>
#include <vw/FileIO/JpegIO.h>
#include <vw/Core/Debugging.h>

static void noop_deleter(vw::uint8*) {}

namespace vw {

class SrcMemoryImageResourceJPEG::Data : public fileio::detail::JpegIODecompress {
    const uint8* m_data;
    size_t m_len;
  protected:
    virtual void bind() { fileio::detail::jpeg_ptr_src(&m_ctx, m_data, m_len); }
  public:
    Data* rewind() const VW_WARN_UNUSED {std::auto_ptr<Data> r(new Data(m_data, m_len)); r->open(); return r.release();}
    Data(const uint8* buffer, size_t len) : m_data(buffer), m_len(len) {}
};

SrcMemoryImageResourceJPEG::SrcMemoryImageResourceJPEG(const uint8* buffer, size_t len)
  : m_data(new Data(buffer, len)) {
    m_data->open();
}

void SrcMemoryImageResourceJPEG::read( ImageBuffer const& dst, BBox2i const& bbox ) const {
  VW_ASSERT( dst.format.cols == size_t(bbox.width()) && dst.format.rows == size_t(bbox.height()),
             ArgumentErr() << VW_CURRENT_FUNCTION << ": Destination buffer has wrong dimensions!" );
  VW_ASSERT( dst.format.cols == size_t(cols()) && dst.format.rows == size_t(rows()),
             ArgumentErr() << VW_CURRENT_FUNCTION << ": Partial reads are not supported");

  if (!m_data->ready())
    m_data.reset(m_data->rewind());

  // shared rather than scoped so we get a deleter.
  boost::shared_array<uint8> buf;

  // no conversion necessary?
  bool simple = m_data->fmt().simple_convert(dst.format);

  size_t bufsize = m_data->line_bytes() * bbox.height();

  // If we don't need to convert, we read directly into the dst buffer (using a
  // noop_deleter, so the destructor doesn't try to delete it)
  if (simple)
    buf.reset( reinterpret_cast<uint8*>(const_cast<void*>(dst.data)), noop_deleter );
  else
    buf.reset( new uint8[bufsize] );

  m_data->read(buf.get(), bufsize);

  if (simple)
    return;

  ImageFormat src_fmt(m_data->fmt());
  src_fmt.rows = bbox.height();
  src_fmt.cols = bbox.width();

  ImageBuffer src(src_fmt, buf.get());
  convert(dst, src, true);
}

int32 SrcMemoryImageResourceJPEG::cols() const {
  return m_data->fmt().cols;
}

int32 SrcMemoryImageResourceJPEG::rows() const {
  return m_data->fmt().rows;
}

int32 SrcMemoryImageResourceJPEG::planes() const {
  return m_data->fmt().planes;
}

PixelFormatEnum SrcMemoryImageResourceJPEG::pixel_format() const {
  return m_data->fmt().pixel_format;
}

ChannelTypeEnum SrcMemoryImageResourceJPEG::channel_type() const {
  return m_data->fmt().channel_type;
}

class DstMemoryImageResourceJPEG::Data : public fileio::detail::JpegIOCompress {
  std::vector<uint8>* m_data;

  protected:
    virtual void bind() { fileio::detail::jpeg_vector_dest(&m_ctx, m_data); }
  public:
    Data(std::vector<uint8>* buffer, const ImageFormat &fmt) : JpegIOCompress(fmt), m_data(buffer) {}
};


DstMemoryImageResourceJPEG::DstMemoryImageResourceJPEG(std::vector<uint8>* buffer, const ImageFormat& fmt)
  : m_data(new Data(buffer, fmt))
{
  m_data->open();
}

void DstMemoryImageResourceJPEG::write( ImageBuffer const& src, BBox2i const& bbox ) {
  size_t width = bbox.width(), height = bbox.height();
  VW_ASSERT( src.format.cols == width && src.format.rows == height,
             ArgumentErr() << VW_CURRENT_FUNCTION << ": partial writes not supported." );
  VW_ASSERT(m_data->ready(), LogicErr() << "Multiple writes to one DstMemoryImageResourceJPEG. Probably a bug?");

  // shared rather than scoped so we get a deleter.
  boost::shared_array<uint8> buf;

  // no conversion necessary?
  bool simple = src.format.simple_convert(m_data->fmt());

  size_t bufsize = m_data->chan_bytes() * width * height;

  // If we don't need to convert, we write directly from the src buffer (using a
  // noop_deleter, so the destructor doesn't try to delete it)
  if (simple)
    buf.reset( reinterpret_cast<uint8*>(const_cast<void*>(src.data)), noop_deleter );
  else {
    buf.reset( new uint8[bufsize] );

    ImageFormat dst_fmt(m_data->fmt());
    dst_fmt.rows = height;
    dst_fmt.cols = width;

    ImageBuffer dst(dst_fmt, buf.get());
    convert(dst, src, true);
  }

  m_data->write(buf.get(), width, height, bufsize);
}

} // namespace vw
