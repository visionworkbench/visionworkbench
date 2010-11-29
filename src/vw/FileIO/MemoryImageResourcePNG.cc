// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/FileIO/MemoryImageResourcePNG.h>
#include <vw/FileIO/PngIO.h>
#include <vw/Core/Debugging.h>

static void noop_deleter(vw::uint8*) {}

namespace vw {

class SrcMemoryImageResourcePNG::Data : public fileio::detail::PngIODecompress {
    const uint8 * const m_begin, *m_cur;
    const uint8 *const m_end;
  protected:
    virtual void bind() { png_set_read_fn(m_ctx, reinterpret_cast<voidp>(this), &SrcMemoryImageResourcePNG::Data::read_fn); }
  public:
    Data* rewind() const VW_WARN_UNUSED {std::auto_ptr<Data> r(new Data(m_begin, m_end-m_begin)); r->open(); return r.release();}

    static void read_fn( png_structp ctx, png_bytep data, png_size_t len ) {
      Data *mgr = reinterpret_cast<Data*>(png_get_io_ptr(ctx));
      VW_ASSERT(mgr->m_cur + len <= mgr->m_end, LogicErr() << "No more data!");
      std::copy(mgr->m_cur, mgr->m_cur+len, data);
      mgr->m_cur += len;
    }
    Data(const uint8* buffer, size_t len) : m_begin(buffer), m_cur(buffer), m_end(buffer+len) {}
};

SrcMemoryImageResourcePNG::SrcMemoryImageResourcePNG(const uint8* buffer, size_t len)
  : m_data(new Data(buffer, len)) {
    m_data->open();
}

void SrcMemoryImageResourcePNG::read( ImageBuffer const& dst, BBox2i const& bbox ) const {
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

int32 SrcMemoryImageResourcePNG::cols() const {
  return m_data->fmt().cols;
}

int32 SrcMemoryImageResourcePNG::rows() const {
  return m_data->fmt().rows;
}

int32 SrcMemoryImageResourcePNG::planes() const {
  return m_data->fmt().planes;
}

PixelFormatEnum SrcMemoryImageResourcePNG::pixel_format() const {
  return m_data->fmt().pixel_format;
}

ChannelTypeEnum SrcMemoryImageResourcePNG::channel_type() const {
  return m_data->fmt().channel_type;
}

class DstMemoryImageResourcePNG::Data : public fileio::detail::PngIOCompress {
  std::vector<uint8>* m_data;
  typedef DstMemoryImageResourcePNG::Data this_type;

  protected:
    virtual void bind() {
      png_set_write_fn(m_ctx, reinterpret_cast<voidp>(this), &this_type::write_fn, &this_type::flush_fn);
    }
  public:
    Data(std::vector<uint8>* buffer, const ImageFormat &fmt) : PngIOCompress(fmt), m_data(buffer) {}

    static void write_fn( png_structp ctx, png_bytep data, png_size_t length )
    {
      Data *mgr = reinterpret_cast<Data*>(png_get_io_ptr(ctx));
      mgr->m_data->insert(mgr->m_data->end(), data, data+length);
    }

    static void flush_fn( png_structp /*ctx*/) {}
};


DstMemoryImageResourcePNG::DstMemoryImageResourcePNG(std::vector<uint8>* buffer, const ImageFormat& fmt)
  : m_data(new Data(buffer, fmt))
{
  m_data->open();
}

void DstMemoryImageResourcePNG::write( ImageBuffer const& src, BBox2i const& bbox ) {
  size_t width = bbox.width(), height = bbox.height();
  VW_ASSERT( src.format.cols == width && src.format.rows == height,
             ArgumentErr() << VW_CURRENT_FUNCTION << ": partial writes not supported." );
  VW_ASSERT(m_data->ready(), LogicErr() << "Multiple writes to one DstMemoryImageResourcePNG. Probably a bug?");

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
