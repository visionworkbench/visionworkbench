// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/FileIO/EncodedImageResourceJPEG.h>
#include <vw/FileIO/JpegIO.h>
#include <vw/Core/Debugging.h>

static void noop_deleter(vw::uint8*) {}

namespace vw {

class SrcEncodedImageResourceJPEG::Data : public fileio::detail::JpegIODecompress {
    const uint8* m_data;
    size_t m_len;
  protected:
    virtual void bind() {
      jpeg_mem_src(&m_ctx, reinterpret_cast<unsigned char*>(const_cast<uint8*>(m_data)), m_len);
    }
  public:
    Data(const uint8* buffer, size_t len) : m_data(buffer), m_len(len) {}
};

void SrcEncodedImageResourceJPEG::set_encoded_data(const uint8* buffer, size_t len) {
  m_data.reset(new Data(buffer, len));
  m_data->open();
}

void SrcEncodedImageResourceJPEG::read( ImageBuffer const& dst, BBox2i const& bbox ) const {
  VW_ASSERT( dst.format.cols == size_t(bbox.width()) && dst.format.rows == size_t(bbox.height()),
             ArgumentErr() << VW_CURRENT_FUNCTION << ": Destination buffer has wrong dimensions!" );
  VW_ASSERT( dst.format.cols == size_t(cols()), ArgumentErr() << VW_CURRENT_FUNCTION << ": Read width must match image width");

  size_t start_line = boost::numeric_cast<size_t>(bbox.min().y()),
         stop_line  = boost::numeric_cast<size_t>(bbox.max().y());

  // If we're already deeper in the file than the read wants, we have to reopen
  // the data to seek to the beginning.
  if (start_line < m_data->line())
    m_data->reopen();

  // shared rather than scoped so we get a deleter.
  boost::shared_array<uint8> buf;

  // no conversion necessary?
  bool simple = m_data->fmt().simple_convert(dst.format);

  // If we don't need to convert, we read directly into the dst buffer (using a
  // noop_deleter, so the destructor doesn't try to delete it)
  if (simple)
    buf.reset( reinterpret_cast<uint8*>(const_cast<void*>(dst.data)), noop_deleter );
  else
    buf.reset( new uint8[m_data->line_bytes() * bbox.height()] );

  // skip lines
  for (size_t i = m_data->line(); i < start_line; ++i)
    m_data->readline(buf.get());

  uint8* current = buf.get();
  for (size_t i = start_line; i < stop_line; ++i) {
    m_data->readline(current);
    current += m_data->line_bytes();
  }

  m_data->flush();

  if (simple)
    return;

  ImageFormat src_fmt(m_data->fmt());
  src_fmt.rows = bbox.height();
  src_fmt.cols = bbox.width();

  ImageBuffer src(src_fmt, buf.get());
  convert( dst, src );
}

int32 SrcEncodedImageResourceJPEG::cols() const {
  return m_data->fmt().cols;
}

int32 SrcEncodedImageResourceJPEG::rows() const {
  return m_data->fmt().rows;
}

int32 SrcEncodedImageResourceJPEG::planes() const {
  return m_data->fmt().planes;
}

PixelFormatEnum SrcEncodedImageResourceJPEG::pixel_format() const {
  return m_data->fmt().pixel_format;
}

ChannelTypeEnum SrcEncodedImageResourceJPEG::channel_type() const {
  return m_data->fmt().channel_type;
}

class DstEncodedImageResourceJPEG::Data : public fileio::detail::JpegIOCompress {
  std::vector<uint8>* m_data;

  struct vector_dest_mgr {
    struct jpeg_destination_mgr pub;
    // Make a guess that we're going to output more than 4k.
    static const size_t BLOCK_SIZE = 4096;
    std::vector<uint8>* vec;

    static void init_destination (j_compress_ptr cinfo) {
      vector_dest_mgr *dest = reinterpret_cast<vector_dest_mgr*>(cinfo->dest);
      dest->vec->resize(BLOCK_SIZE);
      dest->pub.free_in_buffer = BLOCK_SIZE;
      dest->pub.next_output_byte = &(dest->vec->operator[](0));
    }
    static ::boolean empty_output_buffer (j_compress_ptr cinfo) {
      vector_dest_mgr *dest = reinterpret_cast<vector_dest_mgr*>(cinfo->dest);
      // resize gives us exactly what we ask for, but we actually want to do an
      // amortized O(1) insert. So resize to twice the current capacity each
      // time. Normally, size() starts at zero and this wouldn't work, but we
      // initialize it to nonzero above.
      dest->pub.free_in_buffer = dest->vec->size();
      dest->vec->resize(dest->vec->capacity() * 2);
      dest->pub.next_output_byte = &(dest->vec->operator[](dest->pub.free_in_buffer));
      return true;
    }
    static void term_destination (j_compress_ptr cinfo) {
      vector_dest_mgr *dest = reinterpret_cast<vector_dest_mgr*>(cinfo->dest);
      // cut off the default-initialized tail
      dest->vec->resize(dest->vec->size() - dest->pub.free_in_buffer);
    }
  };

  static void std_vector_dest(j_compress_ptr cinfo, std::vector<uint8>* v) {
    VW_ASSERT(v, ArgumentErr() << "std_vector_dest: Expected a non-null vector");

    // allocate the mgr using libjpeg's current memory manager. It has the same
    // lifetime as the j_compress_ptr
    vector_dest_mgr *dest = reinterpret_cast<vector_dest_mgr*>(
        (*cinfo->mem->alloc_small)(reinterpret_cast<j_common_ptr>(cinfo), JPOOL_PERMANENT, sizeof(vector_dest_mgr)));

    dest->pub.init_destination    = &vector_dest_mgr::init_destination;
    dest->pub.empty_output_buffer = &vector_dest_mgr::empty_output_buffer;
    dest->pub.term_destination    = &vector_dest_mgr::term_destination;
    dest->vec                 = v;
    cinfo->dest               = reinterpret_cast<jpeg_destination_mgr*>(dest);
  }


  protected:
    virtual void bind() {
      std_vector_dest(&m_ctx, m_data);
    }
  public:
    Data(std::vector<uint8>* buffer, const ImageFormat &fmt) : JpegIOCompress(fmt), m_data(buffer) {}
};


DstEncodedImageResourceJPEG::DstEncodedImageResourceJPEG(std::vector<uint8>* buffer, const ImageFormat& fmt)
  : m_data(new Data(buffer, fmt))
{ this->set_decoded_data(buffer); }

void DstEncodedImageResourceJPEG::set_decoded_data( std::vector<uint8>* buffer) {
  m_data.reset(new Data(buffer, m_data->fmt()));
  m_data->open();
}

void DstEncodedImageResourceJPEG::write( ImageBuffer const& src, BBox2i const& bbox ) {
  VW_ASSERT( src.format.cols == size_t(bbox.width()) && src.format.rows == size_t(bbox.height()),
             ArgumentErr() << VW_CURRENT_FUNCTION << ": partial writes not supported." );

  m_data->startwrite(bbox.width(), bbox.height());

  size_t start_line = boost::numeric_cast<size_t>(bbox.min().y()),
         stop_line  = boost::numeric_cast<size_t>(bbox.max().y());

  // shared rather than scoped so we get a deleter.
  boost::shared_array<uint8> buf;

  // no conversion necessary?
  bool simple = src.format.simple_convert(m_data->fmt());

  // If we don't need to convert, we write directly from the src buffer (using a
  // noop_deleter, so the destructor doesn't try to delete it)
  if (simple)
    buf.reset( reinterpret_cast<uint8*>(const_cast<void*>(src.data)), noop_deleter );
  else {
    buf.reset( new uint8[m_data->line_bytes() * bbox.height()] );

    ImageFormat dst_fmt(m_data->fmt());
    dst_fmt.rows = bbox.height();
    dst_fmt.cols = bbox.width();

    ImageBuffer dst(dst_fmt, buf.get());
    convert( dst, src );
  }

  uint8* current = buf.get();
  for (size_t i = start_line; i < stop_line; ++i) {
    m_data->writeline(current);
    current += m_data->line_bytes();
  }

  m_data->flush();
}

} // namespace vw
