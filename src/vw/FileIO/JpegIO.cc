#include <vw/FileIO/JpegIO.h>
#include <vw/Core/Exception.h>
#include <vw/Core/Log.h>

static void vw_jpeg_error_exit(j_common_ptr cinfo) {
  char buffer[JMSG_LENGTH_MAX];
  (*cinfo->err->format_message)(cinfo, buffer);
  int msg_code = cinfo->err->msg_code;
  jpeg_destroy(cinfo);
  if ( msg_code == JERR_NO_SOI )
    vw::vw_throw( vw::ArgumentErr() << "JpegIO: Cannot open non-jpeg files.\n" );
  else
    vw::vw_throw( vw::IOErr() << "JpegIO error: " << buffer );
}


namespace vw {
namespace fileio {
namespace detail {

////////////////////////////////////////////////////////////////////////////////
// Base class functions
////////////////////////////////////////////////////////////////////////////////
void JpegIO::init_base(jpeg_error_mgr** mgr) {
  *mgr = jpeg_std_error(&m_jerr);
  m_jerr.error_exit = &vw_jpeg_error_exit;
}

////////////////////////////////////////////////////////////////////////////////
// Decompress
////////////////////////////////////////////////////////////////////////////////
JpegIODecompress::JpegIODecompress()
{
  init_base(&m_ctx.err);
  jpeg_create_decompress(&m_ctx);
}

JpegIODecompress::~JpegIODecompress() {
  jpeg_destroy_decompress(&m_ctx);
}

bool JpegIODecompress::ready() const {
  return m_ctx.output_scanline == 0;
}

void JpegIODecompress::read(uint8* buffer, size_t bufsize) {
  VW_ASSERT(this->ready(), LogicErr() << "Cannot reread from a JpegIO reader");

  jpeg_start_decompress(&m_ctx);
  size_t skip = line_bytes();
  VW_ASSERT(bufsize >= m_ctx.output_height * skip, LogicErr() << "Buffer is too small");

  while (m_ctx.output_scanline < m_ctx.output_height) {
    jpeg_read_scanlines(&m_ctx, &buffer, 1);
    buffer += skip;
  }
  jpeg_finish_decompress(&m_ctx);
}

void JpegIODecompress::open() {
  bind();
  jpeg_read_header(&m_ctx, TRUE);
  jpeg_calc_output_dimensions(&m_ctx);

  m_fmt.cols = m_ctx.output_width;
  m_fmt.rows = m_ctx.output_height;
  m_fmt.channel_type = VW_CHANNEL_UINT8;

  switch (m_ctx.output_components)
  {
    case 1:
      m_fmt.pixel_format = VW_PIXEL_GRAY;
      m_fmt.planes=1;
      break;
    case 2:
      m_fmt.pixel_format = VW_PIXEL_GRAYA;
      m_fmt.planes=1;
      break;
    case 3:
      m_fmt.pixel_format = VW_PIXEL_RGB;
      m_fmt.planes=1;
      break;
    case 4:
      m_fmt.pixel_format = VW_PIXEL_RGBA;
      m_fmt.planes=1; break;
    default:
      m_fmt.planes = m_ctx.output_components;
      m_fmt.pixel_format = VW_PIXEL_SCALAR;
      break;
  }

  m_cstride = m_ctx.output_components;
  m_rstride = m_cstride * m_fmt.cols;
}

////////////////////////////////////////////////////////////////////////////////
// Compress
////////////////////////////////////////////////////////////////////////////////
JpegIOCompress::JpegIOCompress(const ImageFormat& fmt)
{
  m_fmt = fmt;
  init_base(&m_ctx.err);
  jpeg_create_compress(&m_ctx);
}

JpegIOCompress::~JpegIOCompress() {
  jpeg_destroy_compress(&m_ctx);
}

bool JpegIOCompress::ready() const {
  return m_ctx.next_scanline == 0;
}

void JpegIOCompress::write(const uint8* data, size_t rows, size_t cols, size_t bufsize) {
  VW_ASSERT(this->ready(), LogicErr() << "Cannot rewrite to a JpegIO writer");

  size_t skip = cols * chan_bytes();
  VW_ASSERT(bufsize >= rows * skip, LogicErr() << "Buffer is too small");

  m_ctx.image_width  = cols;
  m_ctx.image_height = rows;
  jpeg_start_compress(&m_ctx, TRUE);

  JSAMPLE* jdata = reinterpret_cast<JSAMPLE*>(const_cast<uint8*>(data));
  while (m_ctx.next_scanline < m_ctx.image_height) {
    jpeg_write_scanlines(&m_ctx, &jdata, 1);
    jdata += skip;
  }
  jpeg_finish_compress(&m_ctx);
}

void JpegIOCompress::open() {
  bind();

  switch (fmt().pixel_format) {
    case VW_PIXEL_SCALAR:
      m_ctx.input_components = fmt().planes ? fmt().planes : 1;
      m_ctx.in_color_space = JCS_UNKNOWN;
      break;
    case VW_PIXEL_GRAYA:
      vw_out(DebugMessage, "fileio") << "JpegIOCompress: Warning: alpha channel removed." << std::endl;;
    case VW_PIXEL_GRAY:
      m_ctx.input_components = 1;
      m_ctx.in_color_space = JCS_GRAYSCALE;
      break;
    case VW_PIXEL_RGBA:
      vw_out(DebugMessage, "fileio") << "JpegIOCompress: Warning: alpha channel removed." << std::endl;;
    case VW_PIXEL_RGB:
      m_ctx.input_components = 3;
      m_ctx.in_color_space = JCS_RGB;
      break;
    default:
        vw_throw( IOErr() << "JpegIOCompress: Unsupported pixel type (" << fmt().pixel_format << ")." );
  }

  m_cstride = m_ctx.input_components;

  jpeg_set_defaults(&m_ctx);
  jpeg_set_quality(&m_ctx, 95, TRUE); // XXX Make this settable
}

////////////////////////////////////////////////////////////////////////////////
// Src/Dest Helpers
////////////////////////////////////////////////////////////////////////////////

struct ptr_src_mgr {
  struct jpeg_source_mgr pub;

  static void init_source(j_decompress_ptr /*cinfo*/) {/* no work */}
  static void term_source(j_decompress_ptr /*cinfo*/) {/* no work */}

  static ::boolean fill_input_buffer(j_decompress_ptr /*cinfo*/) {
#if 0
    // No more data to provide... this is what libjpeg does (insert a fake EOI).
    static JOCTET buf[4];
    buf[0] = (JOCTET)0xFF;
    buf[1] = (JOCTET)JPEG_EOI;
    cinfo->src->next_input_byte = buf;
    cinfo->src->bytes_in_buffer = 2;
    vw_log(WarningMessage) << "Damaged JPEG" << std::endl;
#endif
    // I think we should be more strict.
    vw_throw(IOErr() << "Damaged JPEG. No EOI? Cannot continue.");
    return TRUE;
  }

  static void skip_input_data (j_decompress_ptr cinfo, long num_bytes_l)
  {
    VW_ASSERT(num_bytes_l >= 0, ArgumentErr() << "Cannot skip negative bytes");
    size_t num_bytes = num_bytes_l;

    if (num_bytes == 0)
      return;

    ptr_src_mgr *src = reinterpret_cast<ptr_src_mgr*>(cinfo->src);
    VW_ASSERT(num_bytes <= src->pub.bytes_in_buffer, ArgumentErr() << "Cannot skip more bytes than are left");

    src->pub.next_input_byte += num_bytes;
    src->pub.bytes_in_buffer -= num_bytes;
  }
};

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

void jpeg_ptr_src(j_decompress_ptr cinfo, const uint8* data, size_t size) {
  VW_ASSERT(data, ArgumentErr() << "jpeg_ptr_src: Expected a non-null data ptr");
  VW_ASSERT(size, ArgumentErr() << "jpeg_ptr_src: Expected a non-zero size");

  ptr_src_mgr *src = reinterpret_cast<ptr_src_mgr*>(
      (*cinfo->mem->alloc_small)(reinterpret_cast<j_common_ptr>(cinfo), JPOOL_PERMANENT, sizeof(ptr_src_mgr)));

  src->pub.init_source       = &ptr_src_mgr::init_source;
  src->pub.fill_input_buffer = &ptr_src_mgr::fill_input_buffer;
  src->pub.skip_input_data   = &ptr_src_mgr::skip_input_data;
  src->pub.resync_to_restart = ::jpeg_resync_to_restart;
  src->pub.term_source       = &ptr_src_mgr::term_source;
  src->pub.bytes_in_buffer   = size;
  src->pub.next_input_byte   = reinterpret_cast<JOCTET*>(const_cast<uint8*>(data));
  cinfo->src                 = reinterpret_cast<jpeg_source_mgr*>(src);
}

void jpeg_vector_dest(j_compress_ptr cinfo, std::vector<uint8>* v) {
  VW_ASSERT(v, ArgumentErr() << "std_vector_dest: Expected a non-null vector");

  // allocate the mgr using libjpeg's current memory manager. It has the same
  // lifetime as the j_compress_ptr
  vector_dest_mgr *dest = reinterpret_cast<vector_dest_mgr*>(
      (*cinfo->mem->alloc_small)(reinterpret_cast<j_common_ptr>(cinfo), JPOOL_PERMANENT, sizeof(vector_dest_mgr)));

  dest->pub.init_destination    = &vector_dest_mgr::init_destination;
  dest->pub.empty_output_buffer = &vector_dest_mgr::empty_output_buffer;
  dest->pub.term_destination    = &vector_dest_mgr::term_destination;
  dest->vec                     = v;
  cinfo->dest                   = reinterpret_cast<jpeg_destination_mgr*>(dest);
}

}}} // namespace vw::fileio::detail
