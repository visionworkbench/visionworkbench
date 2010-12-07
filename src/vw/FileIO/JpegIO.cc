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

void JpegIO::reopen() {
  this->close();
  this->open();
}

const ImageFormat& JpegIO::fmt() const {
  return m_fmt;
}

size_t JpegIO::chan_bytes() { return m_cstride; }
size_t JpegIO::line_bytes() { return m_rstride; }


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

void JpegIODecompress::flush() {
  jpeg_finish_decompress(&m_ctx);
}
void JpegIODecompress::close() {
  jpeg_abort_decompress(&m_ctx);
}
size_t JpegIODecompress::line() const {
  return m_ctx.output_scanline;
}
void JpegIODecompress::readline(uint8* data) {
  jpeg_read_scanlines(&m_ctx, &data, 1);
}

void JpegIODecompress::open() {
  bind();
  jpeg_read_header(&m_ctx, TRUE);
  jpeg_start_decompress(&m_ctx);

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

void JpegIOCompress::flush() {
  jpeg_finish_compress(&m_ctx);
}
void JpegIOCompress::close() {
  jpeg_abort_compress(&m_ctx);
}
size_t JpegIOCompress::line() const {
  return m_ctx.next_scanline;
}

void JpegIOCompress::startwrite(uint32 image_width, uint32 image_height) {
  m_ctx.image_width  = image_width;
  m_ctx.image_height = image_height;
  m_rstride = m_cstride * image_width;
  jpeg_start_compress(&m_ctx, TRUE);
}

void JpegIOCompress::writeline(uint8* data) {
  JSAMPLE* jdata = reinterpret_cast<JSAMPLE*>(data);
  jpeg_write_scanlines(&m_ctx, &jdata, 1);
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
  jpeg_set_quality(&m_ctx, 100, TRUE); // XXX Make this settable
}


}}} // namespace vw::fileio::detail
