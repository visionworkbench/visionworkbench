#include <vw/FileIO/PngIO.h>
#include <vw/Core/Exception.h>
#include <vw/Core/Log.h>
#include <vw/Core/Settings.h>

static void png_error_handler(png_structp /*png_ptr*/, png_const_charp error_msg)
{
  vw::vw_throw(vw::IOErr() << "PngIO Error: " << error_msg);
}


namespace vw {
namespace fileio {
namespace detail {

////////////////////////////////////////////////////////////////////////////////
// Decompress
////////////////////////////////////////////////////////////////////////////////
PngIODecompress::PngIODecompress()
  : m_read(false) {}

PngIODecompress::~PngIODecompress() {
  png_destroy_read_struct(&m_ctx, &m_info, NULL);
}

bool PngIODecompress::ready() const {
  return !m_read;
}

void PngIODecompress::read(uint8* buffer, size_t bufsize) {
  VW_ASSERT(this->ready(), LogicErr() << "PngIO: Cannot reread");

  size_t skip = line_bytes();
  VW_ASSERT(bufsize >= m_fmt.rows * skip, LogicErr() << "Buffer is too small");

  boost::scoped_array<png_bytep> rows( new png_bytep[m_fmt.rows] );
  for(size_t i=0; i < m_fmt.rows; ++i)
    rows[i] = reinterpret_cast<png_bytep>(buffer + i * m_rstride);

  m_read = true;
  png_read_image(m_ctx, rows.get());
}

void PngIODecompress::open() {
  m_ctx = png_create_read_struct(PNG_LIBPNG_VER_STRING, 0, png_error_handler, NULL);
  VW_ASSERT(m_ctx, IOErr() << "Failed to create read struct");
  m_info = png_create_info_struct(m_ctx);
  if (!m_info) {
    png_destroy_read_struct(&m_ctx, NULL, NULL);
    vw_throw(IOErr() << "Failed to create info struct");
  }

  this->bind();

  png_uint_32 width, height;
  int bit_depth, color_type, interlace_type, compression_type, filter_type;

  png_read_info(m_ctx, m_info);
  png_get_IHDR(m_ctx, m_info, &width, &height, &bit_depth, &color_type, &interlace_type, &compression_type, &filter_type);

  // unpack pngs <8 bit-depth into whole bytes
  if (bit_depth < 8)
    png_set_packing(m_ctx);

  // Expand paletted images to RGB
  // Expand grayscale images of less than 8-bit depth to 8-bit depth
  // Expand tRNS chunks to alpha channels.
  if (bit_depth < 8 || color_type == PNG_COLOR_TYPE_PALETTE || png_get_valid(m_ctx, m_info, PNG_INFO_tRNS))
    png_set_expand(m_ctx);

  // Make sure 16-bit pngs are LSB
  if (bit_depth == 16)
    png_set_swap(m_ctx);

  // Update info with the transforms we applied
  png_read_update_info(m_ctx, m_info);

  // re-fetch information
  png_get_IHDR(m_ctx, m_info, &width, &height, &bit_depth, &color_type, &interlace_type, &compression_type, &filter_type);

  m_fmt.cols = width;
  m_fmt.rows = height;
  m_fmt.planes = 1;

  switch(bit_depth)
  {
    case 1:
    case 2:
    case 4:  vw_throw(vw::LogicErr() << "png_set_expand should have expanded bit_depth < 8!"); break;
    case 8:  m_fmt.channel_type = VW_CHANNEL_UINT8;  break;
    case 16: m_fmt.channel_type = VW_CHANNEL_UINT16; break;
    default: vw_throw(vw::ArgumentErr() << "Unknown bit depth " << bit_depth); break;
  }

  switch(color_type)
  {
    case PNG_COLOR_TYPE_GRAY:       m_fmt.pixel_format = VW_PIXEL_GRAY;  break;
    case PNG_COLOR_TYPE_GRAY_ALPHA: m_fmt.pixel_format = VW_PIXEL_GRAYA; break;
    case PNG_COLOR_TYPE_PALETTE:    vw_throw(LogicErr() << "png_set_expand should have expanded palette!"); break;
    case PNG_COLOR_TYPE_RGB:        m_fmt.pixel_format = VW_PIXEL_RGB;   break;
    case PNG_COLOR_TYPE_RGB_ALPHA:  m_fmt.pixel_format = VW_PIXEL_RGBA;  break;
    default:
      vw_throw(vw::ArgumentErr() << "Unknown color type in png: " << color_type);
  }

  m_cstride = num_channels(m_fmt.pixel_format) * channel_size(m_fmt.channel_type);
  m_rstride = m_cstride * m_fmt.cols;

  size_t png_row = png_get_rowbytes(m_ctx, m_info);
  VW_ASSERT(png_row == m_rstride, LogicErr() << "libpng and vw do not agree on the size of a row");

  // If we did any sample-reading, clear it here
  png_free_data(m_ctx, m_info, PNG_FREE_ROWS, 0);
}

////////////////////////////////////////////////////////////////////////////////
// Compress
////////////////////////////////////////////////////////////////////////////////
PngIOCompress::PngIOCompress(const ImageFormat& fmt) {
  m_fmt = fmt;
  m_written = false;
}

PngIOCompress::~PngIOCompress() {
  png_destroy_write_struct(&m_ctx, &m_info);
}

bool PngIOCompress::ready() const {
  return !m_written;
}

void PngIOCompress::write(const uint8* data, size_t bufsize, size_t rows, size_t cols, size_t planes) {
  VW_ASSERT(this->ready(), LogicErr() << "PngIO: Cannot rewrite");
  VW_ASSERT(planes == 1, LogicErr() << "PNG does not support multi-plane images");

  size_t skip = cols * chan_bytes();
  VW_ASSERT(bufsize >= rows * skip, LogicErr() << "Buffer is too small");

  png_set_compression_level(m_ctx, Z_BEST_SPEED);

  {
    png_uint_32 width, height;
    int bit_depth, color_type, interlace_type, compression_type, filter_method;
    png_get_IHDR(m_ctx, m_info, &width, &height, &bit_depth, &color_type, &interlace_type, &compression_type, &filter_method);
    png_set_IHDR(m_ctx, m_info, cols, rows, bit_depth, color_type, interlace_type, compression_type, filter_method);
  }
  png_write_info(m_ctx, m_info);
  png_set_swap(m_ctx);

  boost::scoped_array<png_bytep> row_pointers( new png_bytep[rows] );

  png_bytep row = reinterpret_cast<png_bytep>(const_cast<uint8*>(data));
  for(size_t i = 0; i < rows; i++) {
    row_pointers[i] = row;
    row += skip;
  }

  png_write_image(m_ctx, row_pointers.get());
  png_write_end(m_ctx, m_info);
}

void PngIOCompress::open() {

  m_ctx = png_create_write_struct(PNG_LIBPNG_VER_STRING, 0, png_error_handler, NULL);
  VW_ASSERT(m_ctx, IOErr() << "Failed to create read struct");
  m_info = png_create_info_struct(m_ctx);
  if (!m_info) {
    png_destroy_read_struct(&m_ctx, NULL, NULL);
    vw_throw(IOErr() << "Failed to create info struct");
  }

  this->bind();

  int bit_depth, color_type, interlace_type, compression_type, filter_method;

  // anything else will be converted to UINT8
  switch (fmt().channel_type) {
    //case VW_CHANNEL_BOOL:
    //  bit_depth = 1;
    //  break;
    case VW_CHANNEL_INT16:
    case VW_CHANNEL_UINT16:
    case VW_CHANNEL_FLOAT16:
    case VW_CHANNEL_GENERIC_2_BYTE:
      bit_depth = 16;
    default:
      bit_depth = 8;
      break;
  }

  color_type = -1;  // Set a default value to avoid compiler warnings
  switch(fmt().pixel_format) {
    case VW_PIXEL_SCALAR: // fall through
    case VW_PIXEL_GRAY:   color_type = PNG_COLOR_TYPE_GRAY;       break;
    case VW_PIXEL_GRAYA:  color_type = PNG_COLOR_TYPE_GRAY_ALPHA; break;
    case VW_PIXEL_RGB:    color_type = PNG_COLOR_TYPE_RGB;        break;
    case VW_PIXEL_RGBA:   color_type = PNG_COLOR_TYPE_RGBA;       break;
    default: vw_throw(vw::ArgumentErr() << "Unsupported pixel format for png: " << fmt().pixel_format);
  }

  interlace_type   = PNG_INTERLACE_NONE;
  compression_type = PNG_COMPRESSION_TYPE_DEFAULT;
  filter_method    = PNG_FILTER_TYPE_DEFAULT;

  png_set_IHDR(m_ctx, m_info, 4, 4, bit_depth, color_type, interlace_type, compression_type, filter_method);

  m_cstride = num_channels(fmt().pixel_format) * channel_size(fmt().channel_type);
}

}}} // namespace vw::fileio::detail
