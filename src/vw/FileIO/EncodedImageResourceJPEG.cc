// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/FileIO/EncodedImageResourceJPEG.h>
#include <vw/Core/Debugging.h>

extern "C" {
#include <jpeglib.h>
#include <jerror.h>
}

static void vw_jpeg_error_exit(j_common_ptr cinfo) {
  char buffer[JMSG_LENGTH_MAX];
  (*cinfo->err->format_message)(cinfo, buffer);
  int msg_code = cinfo->err->msg_code;
  jpeg_destroy(cinfo);
  if ( msg_code == JERR_NO_SOI )
    vw::vw_throw( vw::ArgumentErr() << "EncodedImageResourceJPEG: Cannot open non-jpeg files.\n" );
  else
    vw::vw_throw( vw::IOErr() << "EncodedImageResourceJPEG error: " << buffer );
}

static void null_deleter(vw::uint8*) {}

namespace vw {

struct SrcEncodedImageResourceJPEG::Data {
  private:
    jpeg_error_mgr jerr;
    size_t cstride;
    jpeg_decompress_struct ctx;

    void init() {
      jpeg_read_header(&ctx, TRUE);
      jpeg_start_decompress(&ctx);

      // Set the parent class's format.
      fmt.cols = ctx.output_width;
      fmt.rows = ctx.output_height;
      fmt.channel_type = VW_CHANNEL_UINT8;

      switch (ctx.output_components)
      {
        case 1:
          fmt.pixel_format = VW_PIXEL_GRAY;
          fmt.planes=1;
          break;
        case 2:
          fmt.pixel_format = VW_PIXEL_GRAYA;
          fmt.planes=1;
          break;
        case 3:
          fmt.pixel_format = VW_PIXEL_RGB;
          fmt.planes=1;
          break;
        case 4:
          fmt.pixel_format = VW_PIXEL_RGBA;
          fmt.planes=1; break;
        default:
          fmt.planes = ctx.output_components;
          fmt.pixel_format = VW_PIXEL_SCALAR;
          break;
      }

      cstride = ctx.output_components;
      line = 0;
    }

    void fini() {
      jpeg_abort_decompress(&ctx);
    }

  public:
    ImageFormat fmt;
    size_t line;

    Data(const uint8* data, size_t len) {
      ctx.err = jpeg_std_error(&jerr);
      jerr.error_exit = &vw_jpeg_error_exit;

      jpeg_create_decompress(&ctx);
      jpeg_mem_src(&ctx, reinterpret_cast<unsigned char*>(const_cast<uint8*>(data)), len);
      init();
    }

    size_t byte_size(size_t pixels) {
      return cstride * pixels;
    }

    void readline(uint8* data) {
      jpeg_read_scanlines(&ctx, &data, 1);
      line++;
    }

    void restart() {
      fini();
      init();
    }

    ~Data() {
      jpeg_destroy_decompress(&ctx);
    }
};

void SrcEncodedImageResourceJPEG::set_encoded_data(const uint8* data, size_t len) {
  m_encoded = data;
  m_size = len;
  m_data.reset(new Data(data, len));
}

void SrcEncodedImageResourceJPEG::read( ImageBuffer const& dst, BBox2i const& bbox ) const {
  VW_ASSERT( dst.format.cols == size_t(bbox.width()) && dst.format.rows == size_t(bbox.height()),
             ArgumentErr() << VW_CURRENT_FUNCTION << ": Destination buffer has wrong dimensions!" );
  VW_ASSERT( dst.format.cols == size_t(cols()), ArgumentErr() << VW_CURRENT_FUNCTION << ": Read width must match image width");

  size_t start_line = boost::numeric_cast<size_t>(bbox.min().y()),
         stop_line  = boost::numeric_cast<size_t>(bbox.max().y());

  if (start_line > m_data->line)
    m_data->restart();

  // shared rather than scoped so we get a deleter.
  boost::shared_array<uint8> buf;

  // no conversion necessary?
  bool simple = m_data->fmt.simple_convert(dst.format);

  // If we don't need to convert, we read directly into the dst buffer
  if (simple)
    buf.reset( reinterpret_cast<uint8*>(const_cast<void*>(dst.data)), null_deleter );
  else
    buf.reset( new uint8[m_data->byte_size(1) * bbox.width() * bbox.height()] );

  // skip lines
  for (size_t i = m_data->line; i < start_line; ++i)
    m_data->readline(buf.get());

  uint8* current = buf.get();
  for (size_t i = start_line; i < stop_line; ++i) {
    m_data->readline(current);
    current += m_data->byte_size(bbox.width());
  }

  if (simple)
    return;

  ImageFormat src_fmt(m_data->fmt);
  src_fmt.rows = bbox.height();
  src_fmt.cols = bbox.width();

  ImageBuffer src(src_fmt, buf.get());
  convert( dst, src );
}

int32 SrcEncodedImageResourceJPEG::cols() const {
  return m_data->fmt.cols;
}

int32 SrcEncodedImageResourceJPEG::rows() const {
  return m_data->fmt.rows;
}

int32 SrcEncodedImageResourceJPEG::planes() const {
  return m_data->fmt.planes;
}

PixelFormatEnum SrcEncodedImageResourceJPEG::pixel_format() const {
  return m_data->fmt.pixel_format;
}

ChannelTypeEnum SrcEncodedImageResourceJPEG::channel_type() const {
  return m_data->fmt.channel_type;
}

#if 0
struct SrcEncodedImageResourceJPEG::Data {

  private:
    jpeg_error_mgr jerr;
    size_t cstride;
    jpeg_decompress_struct ctx;

    void init() {
      jpeg_read_header(&ctx, TRUE);
      jpeg_start_decompress(&ctx);

      fmt.cols = ctx.output_width;
      fmt.rows = ctx.output_height;
      fmt.channel_type = VW_CHANNEL_UINT8;

      switch (ctx.output_components)
      {
        case 1:
          fmt.pixel_format = VW_PIXEL_GRAY;
          fmt.planes=1;
          break;
        case 2:
          fmt.pixel_format = VW_PIXEL_GRAYA;
          fmt.planes=1;
          break;
        case 3:
          fmt.pixel_format = VW_PIXEL_RGB;
          fmt.planes=1;
          break;
        case 4:
          fmt.pixel_format = VW_PIXEL_RGBA;
          fmt.planes=1; break;
        default:
          fmt.planes = ctx.output_components;
          fmt.pixel_format = VW_PIXEL_SCALAR;
          break;
      }

      cstride = ctx.output_components;
      line = 0;
    }

    void fini() {
      jpeg_abort_decompress(&ctx);
    }

  public:
    ImageFormat fmt;
    size_t line;

    Data(const uint8* data, size_t len) {
      ctx.err = jpeg_std_error(&jerr);
      jerr.error_exit = &vw_jpeg_error_exit;

      jpeg_create_decompress(&ctx);
      jpeg_mem_src(&ctx, reinterpret_cast<unsigned char*>(const_cast<uint8*>(data)), len);
      init();
    }

    size_t byte_size(size_t pixels) {
      return cstride * pixels;
    }

    void readline(uint8* data) {
      jpeg_read_scanlines(&ctx, &data, 1);
      line++;
    }

    void restart() {
      fini();
      init();
    }

    ~Data() {
      jpeg_destroy_decompress(&ctx);
    }
};

void DstEncodedImageResourceJPEG::set_decoded_data( VarArray<uint8>& data);
void DstEncodedImageResourceJPEG::write( ImageBuffer const& buf, BBox2i const& bbox );
void DstEncodedImageResourceJPEG::flush();
#endif

}
