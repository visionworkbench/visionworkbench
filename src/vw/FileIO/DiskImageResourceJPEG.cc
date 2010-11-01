// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file DiskImageResourceJPEG.cc
///
/// Provides support for the JPEG file format.
///

#ifdef _MSC_VER
#pragma warning(disable:4005)
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <vw/FileIO/DiskImageResourceJPEG.h>
#include <vw/Core/Exception.h>

#include <vector>
#include <boost/scoped_array.hpp>

using namespace vw;

int DiskImageResourceJPEG::default_subsampling_factor = 1;
float DiskImageResourceJPEG::default_quality = 0.95f;

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
    vw_throw( ArgumentErr() << "DiskImageResourceJPEG: Cannot open non-jpeg files.\n" );
  else
    vw_throw( IOErr() << "DiskImageResourceJPEG error: " << buffer );
}

/* A struct to handle the data for jpeglib's internals. Have to use it so
 * we can hide the various types from our own header file, which is not
 * allowed to depend on libjpeg header files.
*/
class DiskImageResourceJPEG::vw_jpeg_decompress_context
{
  // The enclosing class. Need the pointer so we can access its member
  // variables.
  DiskImageResourceJPEG *outer;

  jpeg_error_mgr jerr;

public:
  int current_line;
  jpeg_decompress_struct decompress_ctx;
  JSAMPARRAY scanline;
  int cstride;

  vw_jpeg_decompress_context(DiskImageResourceJPEG *outer) : outer(outer)
  {
    current_line = -1;

    // Set the position of the input stream at zero. Create a new decompress
    // context.
    fseek(outer->m_file_ptr, outer->m_byte_offset, SEEK_SET);

    decompress_ctx.err = jpeg_std_error(&jerr);
    jerr.error_exit = &vw_jpeg_error_exit;

    jpeg_create_decompress(&decompress_ctx);
    jpeg_stdio_src(&decompress_ctx, outer->m_file_ptr);

    // Read the header information and parse it
    jpeg_read_header(&decompress_ctx, TRUE);

    // If the user has requested a grayscale image, we can speed up the
    // decoding process by forcing the jpeg driver to decode the image
    // in grayscale (since the color components need not be processed).
    // This is commented out for now because it seems to be involved
    // in some screwey behavior on Linux.
#if 0
    if (dest.format.pixel_format == VW_PIXEL_GRAY)
      decompress_ctx.out_color_space = JCS_GRAYSCALE;
#endif

    // If the user has requested a subsampled version of the original
    // image, we pass this request along to the JPEG driver here.
    decompress_ctx.scale_num = 1;
    decompress_ctx.scale_denom = outer->m_subsample_factor;

    // Start the decompression
    jpeg_start_decompress(&decompress_ctx);

    // Set the parent class's format.
    outer->m_format.cols = decompress_ctx.output_width;
    outer->m_format.rows = decompress_ctx.output_height;
    outer->m_format.channel_type = VW_CHANNEL_UINT8;

    switch (decompress_ctx.output_components)
    {
      case 1:
        outer->m_format.pixel_format = VW_PIXEL_GRAY;
        outer->m_format.planes=1;
        break;
      case 2:
        outer->m_format.pixel_format = VW_PIXEL_GRAYA;
        outer->m_format.planes=1;
        break;
      case 3:
        outer->m_format.pixel_format = VW_PIXEL_RGB;
        outer->m_format.planes=1;
        break;
      case 4:
        outer->m_format.pixel_format = VW_PIXEL_RGBA;
        outer->m_format.planes=1; break;
      default:
        outer->m_format.planes = decompress_ctx.output_components;
        outer->m_format.pixel_format = VW_PIXEL_SCALAR;
        break;
    }

    // Allocate the scanline
    cstride = decompress_ctx.output_components;
    scanline = (*(decompress_ctx.mem->alloc_sarray))
      ((j_common_ptr) &decompress_ctx, JPOOL_IMAGE, cstride * decompress_ctx.output_width, 1);

    // Finally, we're now at line 0, so set it.
    current_line = 0;
  }

  ~vw_jpeg_decompress_context()
  {
    if(current_line >= 0) {
      jpeg_abort_decompress(&decompress_ctx);
      jpeg_destroy_decompress(&decompress_ctx);
    }
  }


  /* Reads the next line into scanline. Advances current_line appropriately.
   *
   * Warning!
   * We set up error handling here using the struct's error context. This
   * will clobber any error handling you already have, so be sure to
   * reset it!
   */
  void readline() {
    advance(1);
  }

  /* Advances the current point in the file "lines" number of lines.
   * Sets current_line accordingly.
   *
   * Warning!
   * We set up error handling here using the struct's error context. This
   * will clobber any error handling you already have, so be sure to
   * reset it!
  */
  void advance(size_t lines)
  {
    for(size_t i=0; i < lines; i++)
    {
      jpeg_read_scanlines(&decompress_ctx, scanline, 1);
      current_line++;
    }
  }

};

/// Close the JPEG file when the object is destroyed
DiskImageResourceJPEG::~DiskImageResourceJPEG() {
  this->flush();
}

/// Flush the buffered data to disk
void DiskImageResourceJPEG::flush()
{

  if (m_file_ptr) {
    fclose((FILE*)m_file_ptr);
    m_file_ptr = NULL;
  }
}

/// Bind the resource to a file for reading.  Confirm that we can open
/// the file and that it has a sane pixel format.
void DiskImageResourceJPEG::open( std::string const& filename, int subsample_factor, size_t byte_offset )
{
  if (subsample_factor == 1 || subsample_factor == 2 ||
      subsample_factor == 4 || subsample_factor == 8)
  {
    m_subsample_factor = subsample_factor;
  } else {
    vw_throw( ArgumentErr() << "DiskImageResourceJPEG: subsample_factor must be 1, 2, 4, or 8" );
  }

  if(m_file_ptr)
    vw_throw( IOErr() << "DiskImageResourceJPEG: A file is already open." );

  // Open the file on disk
  FILE * infile;
  if ((infile = fopen(filename.c_str(), "rb")) == NULL) {
    vw_throw( ArgumentErr() << "Failed to open \"" << filename << "\" using libJPEG." );
  }
  if (byte_offset) fseek(infile, byte_offset, SEEK_SET);

  m_byte_offset = byte_offset;
  m_filename = filename;
  m_file_ptr = infile;

  // Now that all the needed members are set up, initialize the
  // decompress context.
  ctx = boost::shared_ptr<DiskImageResourceJPEG::vw_jpeg_decompress_context>(new DiskImageResourceJPEG::vw_jpeg_decompress_context(this));
}

/// Bind the resource to a file for writing.
void DiskImageResourceJPEG::create( std::string const& filename,
                                        ImageFormat const& format )
{
  if( format.planes!=1 && format.pixel_format!=VW_PIXEL_SCALAR )
    vw_throw( NoImplErr() << "JPEG doesn't support multi-plane images with compound pixel types." );
  if(m_file_ptr)
    vw_throw( IOErr() << "DiskImageResourceJPEG: A file is already open." );

  // Open the file on disk
  FILE * outfile;
  if ((outfile = fopen(filename.c_str(), "wb")) == NULL) {
    vw_throw( IOErr() << "Failed to open \"" << filename << "\" using libJPEG." );
  }
  m_filename = filename;
  m_format = format;
  m_file_ptr = outfile;

  // The JPEG file format only supports 8-bit channel types, so we
  // force that setting here.
  m_format.channel_type = VW_CHANNEL_UINT8;

  // The JPEG file format only supports color spaces without an alpha
  // channel, so we truncate the alpha channel here and print a
  // warning message.
  if (format.pixel_format == VW_PIXEL_GRAYA) {
    m_format.pixel_format = VW_PIXEL_GRAY;
    vw_out(DebugMessage, "fileio") << "DiskImageResourceJPEG: Warning. alpha channel removed.  ";
  } else if (format.pixel_format == VW_PIXEL_RGBA) {
    m_format.pixel_format = VW_PIXEL_RGB;
    vw_out(DebugMessage, "fileio") << "DiskImageResourceJPEG: Warning. alpha channel removed.  ";
  }
}

/* Read the disk image, or a number of scanlines past the current point,
 * into the given buffer. If you're reading lines, it reads them in
 * starting at the top of the ImageBuffer given. It detects whether to read
 * lines automatically depending on whether bbox represents the whole
 * image or not.
*/
void DiskImageResourceJPEG::read( ImageBuffer const& dest, BBox2i const& bbox) const
{
  VW_ASSERT( int(dest.format.cols)==bbox.width() && int(dest.format.rows)==bbox.height(),
             ArgumentErr() << "DiskImageResourceJPEG (read) Error: Destination buffer has wrong dimensions!" );

  const int start_row = bbox.min().y();
  const unsigned int end_row = bbox.max().y();

  // If we're starting from the beginning, or if we don't have an open
  // decompress context, restart from the beginning of the file.
  if(ctx->current_line != start_row)
  {
    if(ctx->current_line > start_row)
      read_reset();
    ctx->advance(start_row - ctx->current_line);
  }

  // Now read.
  boost::scoped_array<uint8> buf( new uint8[ctx->cstride * bbox.width() * bbox.height()] );

  int32 offset = 0;
  while ( ctx->decompress_ctx.output_scanline < end_row)
  {
    ctx->readline();

    // Copy the data over into a contiguous buffer, which is what
    // convert() expects when we call it below.  This extra copy of
    // the data is undesirable, but probably not a huge performance
    // hit compared to reading the data from disk.
      std::memcpy(buf.get() + offset,
                  ctx->scanline[0] + ctx->cstride * bbox.min().x(),
                  bbox.width() * ctx->cstride);
      offset += bbox.width() * ctx->cstride;

    //for (int i = 0; i < ctx->scanline_size; i++)
    //  buf[ctx->scanline_size * (ctx->decompress_ctx.output_scanline - start_row - 1) + i]
    //    = (uint8)(((JSAMPARRAY)(ctx->scanline))[0][i]);
  }

  // Set up an image buffer around the tdata_t buf object that
  // jpeg used to copy it's data.
  ImageBuffer src;
  src.data = buf.get();
  src.format = m_format;

  // The number of rows we're reading may be different, seeing as how we
  // can read sequential lines.
  src.format.rows = bbox.height();
  src.format.cols = bbox.width();

  // This is part of the grayscale optimization that we remove due to
  // wonkiness on Linux.
#if 0
  if (dest.format.pixel_format == VW_PIXEL_GRAY)
    src.format.pixel_format = VW_PIXEL_GRAY;
#endif

  src.cstride = ctx->cstride;
  src.rstride = src.cstride * src.format.cols;
  src.pstride = src.rstride * src.format.rows;
  convert( dest, src, m_rescale );

}

void DiskImageResourceJPEG::read_reset() const {
  ctx = boost::shared_ptr<DiskImageResourceJPEG::vw_jpeg_decompress_context>(new DiskImageResourceJPEG::vw_jpeg_decompress_context(const_cast<DiskImageResourceJPEG*>(this)));
}

// Write the given buffer into the disk image.
void DiskImageResourceJPEG::write( ImageBuffer const& src, BBox2i const& bbox )
{
  VW_ASSERT( bbox.width()==int(cols()) && bbox.height()==int(rows()),
             NoImplErr() << "DiskImageResourceJPEG does not support partial writes." );
  VW_ASSERT( src.format.cols==uint32(cols()) && src.format.rows==uint32(rows()),
             IOErr() << "Buffer has wrong dimensions in JPEG write." );

  // Set up the JPEG data structures
  jpeg_compress_struct cinfo;
  jpeg_error_mgr jerr;

  cinfo.err = jpeg_std_error(&jerr);
  jerr.error_exit = &vw_jpeg_error_exit;

  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, (FILE*)m_file_ptr);

  cinfo.image_width = m_format.cols;
  cinfo.image_height = m_format.rows;

  switch (m_format.pixel_format)
  {
    case VW_PIXEL_SCALAR:
      cinfo.input_components = m_format.planes;
      cinfo.in_color_space = JCS_UNKNOWN;
      break;
    case VW_PIXEL_GRAY:
      cinfo.input_components = 1;
      cinfo.in_color_space = JCS_GRAYSCALE;
      break;
    case VW_PIXEL_RGB:
      cinfo.input_components = 3;
      cinfo.in_color_space = JCS_RGB;
      break;
    default:
      vw_throw( IOErr() << "DiskImageResourceJPEG: Unsupported pixel type (" << m_format.pixel_format << ")." );
      break;
  }

  // Set up the default values for the header and set the compression
  // quality
  jpeg_set_defaults(&cinfo);
  jpeg_set_quality(&cinfo, (int)(100*m_quality), TRUE); // limit to baseline-JPEG values

  // Set up the image buffer and convert the data into this buffer
  boost::scoped_array<uint8> buf( new uint8[cinfo.image_width*cinfo.input_components*cinfo.image_height] );
  ImageBuffer dst(m_format, buf.get());

  convert( dst, src, m_rescale );

  jpeg_start_compress(&cinfo, TRUE);
  int row_stride = cinfo.image_width*cinfo.input_components;

  // Write the image data to disk.
  JSAMPROW row_pointer[1];
  while (cinfo.next_scanline < cinfo.image_height) {
    row_pointer[0] = &(((uint8*)dst.data)[cinfo.next_scanline * row_stride]);
    jpeg_write_scanlines(&cinfo, row_pointer, 1);
  }

  // Clean up
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);
}

// A FileIO hook to open a file for reading
DiskImageResource* DiskImageResourceJPEG::construct_open( std::string const& filename ) {
  return new DiskImageResourceJPEG( filename );
}

// A FileIO hook to open a file for writing
DiskImageResource* DiskImageResourceJPEG::construct_create( std::string const& filename,
                                                                    ImageFormat const& format ) {
  return new DiskImageResourceJPEG( filename, format );
}
