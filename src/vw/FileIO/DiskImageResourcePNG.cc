// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
//
// Copyright 2006 Carnegie Mellon University. All rights reserved.
//
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
//
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
//
// __END_LICENSE__

/// \file DiskImageResourcePNG.cc
///
/// Provides support for the PNG file format.
///
/// FIXME Currently we do not support 1-bit images, though we could.
///

#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <vector>
#include <fstream>

#include <png.h>

// Non-local error return
// NOTE: To maintain the stack in the presence of
// setjmp/longjump and C++, You cannot create something (with a destructor) on
// the stack and then allow it to longjmp. Therefore, every stack allocation
// needs a setjmp after it. It also makes sense, for sane backtrace purposes,
// to put a setjmp in every function that can longjmp, so the backtrace
// originates from that function.
#include <csetjmp>

#include <vw/Core/Exception.h>
#include <vw/Image/Manipulation.h>
#include <vw/FileIO/DiskImageResourcePNG.h>

using namespace vw;

static const size_t ERROR_MSG_SIZE = 256;

extern "C" {
  struct vw_png_err_mgr {
    jmp_buf error_return;
    char error_msg[ERROR_MSG_SIZE];
  };
}

static void png_error_handler(png_structp png_ptr, png_const_charp error_msg)
{
  vw_png_err_mgr *mgr;
  if (!(mgr = static_cast<vw_png_err_mgr*>(png_get_error_ptr(png_ptr))))
  {
    // We're in big trouble. libpng expects us not to return. Stack could be
    // trashed, and we have nowhere to go. It's not safe to throw here. Bail out.
    vw_out(ErrorMessage, "fileio")
      << "Error while recovering from error in PNG handler: "
      << error_msg << std::endl;
    abort();
  }

  mgr->error_msg[0] = 0; // prep for strncat
  strncat(mgr->error_msg, error_msg, ERROR_MSG_SIZE);

  longjmp(mgr->error_return, 1);
}

/************************************************************************
 ********************** PNG CONTEXT STRUCTURES **************************
 **** These things encapsulate the read/write to the PNG file itself. ***
************************************************************************/
// Common stuff for both the read and write contexts.
struct DiskImageResourcePNG::vw_png_context {
  // The PNG comments
  std::vector<DiskImageResourcePNG::Comment> comments;

  vw_png_context(DiskImageResourcePNG *outer) {
    this->outer = outer;
  }

  virtual ~vw_png_context() { }

protected:
  // Pointer to the containing class, necessary to access some of its
  // members.
  DiskImageResourcePNG *outer;

  // Structures from libpng.
  png_structp png_ptr;
  png_infop info_ptr;

  // Error manager to prevent libpng from calling abort()
  mutable vw_png_err_mgr err_mgr;

  // Pointer to actual file on disk. shared_ptr instead of std::auto_ptr
  // because of problems with auto_ptr's deletion.
  boost::shared_ptr<std::fstream> m_file;
};

// Context for reading.
struct DiskImageResourcePNG::vw_png_read_context:
  public DiskImageResourcePNG::vw_png_context
{

  // Current line we're at in the image.
  int current_line;
  boost::shared_array<uint8> scanline;
  int scanline_size;

  // Other image data.
  int bytes_per_channel;
  int channels;
  bool interlaced;

  // Open a PNG context from a file, for reading.
  vw_png_read_context(DiskImageResourcePNG *outer):
    vw_png_context(outer)
  {
    m_file = boost::shared_ptr<std::fstream>( new std::fstream( outer->m_filename.c_str(), std::ios_base::in | std::ios_base::binary ) );
    if(!m_file)
      vw_throw(IOErr() << "DiskImageResourcePNG: Unable to open input file " << outer->m_filename << ".");

    // Read the first 8 bytes to make sure it's actually a PNG.
    char sig[8];
    m_file->read(sig, 8);
    bool is_png = !png_sig_cmp((png_byte*)sig, 0, 8);
    if(!is_png)
      vw_throw(IOErr() << "DiskImageResourcePNG: Input file " << outer->m_filename << " is not a valid PNG file.");

    // Allocate the structures.
    png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,
                                     &err_mgr, png_error_handler, NULL);
    if(!png_ptr)
      vw_throw(IOErr() << "DiskImageResourcePNG: Failure to create read structure for file " << outer->m_filename << ".");

    if (setjmp(err_mgr.error_return))
      vw_throw( vw::IOErr() << "DiskImageResourcePNG: A libpng error occurred. " << err_mgr.error_msg );

    info_ptr = png_create_info_struct(png_ptr);
    if(!info_ptr) {
      png_destroy_read_struct(&png_ptr, NULL, NULL);
      vw_throw(IOErr() << "DiskImageResourcePNG: Failure to create info structure for file " << outer->m_filename << ".");
    }

    endinfo_ptr = png_create_info_struct(png_ptr);
    if(!endinfo_ptr) {
      png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
      vw_throw(IOErr() << "DiskImageResourcePNG: Failure to create info structure for file " << outer->m_filename << ".");
    }

    // Must call this as we're using fstream and not FILE*
    png_set_read_fn(png_ptr, reinterpret_cast<voidp>(m_file.get()), read_data);

    // Rewind to the beginning of the file.
    png_set_sig_bytes(png_ptr, 8);

    // Read in the info pointer (some stuff will get changed, and we run
    // png_read_update_info).
    png_read_info(png_ptr, info_ptr);

    // Fetch the comments.
    read_comments();

    // Set up expansion. Palette images get expanded to RGB, grayscale
    // images of less than 8 bits per channel are expanded to 8 bits per
    // channel, and tRNS chunks are expanded to alpha channels.
    png_set_expand(png_ptr);

    // png_uint_32 is usually a long, which could be 4 or 8 bytes,
    // depending on platform. be careful.
    png_uint_32 cols, rows;

    int bit_depth;
    int color_type;
    int interlace_type;
    int channels;
    int filter_method;
    int compression_type;
//    int num_passes; // For interlacing :-(

    // Read up to the start of the data, and set some values.
    png_get_IHDR(png_ptr, info_ptr, &cols, &rows, &bit_depth, &color_type,
                 &interlace_type, &compression_type, &filter_method);

    switch(bit_depth) {
      case 1:
      case 2:
      case 4:
      case 8:
        bytes_per_channel = 1;
        outer->m_format.channel_type = VW_CHANNEL_UINT8;
        break;
      case 16:
        bytes_per_channel = 2;
        outer->m_format.channel_type = VW_CHANNEL_UINT16;
        break;
      default:
        // Unreachable
        bytes_per_channel = 1;
        outer->m_format.channel_type = VW_CHANNEL_UINT8;
        break;
    }
    switch(color_type) {
      case PNG_COLOR_TYPE_GRAY:
        channels = 1;
        outer->m_format.pixel_format = VW_PIXEL_GRAY;
        break;
      case PNG_COLOR_TYPE_GRAY_ALPHA:
        channels = 2;
        outer->m_format.pixel_format = VW_PIXEL_GRAYA;
        break;
      case PNG_COLOR_TYPE_PALETTE:
        channels = 3;
        outer->m_format.pixel_format = VW_PIXEL_RGB;
        if( png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS) ) {
          channels++;
          outer->m_format.pixel_format = VW_PIXEL_RGBA;
        }
        break;
      case PNG_COLOR_TYPE_RGB:
        channels = 3;
        outer->m_format.pixel_format = VW_PIXEL_RGB;
        break;
      case PNG_COLOR_TYPE_RGB_ALPHA:
        channels = 4;
        outer->m_format.pixel_format = VW_PIXEL_RGBA;
        break;
      default:
        // Unreachable
        channels = 4;
        outer->m_format.pixel_format = VW_PIXEL_RGBA;
        break;
    }
    if(interlace_type != PNG_INTERLACE_NONE)
      interlaced = true;
    else
      interlaced = false;

    png_read_update_info(png_ptr, info_ptr);
    png_get_IHDR(png_ptr, info_ptr, &cols, &rows, &bit_depth, &color_type, &interlace_type, &compression_type, &filter_method);

    outer->m_format.cols = cols;
    outer->m_format.rows = rows;
    outer->m_format.planes = 1;

    // Allocate the scanline.
    scanline_size = cols * bytes_per_channel * channels;
    scanline = boost::shared_array<uint8>(new uint8[scanline_size]);

    png_start_read_image(png_ptr);

    current_line = 0;
  }

  virtual ~vw_png_read_context() {
    if (setjmp(err_mgr.error_return))
      vw_throw( vw::IOErr() << "DiskImageResourcePNG: A libpng error occurred. " << err_mgr.error_msg );
    png_destroy_read_struct(&png_ptr, &info_ptr, &endinfo_ptr);
    m_file->close();
  }

  void readline() {
    if (setjmp(err_mgr.error_return))
      vw_throw( vw::IOErr() << "DiskImageResourcePNG: A libpng error occurred. " << err_mgr.error_msg );
    png_read_row(png_ptr, static_cast<png_bytep>(scanline.get()), NULL);
    current_line++;
  }

  void readall(boost::scoped_array<uint8> &dst) {
    if(current_line != 0)
      vw_throw(IOErr() << "DiskImageResourcePNG: cannot read entire file unless line marker set at beginning.");
    png_bytep row_pointers[outer->m_format.rows];
    if (setjmp(err_mgr.error_return))
      vw_throw( vw::IOErr() << "DiskImageResourcePNG: A libpng error occurred. " << err_mgr.error_msg );
    for(int i=0; i < outer->m_format.rows; i++)
      row_pointers[i] = static_cast<png_bytep>(dst.get()) + i*scanline_size;
    png_read_image(png_ptr, row_pointers);
    current_line = outer->m_format.rows;
  }


  // Advances place in the image by line lines.
  void advance(size_t lines) {
    if (setjmp(err_mgr.error_return))
      vw_throw( vw::IOErr() << "DiskImageResourcePNG: A libpng error occurred. " << err_mgr.error_msg );
    for(size_t i = 0; i < lines; i++) {
      png_read_row(png_ptr, NULL, NULL);
      current_line++;
    }
  }

private:
  // Structures from libpng.
  png_infop endinfo_ptr;

  // Function for reading data, given to PNG.
  static void read_data( png_structp png_ptr, png_bytep data, png_size_t length ) {
    std::fstream *fs = static_cast<std::fstream*>(png_get_io_ptr(png_ptr));;
    fs->read( reinterpret_cast<char*>(data), length );
  }

  // Fetches the comments out of the PNG when we first open it.
  void read_comments() {
    png_text *text_ptr;
    int num_comments = png_get_text(png_ptr, info_ptr, &text_ptr, 0);
    comments.clear();
    for ( int i=0; i<num_comments; ++i ) {
      DiskImageResourcePNG::Comment c;
      c.key = text_ptr[i].key;
      c.text = text_ptr[i].text;
#ifdef PNG_iTXt_SUPPORTED
      c.lang = text_ptr[i].lang;
      c.lang_key = text_ptr[i].lang_key;
#endif
      switch( text_ptr[i].compression ) {
        case PNG_TEXT_COMPRESSION_NONE:
          c.utf8 = false;
          c.compressed = false;
          break;
        case PNG_TEXT_COMPRESSION_zTXt:
          c.utf8 = false;
          c.compressed = true;
          break;
#ifdef PNG_iTXt_SUPPORTED
        case PNG_ITXT_COMPRESSION_NONE:
          c.utf8 = true;
          c.compressed = false;
          break;
        case PNG_ITXT_COMPRESSION_zTXt:
          c.utf8 = true;
          c.compressed = true;
          break;
#endif
        default:
          vw_out(WarningMessage, "fileio") << "Unsupported PNG comment type in PNG read!" << std::endl;
          continue;
      }
      comments.push_back( c );
    }
  }
};

// Context for writing
struct DiskImageResourcePNG::vw_png_write_context:
  public DiskImageResourcePNG::vw_png_context
{
  int scanline_size;

  vw_png_write_context(DiskImageResourcePNG *outer, const DiskImageResourcePNG::Options &options):
    vw_png_context(outer)
  {
    m_file = boost::shared_ptr<std::fstream>( new std::fstream( outer->m_filename.c_str(), std::ios_base::out | std::ios_base::binary ) );
    if(!m_file)
      vw_throw(IOErr() << "DiskImageResourcePNG: Unable to open output file " << outer->m_filename << ".");

    // Allocate the structures.
    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                                     &err_mgr, png_error_handler, NULL);
    if(!png_ptr)
      vw_throw(IOErr() << "DiskImageResourcePNG: Failure to create read structure for file " << outer->m_filename << ".");

    if (setjmp(err_mgr.error_return))
      vw_throw( vw::IOErr() << "DiskImageResourcePNG: A libpng error occurred. " << err_mgr.error_msg );

    info_ptr = png_create_info_struct(png_ptr);
    if(!info_ptr) {
      png_destroy_read_struct(&png_ptr, NULL, NULL);
      vw_throw(IOErr() << "DiskImageResourcePNG: Failure to create info structure for file " << outer->m_filename << ".");
    }

    // Must call this as we're using fstream and not FILE*
    png_set_write_fn(png_ptr, reinterpret_cast<voidp>(m_file.get()), write_data, flush_data);

    // Set some needed values.
    int width     = outer->m_format.cols;
    int height    = outer->m_format.rows;
    int channels  = num_channels(outer->m_format.pixel_format);
    int bit_depth = outer->m_format.channel_type == VW_CHANNEL_UINT16 ? 16 : 8;

    int color_type;
    switch(outer->m_format.pixel_format) {
      case VW_PIXEL_GRAY:   color_type = PNG_COLOR_TYPE_GRAY;       break;
      case VW_PIXEL_GRAYA:  color_type = PNG_COLOR_TYPE_GRAY_ALPHA; break;
      case VW_PIXEL_RGB:    color_type = PNG_COLOR_TYPE_RGB;        break;
      case VW_PIXEL_RGBA:   color_type = PNG_COLOR_TYPE_RGBA;       break;
      default:              color_type = PNG_COLOR_TYPE_RGBA;       break;
    }
    if(options.using_palette) {
      color_type = PNG_COLOR_TYPE_PALETTE;
      channels = 3;
    }

    int interlace_type;
    if(options.using_interlace)
      interlace_type = PNG_INTERLACE_ADAM7;
    else
      interlace_type = PNG_INTERLACE_NONE;

    int compression_type = PNG_COMPRESSION_TYPE_DEFAULT;
    int filter_method    = PNG_FILTER_TYPE_DEFAULT;

    png_set_IHDR(png_ptr, info_ptr, width, height, bit_depth, color_type, interlace_type, compression_type, filter_method);

    if(options.using_palette && options.using_palette_indices) {
      png_colorp palette = reinterpret_cast<png_colorp>(png_malloc( png_ptr, options.palette.cols() * sizeof(png_color) ));
      png_bytep alpha    = reinterpret_cast<png_bytep>(png_malloc( png_ptr, options.palette.cols() * sizeof(png_byte) ));
      if (setjmp(err_mgr.error_return))
        vw_throw( vw::IOErr() << "DiskImageResourcePNG: A libpng error occurred. " << err_mgr.error_msg );
      for ( int i=0; i < int(options.palette.cols()); ++i ) {
        palette[i].red   = options.palette(i,0).r();
        palette[i].green = options.palette(i,0).g();
        palette[i].blue  = options.palette(i,0).b();
        alpha[i]         = options.palette(i,0).a();
      }
      png_set_PLTE( png_ptr, info_ptr, palette, options.palette.cols() );
      if(options.using_palette_alpha) {
        png_set_tRNS( png_ptr, info_ptr, alpha, options.palette.cols(), 0 );
        channels++;
      }
    }

    png_set_compression_level(png_ptr, options.compression_level);

    // Set up the scanline for writing.
    scanline_size = (bit_depth / 8) * width * channels;

    // Finally, write the info out.
    png_write_info(png_ptr, info_ptr);
  }

  virtual ~vw_png_write_context() {
    if (setjmp(err_mgr.error_return))
      vw_throw( vw::IOErr() << "DiskImageResourcePNG: A libpng error occurred. " << err_mgr.error_msg );
    png_write_end(png_ptr, info_ptr);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    m_file->close();
  }

  // Writes the given ImageBuffer (with the same dimensions as m_format)
  // to the file. Closing happens when the context is destroyed.
  void write(const ImageBuffer &buf) const {
    png_bytep row_pointers[outer->m_format.rows];
    if (setjmp(err_mgr.error_return))
      vw_throw( vw::IOErr() << "DiskImageResourcePNG: A libpng error occurred. " << err_mgr.error_msg );
    for(int i=0; i < outer->m_format.rows; i++)
      row_pointers[i] = reinterpret_cast<uint8*>(buf.data) + i * scanline_size;

    png_write_image(png_ptr, row_pointers);
  }

private:
  // Function for libpng to use to write as we're not using FILE*.
  static void write_data( png_structp png_ptr, png_bytep data, png_size_t length ) {
    std::fstream *fs = static_cast<std::fstream*>(png_get_io_ptr(png_ptr));
    fs->write( (char*)data, length );
  }

  static void flush_data( png_structp png_ptr) {
    std::fstream *fs = static_cast<std::fstream*>(png_get_io_ptr(png_ptr));
    fs->flush();
  }
};

// *********************************************************************
// Public interface function definitions
// *********************************************************************

/***********************************************************************
 *********************** CONSTRUCTORS & DESTRUCTORS ********************
***********************************************************************/

DiskImageResourcePNG::DiskImageResourcePNG( std::string const& filename ): DiskImageResource(filename)
{
  open(filename);
}

DiskImageResource* DiskImageResourcePNG::construct_open( std::string const& filename ) {
  return new DiskImageResourcePNG( filename );
}

vw::DiskImageResourcePNG::~DiskImageResourcePNG() {
}

DiskImageResourcePNG::DiskImageResourcePNG( std::string const& filename, ImageFormat const& format ):
  DiskImageResource( filename )
{
  create( filename, format );
}

DiskImageResourcePNG::DiskImageResourcePNG( std::string const& filename, ImageFormat const& format, DiskImageResourcePNG::Options const& options):
  DiskImageResource(filename)
{
  create( filename, format, options );
}

DiskImageResource* DiskImageResourcePNG::construct_create( std::string const& filename,
                                                                   ImageFormat const& format ) {
  return new DiskImageResourcePNG( filename, format );
}

/***********************************************************************
 ************************ STUFF FOR READING ****************************
***********************************************************************/

void DiskImageResourcePNG::open( std::string const& filename ) {
  m_ctx = boost::shared_ptr<vw_png_context>( new vw_png_read_context( const_cast<DiskImageResourcePNG *>(this) ) );
}

void DiskImageResourcePNG::read( ImageBuffer const& dest, BBox2i const& bbox ) const {
  vw_png_read_context *ctx = dynamic_cast<vw_png_read_context *>(m_ctx.get());
  const int start_line = bbox.min().y();
  const int end_line = bbox.max().y();
  // Error checking.
  VW_ASSERT( dest.format.cols == cols(), IOErr() << "DiskImageResourcePNG: Destination buffer has different columns than file.");
  VW_ASSERT( dest.format.rows <= rows(), IOErr() << "DiskImageResourcePNG: Destination buffer has too many rows.");
  VW_ASSERT(bbox.height() == dest.format.rows, ArgumentErr() << "DiskImageResourcePNG: Destination buffer and requested bounding box have different height.");

  boost::scoped_array<uint8> buf( new uint8[ctx->scanline_size * bbox.height()] );
  // Interlacing is causing problems when read line-by-line...I think it's
  // a bug in libpng.
  if( ctx->interlaced ) {
    if( bbox.height() != rows() )
      vw_throw( NoImplErr() << "DiskImageResourcePNG: Reading interlaced files line-by-line is currently unsupported." );

    ctx->readall(buf);
  } else {
    // FIXME: Normal operation. Make this what happens all the time when
    // the libpng bug gets fixed. The bug in question causes the final
    // parts of the expanded, interlaced image (the ones we care about) to
    // seemingly 'lose' a row.

    // If our start line is a spot before the current line, we need to reopen the file.
    if(start_line < ctx->current_line)
      read_reset();
    if(start_line > ctx->current_line)
      ctx->advance(start_line - ctx->current_line);

    // Now read.
    while(ctx->current_line < end_line) {
      ctx->readline();

      // Copy the data over into a contiguous buffer, which is what
      // convert() expects when we call it below. This extra copy of the
      // data is undesirable, but probably not a huge performance hit
      // compared to reading the data from disk.
      for(int i=0; i < ctx->scanline_size; i++)
        buf[ctx->scanline_size * (ctx->current_line - start_line - 1) + i] =
          ctx->scanline[i];
    }

  }

  // Set up an image buffer for convert().
  ImageBuffer src;
  src.data = buf.get();
  src.format = m_format;
  // The number of rows we're reading may be different, seeing as how we
  // can read sequential lines.
  src.format.rows = bbox.height();
  src.cstride = ctx->scanline_size / m_format.cols;
  src.rstride = ctx->scanline_size;
  src.pstride = ctx->scanline_size * bbox.height();
  convert(dest, src);
}

void DiskImageResourcePNG::read_reset() const {
  m_ctx = boost::shared_ptr<DiskImageResourcePNG::vw_png_read_context>( new DiskImageResourcePNG::vw_png_read_context( const_cast<DiskImageResourcePNG *>(this) ) );
}

unsigned DiskImageResourcePNG::num_comments() const {
  return m_ctx->comments.size();
}

DiskImageResourcePNG::Comment const& DiskImageResourcePNG::get_comment( unsigned i ) const {
  return m_ctx->comments[i];
}

std::string const& DiskImageResourcePNG::get_comment_key( unsigned i ) const {
  return get_comment(i).key;
}

std::string const& DiskImageResourcePNG::get_comment_value( unsigned i ) const {
  return get_comment(i).text;
}

/***********************************************************************
 ************************ STUFF FOR WRITING ****************************
***********************************************************************/

void DiskImageResourcePNG::create( std::string const& filename, ImageFormat const& format ) {
  DiskImageResourcePNG::create( filename, format, DiskImageResourcePNG::Options() );
}

void DiskImageResourcePNG::create( std::string const& filename, ImageFormat const& format, DiskImageResourcePNG::Options const& options ) {
  if(m_ctx.get() != NULL)
    vw_throw( IOErr() << "DiskImageResourcePNG: A file is already open." );

  m_filename = filename;
  m_format = format;

  m_ctx = boost::shared_ptr<vw_png_context>( new vw_png_write_context( const_cast<DiskImageResourcePNG *>(this), options ) );
}

void DiskImageResourcePNG::write( ImageBuffer const& src, BBox2i const& bbox ) {
  vw_png_write_context *ctx = dynamic_cast<vw_png_write_context *>( m_ctx.get() );

  VW_ASSERT( bbox.width()==int(cols()) && bbox.height()==int(rows()),
             NoImplErr() << "DiskImageResourcePNG does not support partial writes." );
  VW_ASSERT( src.format.cols==cols() && src.format.rows==rows(),
             ArgumentErr() << "DiskImageResourcePNG: Buffer has wrong dimensions in PNG write." );

  // Set up the image buffer and convert the data into this buffer.
  ImageBuffer dst;
  boost::scoped_array<uint8> buf(new uint8[ctx->scanline_size * bbox.height()]);

  dst.data = buf.get();
  dst.format = m_format;
  dst.format.rows = bbox.height();

  if (dst.format.channel_type != VW_CHANNEL_UINT16)
    dst.format.channel_type = VW_CHANNEL_UINT8;

  dst.cstride = num_channels(dst.format.pixel_format) * channel_size(dst.format.channel_type);
  dst.rstride = dst.cstride * dst.format.cols;
  dst.pstride = dst.rstride * dst.format.rows;

  convert(dst, src);

  // Write.
  ctx->write(dst);
}
