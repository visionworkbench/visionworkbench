// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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

#include <vw/FileIO/DiskImageResourcePNG.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Exception.h>
#include <vw/Image/Manipulation.h>

#include <png.h>
#include <vector>
#include <fstream>

using namespace vw;

static void png_error_handler(png_structp /*png_ptr*/, png_const_charp error_msg)
{
  vw_throw(IOErr() << "DiskImageResourcePNG: " << error_msg);
}

// Default PNG file creation options, and related functions.
namespace {
  static int default_compression_level = Z_BEST_COMPRESSION;
}

DiskImageResourcePNG::Options::Options() {
  compression_level = default_compression_level;
  using_interlace = false;
  using_palette = false;
  using_palette_indices = false;
  using_palette_alpha = false;
}

void DiskImageResourcePNG::set_default_compression_level( int level ) {
  default_compression_level = level;
}


/************************************************************************
 ********************** PNG CONTEXT STRUCTURES **************************
 **** These things encapsulate the read/write to the PNG file itself. ***
************************************************************************/

// noncopyable because destructor frees the png_structp
struct png_context_t : private boost::noncopyable {
  png_structp ptr;
  png_infop info, info_end;
  boost::shared_ptr<std::fstream> file;

  enum Mode {
    PNG_UNINIT,
    PNG_READ,
    PNG_WRITE
  };

  png_context_t()
    : ptr(0), info(0), info_end(0), m_mode(PNG_UNINIT) {}

  explicit png_context_t(const char* filename, Mode m)
    : ptr(0), info(0), info_end(0), m_mode(m)
  {
    VW_ASSERT(filename, ArgumentErr() << "Filename cannot be null");
    VW_ASSERT(m != PNG_UNINIT, ArgumentErr() << "png_context_t constructed with uninitialized argument");

    file.reset(new std::fstream(filename, (m == PNG_READ ? std::ios::in : std::ios::out) | std::ios_base::binary ));

    if(!file || !file->is_open())
      vw_throw(ArgumentErr() << "DiskImageResourcePNG: Unable to open file " << filename << ".");

    ptr = create();
    if(!ptr)
      vw_throw(IOErr() << "DiskImageResourcePNG: Failed to create context struct for " << (m_mode == PNG_READ ? "read." : "write."));

    info = png_create_info_struct(ptr);
    if(!info) {
      destroy(&ptr, 0);
      vw_throw(IOErr() << "DiskImageResourcePNG: Failed to create info struct for " << (m_mode == PNG_READ ? "read." : "write."));
    }

    if (m == PNG_READ) {
      info_end = png_create_info_struct(ptr);
      if(!info_end) {
        destroy(&ptr, &info, 0);
        vw_throw(IOErr() << "DiskImageResourcePNG: Failed to create end info struct for " << (m_mode == PNG_READ ? "read." : "write."));
      }
    }
  }

  ~png_context_t() VW_NOTHROW {
    if (m_mode == PNG_UNINIT)
      return;
    destroy(&ptr, &info, &info_end);
    if (file->is_open())
      file->close();
  }
  private:
  Mode m_mode;

  png_structp create() const {
    if (m_mode == PNG_READ)
      return png_create_read_struct(PNG_LIBPNG_VER_STRING, 0, png_error_handler, NULL);
    else
      return png_create_write_struct(PNG_LIBPNG_VER_STRING, 0, png_error_handler, NULL);
  }

  void destroy(png_structpp p, png_infopp info, png_infopp end_info = 0) const {
    if (m_mode == PNG_READ)
      png_destroy_read_struct(p, info, end_info);
    else
      png_destroy_write_struct(p, info);
  }
};

// Common stuff for both the read and write contexts.
struct DiskImageResourcePNG::vw_png_context
{
  vw_png_context(DiskImageResourcePNG *outer)
    : outer(outer) { };

  virtual ~vw_png_context() { }

  std::vector<DiskImageResourcePNG::Comment> comments;  // The PNG comments

  virtual void read_comments() {}

  // Cache the column stride for read/writes in bounding boxes
  int32 cstride;

protected:
  // Pointer to containing class, necessary to access some of its members.
  DiskImageResourcePNG *outer;
};

// Context for reading.
struct DiskImageResourcePNG::vw_png_read_context:
  public DiskImageResourcePNG::vw_png_context
{
  png_context_t ctx;

  // Current line we're at in the image.
  int current_line;
  boost::shared_array<uint8> scanline;
  bool comments_loaded;

  // Other image data.
  int bytes_per_channel;
  int channels;
  bool interlaced;

  // Open a PNG context from a file, for reading.
  vw_png_read_context(DiskImageResourcePNG *outer) :
    vw_png_context(outer), ctx(outer->m_filename.c_str(), png_context_t::PNG_READ), current_line(0), comments_loaded(false)
  {
    // Read the first 8 bytes to make sure it's actually a PNG.
    char sig[8];
    ctx.file->read(sig, 8);
    if(png_sig_cmp((png_byte*)sig, 0, 8))
      vw_throw(ArgumentErr() << "DiskImageResourcePNG: Input file "
               << outer->m_filename << " is not a valid PNG file.");

    // Must call this as we're using fstream and not FILE*
    png_set_read_fn(ctx.ptr, reinterpret_cast<voidp>(ctx.file.get()), read_data);

    // Rewind to the beginning of the file.
    png_set_sig_bytes(ctx.ptr, 8);

    // Read in the info pointer (some stuff will get changed, and we run
    // png_read_update_info).
    png_read_info(ctx.ptr, ctx.info);

    // Set up expansion. Palette images get expanded to RGB, grayscale
    // images of less than 8 bits per channel are expanded to 8 bits per
    // channel, and tRNS chunks are expanded to alpha channels.
    png_set_expand(ctx.ptr);

    // png_uint_32 is usually a long, which could be 4 or 8 bytes,
    // depending on platform. be careful.
    png_uint_32 cols, rows;

    int bit_depth;
    int color_type;
    int interlace_type;
    int channels = 1;  // Set a default value to avoid compiler warnings.
    int filter_method;
    int compression_type;
//    int num_passes; // For interlacing :-(

    // Read up to the start of the data, and set some values.
    png_get_IHDR(ctx.ptr, ctx.info, &cols, &rows, &bit_depth, &color_type,
                 &interlace_type, &compression_type, &filter_method);

    switch(bit_depth)
    {
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

        // If this is a little endian machine, we need to instruct PNG to
        // swap the endianness as it reads the file.
        // (Note: this call must come AFTER png_read_info)
#       if VW_BYTE_ORDER == VW_LITTLE_ENDIAN
          png_set_swap(ctx.ptr);
#       endif

        break;
      default:
        vw_throw(vw::ArgumentErr() << "Unknown bit depth " << bit_depth);
        break;
    }

    switch(color_type)
    {
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
        if( png_get_valid(ctx.ptr, ctx.info, PNG_INFO_tRNS) ) {
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
        vw_throw(vw::ArgumentErr() << "Unknown color type in png");
    }
    if(interlace_type != PNG_INTERLACE_NONE)
      interlaced = true;
    else
      interlaced = false;

    png_read_update_info(ctx.ptr, ctx.info);
    png_get_IHDR(ctx.ptr, ctx.info, &cols, &rows, &bit_depth, &color_type, &interlace_type, &compression_type, &filter_method);

    outer->m_format.cols = cols;
    outer->m_format.rows = rows;
    outer->m_format.planes = 1;

    // Allocate the scanline.
    cstride = bytes_per_channel * channels;
    scanline = boost::shared_array<uint8>(new uint8[cstride * cols]);

    png_start_read_image(ctx.ptr);
  }

  void readline()
  {
    png_read_row(ctx.ptr, static_cast<png_bytep>(scanline.get()), NULL);
    current_line++;
  }

  void readall(boost::scoped_array<uint8> &dst)
  {
    if(current_line != 0)
      vw_throw(IOErr() << "DiskImageResourcePNG: cannot read entire file unless line marker set at beginning.");

    boost::scoped_array<png_bytep> row_pointers( new png_bytep[outer->m_format.rows] );
    for(size_t i=0; i < outer->m_format.rows; i++)
      row_pointers[i] = static_cast<png_bytep>(dst.get()) + i * outer->m_format.cols * cstride;
    png_read_image(ctx.ptr, row_pointers.get());
    current_line = outer->m_format.rows;
  }


  // Advances place in the image by line lines.
  void advance(size_t lines)
  {
    for(size_t i = 0; i < lines; i++) {
      png_read_row(ctx.ptr, NULL, NULL);
      current_line++;
    }
  }

  // Fetches the comments out of the PNG when we first open it.
  void read_comments()
  {
    if( comments_loaded ) return;
    advance(outer->rows()-current_line);
    png_read_end(ctx.ptr, ctx.info_end);
    comments_loaded = true;

    png_text *text_ptr;
    int num_comments = png_get_text(ctx.ptr, ctx.info_end, &text_ptr, 0);
    comments.clear();
    for ( int i=0; i<num_comments; ++i )
    {
      DiskImageResourcePNG::Comment c;
      c.key = text_ptr[i].key;
      c.text = text_ptr[i].text;
#ifdef PNG_iTXt_SUPPORTED
      c.lang = text_ptr[i].lang;
      c.lang_key = text_ptr[i].lang_key;
#endif
      switch( text_ptr[i].compression )
      {
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

private:
  // Function for reading data, given to PNG.
  static void read_data( png_structp png_ptr, png_bytep data, png_size_t length )
  {
    std::fstream *fs = static_cast<std::fstream*>(png_get_io_ptr(png_ptr));;
    fs->read( reinterpret_cast<char*>(data), length );
  }
};

// Context for writing
struct DiskImageResourcePNG::vw_png_write_context:
  public DiskImageResourcePNG::vw_png_context
{
  png_context_t ctx;

  vw_png_write_context(DiskImageResourcePNG *outer, const DiskImageResourcePNG::Options &options):
    vw_png_context(outer), ctx(outer->m_filename.c_str(), png_context_t::PNG_WRITE)
  {
    // Set some needed values.
    int width     = outer->m_format.cols;
    int height    = outer->m_format.rows;
    int channels  = num_channels(outer->m_format.pixel_format);
    int bit_depth;

    png_set_compression_level(ctx.ptr, Z_BEST_SPEED);

    // Must call this as we're using fstream and not FILE*
    png_set_write_fn(ctx.ptr, reinterpret_cast<voidp>(ctx.file.get()), write_data, flush_data);

    // anything else will be converted to UINT8
    switch (outer->m_format.channel_type) {
      case VW_CHANNEL_INT16:
      case VW_CHANNEL_UINT16:
      case VW_CHANNEL_FLOAT16:
      case VW_CHANNEL_GENERIC_2_BYTE:
        bit_depth = 16;
        break;
      default:
        bit_depth = 8;
        break;
    }

    int color_type = PNG_COLOR_TYPE_GRAY;  // Set a default value to avoid compiler warnings
    switch(outer->m_format.pixel_format) {
      case VW_PIXEL_SCALAR: // fall through
      case VW_PIXEL_GRAY:   color_type = PNG_COLOR_TYPE_GRAY;       break;
      case VW_PIXEL_GRAYA:  color_type = PNG_COLOR_TYPE_GRAY_ALPHA; break;
      case VW_PIXEL_RGB:    color_type = PNG_COLOR_TYPE_RGB;        break;
      case VW_PIXEL_RGBA:   color_type = PNG_COLOR_TYPE_RGBA;       break;
      default: vw_throw(vw::ArgumentErr() << "Unsupported pixel format for png: "
                                          << outer->m_format.pixel_format);
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

    png_set_IHDR(ctx.ptr, ctx.info, width, height, bit_depth, color_type, interlace_type, compression_type, filter_method);

    if(options.using_palette && options.using_palette_indices)
    {
      png_colorp palette = reinterpret_cast<png_colorp>(png_malloc( ctx.ptr, options.palette.cols() * sizeof(png_color) ));
      png_bytep alpha    = reinterpret_cast<png_bytep>(png_malloc( ctx.ptr, options.palette.cols() * sizeof(png_byte) ));

      for ( int i=0; i < int(options.palette.cols()); ++i )
      {
        palette[i].red   = options.palette(i,0).r();
        palette[i].green = options.palette(i,0).g();
        palette[i].blue  = options.palette(i,0).b();
        alpha[i]         = options.palette(i,0).a();
      }
      png_set_PLTE( ctx.ptr, ctx.info, palette, options.palette.cols() );
      if(options.using_palette_alpha) {
        png_set_tRNS( ctx.ptr, ctx.info, alpha, options.palette.cols(), 0 );
        channels++;
      }
    }

    png_set_compression_level(ctx.ptr, options.compression_level);

    // Set up the scanline for writing.
    cstride = (bit_depth / 8) * channels;

    // Finally, write the info out.
    png_write_info(ctx.ptr, ctx.info);

    // If this is a little endian machine, we need to instruct PNG to
    // swap the endianness as it writes the file.
    // libpng checks bit_depth in the set_swap function, so don't bother
    // to check it here. (Note: this call must come AFTER png_write_info)
#   if VW_BYTE_ORDER == VW_LITTLE_ENDIAN
      png_set_swap(ctx.ptr);
#   endif

  }

  // Writes the given ImageBuffer (with the same dimensions as m_format)
  // to the file. Closing happens when the context is destroyed.
  void write(const ImageBuffer &buf) const
  {
    boost::scoped_array<png_bytep> row_pointers( new png_bytep[outer->m_format.rows] );

    for(size_t i=0; i < outer->m_format.rows; i++)
      row_pointers[i] = reinterpret_cast<uint8*>(buf.data) + i * cstride * outer->m_format.cols;

    png_write_image(ctx.ptr, row_pointers.get());
    png_write_end(ctx.ptr, ctx.info);
  }

private:
  // Function for libpng to use to write as we're not using FILE*.
  static void write_data( png_structp png_ptr, png_bytep data, png_size_t length )
  {
    std::fstream *fs = static_cast<std::fstream*>(png_get_io_ptr(png_ptr));
    fs->write( (char*)data, length );
  }

  static void flush_data( png_structp png_ptr)
  {
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

DiskImageResourcePNG::DiskImageResourcePNG( std::string const& filename )
  : DiskImageResource(filename)
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
                                                           ImageFormat const& format )
{
  return new DiskImageResourcePNG( filename, format );
}

/***********************************************************************
 ************************ STUFF FOR READING ****************************
***********************************************************************/

void DiskImageResourcePNG::open( std::string const& /*filename*/ ) {
  m_ctx = boost::shared_ptr<vw_png_context>( new vw_png_read_context( const_cast<DiskImageResourcePNG *>(this) ) );

  // Block reading is supported, we only use it in the event of really large images.
  if ( size_t(cols()*rows()*4*3) > vw_settings().system_cache_size() )
    m_block_size = Vector2i( cols(), 128 ); // 128 seems like a good number.
  else
    m_block_size = Vector2i( cols(), rows() );
}

void DiskImageResourcePNG::read( ImageBuffer const& dest, BBox2i const& bbox ) const
{
  vw_png_read_context *ctx = dynamic_cast<vw_png_read_context *>(m_ctx.get());
  const int start_line = bbox.min().y();
  const int end_line = bbox.max().y();
  // Error checking.
  VW_ASSERT( int(dest.format.cols)==bbox.width() && int(dest.format.rows)==bbox.height(),
             ArgumentErr() << "DiskImageResourcePNG (read) Error: Destination buffer has wrong dimensions!" );

  boost::scoped_array<uint8> buf( new uint8[ctx->cstride * bbox.width() * bbox.height()] );
  // Interlacing is causing problems when read line-by-line...I think it's
  // a bug in libpng.
  if( ctx->interlaced )
  {
    if( bbox.height() != rows() )
      vw_throw( NoImplErr() << "DiskImageResourcePNG: Reading interlaced files line-by-line is currently unsupported." );

    ctx->readall(buf);
  }
  else
  {
    // FIXME: Normal operation. Make this what happens all the time when
    // the libpng bug gets fixed. The bug in question causes the final
    // parts of the expanded, interlaced image (the ones we care about) to
    // seemingly 'lose' a row.

    // If our start line is a spot before the current line, we need to reopen the file.
    if(start_line < ctx->current_line)
      read_reset();
    if(start_line > ctx->current_line)
      ctx->advance(start_line - ctx->current_line);

    // Now read
    int32 offset = 0;
    while(ctx->current_line < end_line) {
      ctx->readline();

      // Copy the data over into a contiguous buffer, which is what
      // convert() expects when we call it below. This extra copy of the
      // data is undesirable, but probably not a huge performance hit
      // compared to reading the data from disk.
      std::memcpy(buf.get() + offset,
                  ctx->scanline.get() + ctx->cstride * bbox.min().x(),
                  bbox.width() * ctx->cstride);
      offset += bbox.width() * ctx->cstride;
    }

  }

  // Set up an image buffer for convert().
  ImageBuffer src;
  src.data = buf.get();
  src.format = m_format;
  // The number of rows we're reading may be different, seeing as how we
  // can read sequential lines.
  src.format.rows = bbox.height();
  src.format.cols = bbox.width();
  src.unpremultiplied = true;

  src.cstride = ctx->cstride;
  src.rstride = src.cstride * src.format.cols;
  src.pstride = src.rstride * src.format.rows;
  convert(dest, src, m_rescale);
}

void DiskImageResourcePNG::read_reset() const {
  m_ctx.reset( new DiskImageResourcePNG::vw_png_read_context( const_cast<DiskImageResourcePNG *>(this) ) );
}

unsigned DiskImageResourcePNG::num_comments() const {
  m_ctx->read_comments();
  return boost::numeric_cast<unsigned>(m_ctx->comments.size());
}

DiskImageResourcePNG::Comment const& DiskImageResourcePNG::get_comment( unsigned i ) const {
  m_ctx->read_comments();
  return m_ctx->comments[i];
}

std::string const& DiskImageResourcePNG::get_comment_key( unsigned i ) const {
  m_ctx->read_comments();
  return get_comment(i).key;
}

std::string const& DiskImageResourcePNG::get_comment_value( unsigned i ) const {
  m_ctx->read_comments();
  return get_comment(i).text;
}

/***********************************************************************
 ************************ STUFF FOR WRITING ****************************
***********************************************************************/

void DiskImageResourcePNG::create( std::string const& filename, ImageFormat const& format ) {
  DiskImageResourcePNG::create( filename, format, DiskImageResourcePNG::Options() );
}

void DiskImageResourcePNG::create( std::string const& filename, ImageFormat const& format, DiskImageResourcePNG::Options const& options )
{
  if(m_ctx.get() != NULL)
    vw_throw( IOErr() << "DiskImageResourcePNG: A file is already open." );

  m_filename = filename;
  m_format = format;

  m_ctx = boost::shared_ptr<vw_png_context>( new vw_png_write_context( const_cast<DiskImageResourcePNG *>(this), options ) );

  // Block writing is not supported
  m_block_size = Vector2i( cols(), rows() );
}

void DiskImageResourcePNG::write( ImageBuffer const& src, BBox2i const& bbox )
{
  vw_png_write_context *ctx = dynamic_cast<vw_png_write_context *>( m_ctx.get() );

  VW_ASSERT( bbox.width()==int(cols()) && bbox.height()==int(rows()),
             NoImplErr() << "DiskImageResourcePNG does not support partial writes." );
  VW_ASSERT( src.format.cols==uint32(cols()) && src.format.rows==uint32(rows()),
             ArgumentErr() << "DiskImageResourcePNG: Buffer has wrong dimensions in PNG write." );

  // Set up the image buffer and convert the data into this buffer.
  ImageBuffer dst;
  boost::scoped_array<uint8> buf(new uint8[ctx->cstride * bbox.width() * bbox.height()]);

  dst.data = buf.get();
  dst.format = m_format;
  dst.format.rows = bbox.height();
  dst.format.cols = bbox.width();
  dst.unpremultiplied = true;

  if (dst.format.channel_type != VW_CHANNEL_UINT16 &&
      dst.format.channel_type != VW_CHANNEL_INT16)
    dst.format.channel_type = VW_CHANNEL_UINT8;

  dst.cstride = num_channels(dst.format.pixel_format) * channel_size(dst.format.channel_type);
  dst.rstride = dst.cstride * dst.format.cols;
  dst.pstride = dst.rstride * dst.format.rows;

  convert(dst, src, m_rescale);

  // Write.
  ctx->write(dst);
}
