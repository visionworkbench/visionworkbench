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

#include <vw/Core/Exception.h>
#include <vw/Image/Manipulation.h>
#include <vw/FileIO/DiskImageResourcePNG.h>

namespace vw {

  // *******************************************************************
  // The PNG library interface class
  // *******************************************************************

  class DiskImageResourceInfoPNG {
    std::string m_filename;
    bool m_readable;
    bool m_palette_based;
    bool m_use_palette_indices;
    ImageFormat &m_format;
    std::vector<DiskImageResourcePNG::Comment> comments;
    ImageView<PixelRGBA<uint8> > m_palette;

    static void read_data( png_structp png_ptr, png_bytep data, png_size_t length ) {
      std::fstream *fs = (std::fstream*) png_get_io_ptr(png_ptr);
      fs->read( (char*)data, length );
    }

    static void write_data( png_structp png_ptr, png_bytep data, png_size_t length ) {
      std::fstream *fs = (std::fstream*) png_get_io_ptr(png_ptr);
      fs->write( (char*)data, length );
    }

    static void flush_data( png_structp png_ptr ) {
        std::fstream *fs = (std::fstream*) png_get_io_ptr(png_ptr);
        fs->flush();
    }
    
    void read_init( std::fstream &fs, png_structp &png_ptr, png_infop &info_ptr, png_infop &end_ptr ) {
      char header[8];
      fs.read( header, 8 );
      if ( !fs || png_sig_cmp((png_bytep)header,0,8) )
        vw_throw( IOErr() << "Input file " << m_filename << " is not a valid PNG file." );
      
      png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,0,0,0);
      if ( !png_ptr ) vw_throw( IOErr() << "Unable to initialize PNG reader." );
      
      info_ptr = png_create_info_struct(png_ptr);
      if ( !info_ptr ) vw_throw( IOErr() << "Unable to initialize PNG reader." );
      
      end_ptr = png_create_info_struct(png_ptr);
      if (!end_ptr) vw_throw( IOErr() << "Unable to initialize PNG reader." );
      
      png_set_read_fn(png_ptr, (voidp)&fs, read_data);
      
      png_set_sig_bytes(png_ptr, 8);
  
      png_read_info(png_ptr, info_ptr);
    }

    void read_comments( png_structp &png_ptr, png_infop &info_ptr ) {
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
          vw_out(DebugMessage) << "Unsupported PNG comment type in PNG read!" << std::endl;
          continue;
        }
        comments.push_back( c );
      }
    }

    void read_cleanup( png_structp &png_ptr, png_infop &info_ptr, png_infop &end_ptr ) const {
      if ( png_ptr ) png_destroy_read_struct( &png_ptr, &info_ptr, &end_ptr );
    }
    
    void write_init( std::fstream& fs, png_structp &png_ptr, png_infop &info_ptr ) {
      png_ptr = 0;
      info_ptr = 0;
      png_ptr = png_create_write_struct( PNG_LIBPNG_VER_STRING, 0, 0, 0 );
      if ( !png_ptr ) vw_throw( IOErr() << "Unable to initialize PNG writer." );
        
      png_set_write_fn( png_ptr, (voidp)&fs, write_data, flush_data );
        
      info_ptr = png_create_info_struct(png_ptr);
      if ( !info_ptr ) {
        png_destroy_write_struct(&png_ptr,&info_ptr);
        vw_throw( IOErr() << "Unable to initialize PNG writer." );
      }

      int color_type = 0;
      if( m_palette_based ) {
        color_type = PNG_COLOR_TYPE_PALETTE;
        if ( m_use_palette_indices ) {
          png_colorp palette = (png_colorp) png_malloc( png_ptr, m_palette.cols() * sizeof(png_color) );
          png_bytep alpha = (png_bytep) png_malloc( png_ptr, m_palette.cols() * sizeof(png_byte) );
          for ( int i=0; i<(int)m_palette.cols(); ++i ) {
            palette[i].red = m_palette(i,0).r();
            palette[i].green = m_palette(i,0).g();
            palette[i].blue = m_palette(i,0).b();
            alpha[i] = m_palette(i,0).a();
          }
          png_set_PLTE( png_ptr, info_ptr, palette, m_palette.cols() );
          png_set_tRNS( png_ptr, info_ptr, alpha, m_palette.cols(), 0 );
        }
      }
      else {
        switch( m_format.pixel_format ) {
        case VW_PIXEL_GRAY:  color_type = PNG_COLOR_TYPE_GRAY;       break;
        case VW_PIXEL_GRAYA: color_type = PNG_COLOR_TYPE_GRAY_ALPHA; break;
        case VW_PIXEL_RGB:   color_type = PNG_COLOR_TYPE_RGB;        break;
        case VW_PIXEL_RGBA:  color_type = PNG_COLOR_TYPE_RGB_ALPHA;  break;
        default:
          png_destroy_write_struct(&png_ptr,&info_ptr);
          vw_throw( LogicErr() << "Unexpected pixel format (" << m_format.pixel_format << ") in PNG write." );
        }
      }
        
      png_set_IHDR( png_ptr, info_ptr, m_format.cols, m_format.rows, (m_format.channel_type==VW_CHANNEL_UINT8) ? 8 : 16,
                    color_type, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT );
      png_write_info( png_ptr, info_ptr );
    }

    void write_cleanup( png_structp &png_ptr, png_infop &info_ptr ) {
      if ( ! png_ptr ) return;
      if( m_palette_based && m_use_palette_indices ) {
        png_colorp palette;
        png_bytep alpha;
        int tmp;
        png_get_PLTE( png_ptr, info_ptr, &palette, &tmp );
      }
      png_destroy_write_struct( &png_ptr, &info_ptr );
    }

  public:
    DiskImageResourceInfoPNG( ImageFormat &format )
      : m_filename(), m_readable(false), m_palette_based(false), m_use_palette_indices(false), m_format(format) {}

    ~DiskImageResourceInfoPNG() { }

    void open( std::string const& filename ) {
      m_readable = false;
      m_filename = filename;
      std::fstream fs( m_filename.c_str(), std::ios_base::in | std::ios_base::binary );
      if ( !fs ) vw_throw( IOErr() << "Unable to open input file " << m_filename << "." );
      png_structp png_ptr;
      png_infop info_ptr, end_ptr;
      read_init( fs, png_ptr, info_ptr, end_ptr );
      read_comments( png_ptr, info_ptr );
      
      m_format.cols = png_get_image_width(png_ptr, info_ptr);
      m_format.rows = png_get_image_height(png_ptr, info_ptr);
      m_format.planes = 1;
      int depth = png_get_bit_depth(png_ptr, info_ptr);
      int type = png_get_color_type(png_ptr, info_ptr);
      
      switch( type ) {
      case PNG_COLOR_TYPE_GRAY:
        m_format.pixel_format = VW_PIXEL_GRAY;
        break;
      case PNG_COLOR_TYPE_GRAY_ALPHA:
        m_format.pixel_format = VW_PIXEL_GRAYA;
        break;
      case PNG_COLOR_TYPE_RGB:
        m_format.pixel_format = VW_PIXEL_RGB;
        break;
      case PNG_COLOR_TYPE_RGB_ALPHA:
        m_format.pixel_format = VW_PIXEL_RGBA;
        break;
      case PNG_COLOR_TYPE_PALETTE: {
        m_palette_based = true;
        if ( png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS) )
          m_format.pixel_format = VW_PIXEL_RGBA;
        else m_format.pixel_format = VW_PIXEL_RGB;
        break;
      }
      default:
        vw_throw( NoImplErr() << "Unsupported PNG pixel format (" << type << ")." );
      }
      
      switch( depth ) {
      case 8:
        m_format.channel_type = VW_CHANNEL_UINT8;
        break;
      case 16:
        m_format.channel_type = VW_CHANNEL_UINT16;
        break;
      default:
        vw_throw( NoImplErr() << "Unsupported PNG pixel depth (" << depth << ")." );
      }
      
      read_cleanup( png_ptr, info_ptr, end_ptr );
      m_readable = true;
    }

    void create( std::string const& filename, ImageFormat const& format ) {
      m_format = format;
      // If we're asked for a multi-plane scalar format, pretend we've
      // been asked for the corresponding single-plane multi-channel format.
      if ( m_format.pixel_format == VW_PIXEL_SCALAR ) {
        switch( m_format.planes ) {
        case 1: m_format.pixel_format = VW_PIXEL_GRAY;  break;
        case 2: m_format.pixel_format = VW_PIXEL_GRAYA; break;
        case 3: m_format.pixel_format = VW_PIXEL_RGB;   break;
        case 4: m_format.pixel_format = VW_PIXEL_RGBA;  break;
        default: vw_throw( ArgumentErr() << "PNG files do not support more than four planes." );
        }
        m_format.planes = 1;
      }
      if ( m_format.planes != 1 ) vw_throw( ArgumentErr() << "PNG files do not support multiple images." );
      if ( m_format.pixel_format!=VW_PIXEL_GRAY && m_format.pixel_format!=VW_PIXEL_GRAYA && 
          m_format.pixel_format!=VW_PIXEL_RGB  && m_format.pixel_format!=VW_PIXEL_RGBA )
        vw_throw( ArgumentErr() << "Unrecognized pixel format (" << m_format.pixel_format << ") for PNG image." );
      if ( m_format.channel_type!=VW_CHANNEL_UINT16 ) m_format.channel_type = VW_CHANNEL_UINT8;
      m_readable = false;
      m_filename = filename;
    }

    void read( ImageBuffer const& dest, BBox2i bbox ) {
      std::fstream fs( m_filename.c_str(), std::ios_base::in | std::ios_base::binary );
      if ( !fs ) vw_throw( IOErr() << "Unable to open input file " << m_filename << "." );
      png_structp png_ptr;
      png_infop info_ptr, end_ptr;
      read_init( fs, png_ptr, info_ptr, end_ptr );
      VW_ASSERT( bbox.min().x()>=0 && bbox.min().y()>=0 && bbox.max().x()<=int(m_format.cols) && bbox.max().y()<=int(m_format.rows),
                 ArgumentErr() << "Requested bounding box extends beyond disk image boundaries in PNG read." );
      VW_ASSERT( int(dest.format.cols)==bbox.width() && int(dest.format.rows)==bbox.height(),
                 ArgumentErr() << "Destination buffer has wrong dimensions in PNG read." );

      ImageBuffer src;
      unsigned Bpp = (m_format.channel_type==VW_CHANNEL_UINT16) ? 2 : 1;
      unsigned channels = png_get_channels( png_ptr, info_ptr );

      // Handle palette-based images
      if( ! m_use_palette_indices ) {
        int type = png_get_color_type(png_ptr, info_ptr);
        if ( type == PNG_COLOR_TYPE_PALETTE ) {
          png_set_palette_to_rgb(png_ptr);
          if ( png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS) ) {
            png_set_tRNS_to_alpha(png_ptr);
            channels = 4;
          }
          else channels = 3;
        }
      }

      // This is a terrible hack for detecting little-endian architectures.
      if ( m_format.channel_type==VW_CHANNEL_UINT16 ) {
        uint16 x = 1;
        if ( *(uint8*)&x == 1 ) png_set_swap( png_ptr );
      }

      src.format = m_format;
      src.cstride = Bpp*channels;
      src.rstride = m_format.cols*Bpp*channels;
      src.pstride = m_format.rows*m_format.cols*Bpp*channels;
      src.unpremultiplied = true;
      src.format.cols = bbox.width();

      if( png_get_interlace_type( png_ptr, info_ptr ) != PNG_INTERLACE_NONE || 
          m_format.cols * m_format.rows * channels * Bpp < (2<<20) ) {
        // Handle interlaced and small images all at once.
        std::vector<uint8> buffer(m_format.cols*m_format.rows*channels*Bpp);
        src.data = &buffer[0] + bbox.min().x()*src.cstride + bbox.min().y()*src.rstride;
        src.format.rows = bbox.height();
        std::vector<png_bytep> row_pointers( m_format.rows );
        for ( int32 i=0; i<m_format.rows; ++i )
          row_pointers[i] = (png_bytep)(&buffer[0]) + i*src.rstride;
        png_read_image( png_ptr, &row_pointers[0] );
        convert( dest, src );
      }
      else {
        // Handle large non-interlaced images one row at a time.
        // FIXME We should really do multiple rows at a time if 
        // there are few columns, to avoid overhead in convert().
        std::vector<uint8> buffer(m_format.cols*channels*Bpp);
        src.data = &buffer[0] + bbox.min().x()*src.cstride;
        src.format.rows = 1;
        ImageBuffer dest_row = dest;
        dest_row.format.rows = 1;
        for( int row=0; row<bbox.min().y(); ++row ) {
          // Skip un-needed rows...
          png_read_row( png_ptr, &buffer[0], 0 );
        }
        for( int row=0; row<bbox.height(); ++row ) {
          png_read_row( png_ptr, &buffer[0], 0 );
          convert( dest_row, src );
          dest_row.data = (uint8*)dest_row.data + dest_row.rstride;
        }
      }

      read_cleanup( png_ptr, info_ptr, end_ptr );
    }

    void write( ImageBuffer const& src ) {
      std::fstream fs( m_filename.c_str(), std::ios::out | std::ios_base::binary );
      if ( !fs ) vw_throw( IOErr() << "Unable to open output file \"" << m_filename << "\"." );
      png_structp png_ptr;
      png_infop info_ptr;
      write_init( fs, png_ptr, info_ptr );
      VW_ASSERT( src.format.cols==m_format.cols && src.format.rows==m_format.rows, IOErr() << "Buffer has wrong dimensions in PNG write." );
      ImageBuffer dst;
      unsigned Bpp = (m_format.channel_type==VW_CHANNEL_UINT16) ? 2 : 1;
      unsigned channels = png_get_channels( png_ptr, info_ptr );
      std::vector<uint8> buffer(m_format.cols*m_format.rows*channels*Bpp);
      dst.data = &buffer[0];
      dst.format = m_format;
      dst.cstride = Bpp*channels;
      dst.rstride = m_format.cols*Bpp*channels;
      dst.pstride = m_format.rows*m_format.cols*Bpp*channels;
      dst.unpremultiplied = true;

      std::vector<png_bytep> row_pointers( m_format.rows );
      for ( int32 i=0; i<m_format.rows; ++i )
        row_pointers[i] = (png_bytep)(dst.data) + i*dst.rstride;

      // This is a terrible hack for detecting little-endian architectures.
      if ( m_format.channel_type==VW_CHANNEL_UINT16 ) {
        uint16 x = 1;
        if ( *(uint8*)&x == 1 ) png_set_swap( png_ptr );
      }

      convert( dst, src );
      png_write_image( png_ptr, &row_pointers[0] );
      png_write_end( png_ptr, info_ptr );
      write_cleanup( png_ptr, info_ptr );
    }

    unsigned num_comments() const {
      return comments.size();
    }

    DiskImageResourcePNG::Comment const& get_comment( unsigned i ) const {
      return comments[i];
    }

    bool is_palette_based() const {
      return m_palette_based;
    }

    ImageView<PixelRGBA<uint8> > get_palette() const {
      return copy( m_palette );
    }

    void set_palette( ImageView<PixelRGBA<uint8> > const& palette ) {
      m_palette = copy( palette );
      m_palette_based = true;
    }

    void set_use_palette_indices() {
      if( ! m_palette_based )
        vw_throw( IOErr() << "PNG file is not palette-based!" );
      m_use_palette_indices = true;
      m_format.pixel_format = VW_PIXEL_SCALAR;
      m_format.channel_type = VW_CHANNEL_UINT8;
    }
  };

} // namespace vw


// *********************************************************************
// Public interface function definitions
// *********************************************************************

vw::DiskImageResourcePNG::DiskImageResourcePNG( std::string const& filename )
  : DiskImageResource( filename ), m_info( new DiskImageResourceInfoPNG(m_format) ) {
  open( filename );
}

vw::DiskImageResourcePNG::DiskImageResourcePNG( std::string const& filename, 
                                                ImageFormat const& format )
  : DiskImageResource( filename ), m_info( new DiskImageResourceInfoPNG(m_format) ) {
  create( filename, format );
}

vw::DiskImageResourcePNG::~DiskImageResourcePNG() {
}

void vw::DiskImageResourcePNG::open( std::string const& filename ) {
  m_info->open( filename );
}

void vw::DiskImageResourcePNG::create( std::string const& filename, 
                                       ImageFormat const& format ) {
  m_info->create( filename, format );
}

void vw::DiskImageResourcePNG::read( ImageBuffer const& dest, BBox2i const& bbox ) const {
  m_info->read( dest, bbox );
}

void vw::DiskImageResourcePNG::write( ImageBuffer const& src, BBox2i const& bbox ) {
  VW_ASSERT( bbox.width()==int(cols()) && bbox.height()==int(rows()),
             NoImplErr() << "DiskImageResourcePNG does not support partial writes." );
  m_info->write( src );
}

vw::DiskImageResource* vw::DiskImageResourcePNG::construct_open( std::string const& filename ) {
  return new DiskImageResourcePNG( filename );
}

vw::DiskImageResource* vw::DiskImageResourcePNG::construct_create( std::string const& filename,
                                                                   ImageFormat const& format ) {
  return new DiskImageResourcePNG( filename, format );
}

unsigned vw::DiskImageResourcePNG::num_comments() const {
  return m_info->num_comments();
}

vw::DiskImageResourcePNG::Comment const& vw::DiskImageResourcePNG::get_comment( unsigned i ) const {
  return m_info->get_comment(i);
}

bool vw::DiskImageResourcePNG::is_palette_based() const {
  return m_info->is_palette_based();
}

vw::ImageView<vw::PixelRGBA<vw::uint8> > vw::DiskImageResourcePNG::get_palette() const {
  return m_info->get_palette();
}

void vw::DiskImageResourcePNG::set_palette( ImageView<PixelRGBA<uint8> > const& palette ) {
  m_info->set_palette( palette );
}

void vw::DiskImageResourcePNG::set_use_palette_indices() {
  m_info->set_use_palette_indices();
}

