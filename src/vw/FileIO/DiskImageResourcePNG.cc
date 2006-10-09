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
/// TODO Currently we do not support 1-bit images, though we could.
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
#include <vw/Image/GenericImageBuffer.h>
#include <vw/FileIO/DiskImageResourcePNG.h>

namespace vw {

  // *******************************************************************
  // The PNG library interface class
  // *******************************************************************

  class DiskImageResourceInfoPNG {
    std::string m_filename;
    std::fstream *m_fs_ptr;
    bool m_readable;
    GenericImageFormat &m_format;

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
    
    void read_init( png_structp &png_ptr, png_infop &info_ptr, png_infop &end_ptr ) const {
      char header[8];
      m_fs_ptr->read( header, 8 );
      if( !*m_fs_ptr || png_sig_cmp((png_bytep)header,0,8) )
        throw IOErr() << "Input file " << m_filename << " is not a valid PNG file.";
      
      png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,0,0,0);
      if( !png_ptr ) throw IOErr() << "Unable to initialize PNG reader.";
      
      info_ptr = png_create_info_struct(png_ptr);
      if( !info_ptr ) throw IOErr() << "Unable to initialize PNG reader.";
      
      end_ptr = png_create_info_struct(png_ptr);
      if (!end_ptr) throw IOErr() << "Unable to initialize PNG reader.";
      
      png_set_read_fn(png_ptr, (voidp)m_fs_ptr, read_data);
      
      png_set_sig_bytes(png_ptr, 8);
  
      png_read_info(png_ptr, info_ptr);
    }

    void read_cleanup( png_structp &png_ptr, png_infop &info_ptr, png_infop &end_ptr ) const {
      if( m_fs_ptr ) m_fs_ptr->seekg( 0 );
      if( png_ptr ) png_destroy_read_struct( &png_ptr, &info_ptr, &end_ptr );
    }
    
    void write_init( png_structp &png_ptr, png_infop &info_ptr ) {
      png_ptr = 0;
      info_ptr = 0;
      try {
        png_ptr = png_create_write_struct( PNG_LIBPNG_VER_STRING, 0, 0, 0 );
        if( !png_ptr ) throw IOErr() << "Unable to initialize PNG writer.";
        
        png_set_write_fn( png_ptr, (voidp)m_fs_ptr, write_data, flush_data );
        
        info_ptr = png_create_info_struct(png_ptr);
        if( !info_ptr ) throw IOErr() << "Unable to initialize PNG reader.";
        
        int color_type = 0;
        switch( m_format.pixel_format ) {
        case VW_PIXEL_GRAY:  color_type = PNG_COLOR_TYPE_GRAY;       break;
        case VW_PIXEL_GRAYA: color_type = PNG_COLOR_TYPE_GRAY_ALPHA; break;
        case VW_PIXEL_RGB:   color_type = PNG_COLOR_TYPE_RGB;        break;
        case VW_PIXEL_RGBA:  color_type = PNG_COLOR_TYPE_RGB_ALPHA;  break;
        default: throw LogicErr() << "Unexpected pixel format in PNG write.";
        }
        
        png_set_IHDR( png_ptr, info_ptr, m_format.cols, m_format.rows, (m_format.channel_type==VW_CHANNEL_UINT8) ? 8 : 16,
                      color_type, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT );
        png_write_info( png_ptr, info_ptr );
      }
      catch(...) {
        if( png_ptr ) png_destroy_write_struct(&png_ptr,&info_ptr);
        throw;
      }
    }

    void write_cleanup( png_structp &png_ptr, png_infop &info_ptr ) {
      if( m_fs_ptr ) m_fs_ptr->seekg( 0 );
      if( png_ptr ) png_destroy_write_struct( &png_ptr, &info_ptr );
    }

  public:
    DiskImageResourceInfoPNG( GenericImageFormat &format )
      : m_filename(), m_fs_ptr(0), m_readable(false), m_format(format) {}

    ~DiskImageResourceInfoPNG() {
      delete m_fs_ptr;
    }

    void open( std::string const& filename ) {
      m_readable = false;
      m_filename = filename;
      m_fs_ptr = new std::fstream( m_filename.c_str(), std::ios_base::in | std::ios_base::binary );
      if( !m_fs_ptr || !*m_fs_ptr )
        throw IOErr() << "Unable to open input file " << m_filename << ".";
      png_structp png_ptr;
      png_infop info_ptr, end_ptr;
      read_init( png_ptr, info_ptr, end_ptr );
      
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
      default:
        throw NoImplErr() << "Unsupported PNG pixel format.";
      }
      
      switch( depth ) {
      case 8:
        m_format.channel_type = VW_CHANNEL_UINT8;
        break;
      case 16:
        m_format.channel_type = VW_CHANNEL_UINT16;
        break;
      default:
        throw NoImplErr() << "Unsupported PNG pixel depth.";
      }
      
      read_cleanup( png_ptr, info_ptr, end_ptr );
      m_readable = true;
    }

    void create( std::string const& filename, GenericImageFormat const& format ) {
      m_format = format;
      // If we're asked for a multi-plane scalar format, pretend we've
      // been asked for the corresponding single-plane multi-channel format.
      if( m_format.pixel_format == VW_PIXEL_SCALAR ) {
        switch( m_format.planes ) {
        case 1: m_format.pixel_format = VW_PIXEL_GRAY;  break;
        case 2: m_format.pixel_format = VW_PIXEL_GRAYA; break;
        case 3: m_format.pixel_format = VW_PIXEL_RGB;   break;
        case 4: m_format.pixel_format = VW_PIXEL_RGBA;  break;
        default: throw ArgumentErr() << "PNG files do not support more than four planes.";
        }
        m_format.planes = 1;
      }
      if( m_format.planes != 1 ) throw ArgumentErr() << "PNG files to not support multiple images.";
      if( m_format.pixel_format!=VW_PIXEL_GRAY && m_format.pixel_format!=VW_PIXEL_GRAYA && 
          m_format.pixel_format!=VW_PIXEL_RGB  && m_format.pixel_format!=VW_PIXEL_RGBA )
        throw ArgumentErr() << "Unrecognized pixel format for PNG image.";
      if( m_format.channel_type!=VW_CHANNEL_UINT16 ) m_format.channel_type = VW_CHANNEL_UINT8;
      m_readable = false;
      m_filename = filename;
      m_fs_ptr = new std::fstream( m_filename.c_str(), std::ios::out | std::ios_base::binary );
      if( !m_fs_ptr || !*m_fs_ptr )
        throw IOErr() << "Unable to open output file \"" << m_filename << "\".";
    }

    void read( GenericImageBuffer const& dest, BBox2i bbox ) const {
      png_structp png_ptr;
      png_infop info_ptr, end_ptr;
      read_init( png_ptr, info_ptr, end_ptr );
      VW_ASSERT( bbox.min().x()>=0 && bbox.min().y()>=0 && bbox.max().x()<=m_format.cols && bbox.max().y()<=m_format.rows,
                 ArgumentErr() << "Requested bounding box extends beyond disk image boundaries in PNG read." );
      VW_ASSERT( dest.format.cols==bbox.width() && dest.format.rows==bbox.height(),
                 ArgumentErr() << "Destination buffer has wrong dimensions in PNG read." );

      GenericImageBuffer src;
      int Bpp = (m_format.channel_type==VW_CHANNEL_UINT16) ? 2 : 1;
      int channels = png_get_channels( png_ptr, info_ptr );
      std::vector<uint8> buffer(m_format.cols*m_format.rows*channels*Bpp);

      src.format = m_format;
      src.data = &buffer[0];
      src.cstride = Bpp*channels;
      src.rstride = m_format.cols*Bpp*channels;
      src.pstride = m_format.rows*m_format.cols*Bpp*channels;
      
      std::vector<png_bytep> row_pointers( m_format.rows );
      for( int i=0; i<m_format.rows; ++i )
        row_pointers[i] = (png_bytep)(src.data) + i*src.rstride;

      src.format.cols = bbox.width();
      src.format.rows = bbox.height();
      src.data = (uint8*)(src.data) + bbox.min().x()*src.cstride + bbox.min().y()*src.rstride;

      png_read_image( png_ptr, &row_pointers[0] );
      convert( dest, src );
      read_cleanup( png_ptr, info_ptr, end_ptr );
    }

    void write( GenericImageBuffer const& src ) {
      png_structp png_ptr;
      png_infop info_ptr;
      write_init( png_ptr, info_ptr );
      VW_ASSERT( src.format.cols==m_format.cols && src.format.rows==m_format.rows, IOErr() << "Buffer has wrong dimensions in PNG write." );
      GenericImageBuffer dst;
      int Bpp = (m_format.channel_type==VW_CHANNEL_UINT16) ? 2 : 1;
      int channels = png_get_channels( png_ptr, info_ptr );
      std::vector<uint8> buffer(m_format.cols*m_format.rows*channels*Bpp);
      dst.data = &buffer[0];
      dst.format = m_format;
      dst.cstride = Bpp*channels;
      dst.rstride = m_format.cols*Bpp*channels;
      dst.pstride = m_format.rows*m_format.cols*Bpp*channels;
      
      std::vector<png_bytep> row_pointers( m_format.rows );
      for( int i=0; i<m_format.rows; ++i )
        row_pointers[i] = (png_bytep)(dst.data) + i*dst.rstride;

      convert( dst, src );
      png_write_image( png_ptr, &row_pointers[0] );
      png_write_end( png_ptr, info_ptr );
      write_cleanup( png_ptr, info_ptr );
    }
  };

} // namespace vw


// *********************************************************************
// Public interface function definitions
// *********************************************************************

vw::DiskImageResourcePNG::DiskImageResourcePNG( std::string const& filename )
  : m_info( new DiskImageResourceInfoPNG(m_format) ) {
  open( filename );
}

vw::DiskImageResourcePNG::DiskImageResourcePNG( std::string const& filename, 
                                                GenericImageFormat const& format )
  : m_info( new DiskImageResourceInfoPNG(m_format) ) {
  create( filename, format );
}

vw::DiskImageResourcePNG::~DiskImageResourcePNG() {
}

void vw::DiskImageResourcePNG::open( std::string const& filename ) {
  m_info->open( filename );
}

void vw::DiskImageResourcePNG::create( std::string const& filename, 
                                       GenericImageFormat const& format ) {
  m_info->create( filename, format );
}

void vw::DiskImageResourcePNG::read( GenericImageBuffer const& dest ) const {
  m_info->read( dest, BBox2i(0,0,cols(),rows()) );
}

void vw::DiskImageResourcePNG::read( GenericImageBuffer const& dest, BBox2i bbox ) const {
  m_info->read( dest, bbox );
}

void vw::DiskImageResourcePNG::write( GenericImageBuffer const& src ) {
  m_info->write( src );
}

vw::DiskImageResource* vw::DiskImageResourcePNG::construct_open( std::string const& filename ) {
  return new DiskImageResourcePNG( filename );
}

vw::DiskImageResource* vw::DiskImageResourcePNG::construct_create( std::string const& filename,
                                                                   GenericImageFormat const& format ) {
  return new DiskImageResourcePNG( filename, format );
}
