// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file DiskImageResourcePBM.cc
///
/// Provides support for the PBM file format.
///

#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <vector>
#include <fstream>
#include <stdio.h>

#include <boost/algorithm/string.hpp>
using namespace boost;

#include <vw/Core/Exception.h>
#include <vw/Core/Debugging.h>
#include <vw/FileIO/DiskImageResourcePBM.h>

using namespace vw;

// Used to skip comment lines found in file
void skip_any_comments( FILE * stream ) {
  char temp;
  temp = fgetc(stream);
  // Sometimes this can land on just being before a new line
  while ( temp == '\n' ) 
    temp = fgetc(stream);
  // Clearing away comment
  while ( temp == '#' ) {
    do {
      temp = fgetc(stream);
    } while ( temp != '\n' );
    temp = fgetc(stream);
  }
  fseek(stream,-1,SEEK_CUR);
}
// Used to normalize an array of uint8s
void normalize( uint8* data, uint32 count, uint8 max_value ) {
  uint8* pointer = data;
  for ( uint32 i = 0; i < count; i++ ) {
    if ( *pointer > max_value )
      *pointer = 255;
    else {
      uint8 temp = *pointer;
      *pointer = uint8(uint16(temp)*255/uint16(max_value));
    }
    pointer++;
  }
}

// Constructors
DiskImageResourcePBM::DiskImageResourcePBM( std::string const& filename ) : DiskImageResource( filename ) {
  open( filename );
}

DiskImageResourcePBM::DiskImageResourcePBM( std::string const& filename, ImageFormat const& format ) : DiskImageResource( filename ) {
  create( filename, format );
}

// Bind the resource to a file for reading. Confirm that we can open
// the file and that it has a sane pixel format.
void DiskImageResourcePBM::open( std::string const& filename ) { 

  FILE* input_file = fopen(filename.c_str(), "r");
  if (!input_file) vw_throw( vw::IOErr() << "Failed to open \"" << filename << "\"." );

  char c_line[2048];

  // Reading version info
  skip_any_comments(input_file);
  fscanf(input_file,"%s",c_line);
  m_magic = std::string(c_line);
  if ( !(m_magic == "P6" || m_magic == "P5" || m_magic == "P4" ||
	 m_magic == "P3" || m_magic == "P2" || m_magic == "P1" ) )
    vw_throw( IOErr() << "DiskImageResourcePBM: unsupported / or incorrect magic number identifer \"" << m_magic << "\"." );

  // Getting image width, height, and max gray value.
  int32 iwidth, iheight;
  skip_any_comments(input_file);
  fscanf(input_file,"%d",&iwidth);
  skip_any_comments(input_file);
  fscanf(input_file,"%d",&iheight);
  if ( m_magic != "P1" && m_magic != "P4" ) {
    skip_any_comments(input_file);
    fscanf(input_file,"%d",&m_max_value);
  } else
    m_max_value = 1; // Binary image
  fgetpos(input_file, &m_image_data_position);
  fclose(input_file);

  // Checking bit sanity
  if (m_max_value <= 0 || m_max_value > 255 )
    vw_throw( IOErr() << "DiskImageResourcePBM: invalid bit type, Netpbm support 8 bit channel types and lower. Max requested: " << m_max_value );

  m_format.cols = iwidth;
  m_format.rows = iheight;
  m_format.planes = 1;

  if ( m_magic == "P1" || m_magic == "P4" ) {
    // Boolean File Type
    m_format.channel_type = VW_CHANNEL_BOOL;
    m_format.pixel_format = VW_PIXEL_GRAY;
  } else if ( m_magic == "P2" || m_magic == "P5" ) {
    // Grayscale image
    m_format.channel_type = VW_CHANNEL_UINT8;
    m_format.pixel_format = VW_PIXEL_GRAY;
  } else if ( m_magic == "P3" || m_magic == "P6" ) {
    // RGB image
    m_format.channel_type = VW_CHANNEL_UINT8;
    m_format.pixel_format = VW_PIXEL_RGB;
  } else
    vw_throw( IOErr() << "DiskImageResourcePBM: how'd you get here? Invalid magic number." );

}

// Read the disk image into the given buffer.
void DiskImageResourcePBM::read( ImageBuffer const& dest, BBox2i const& bbox )  const {

  VW_ASSERT( bbox.width()==int(cols()) && bbox.height()==int(rows()),
	     NoImplErr() << "DiskImageResourcePBM does not support partial reads." );
  VW_ASSERT( dest.format.cols==cols() && dest.format.rows==rows(),
	     IOErr() << "Buffer has wrong dimensions in PBM read." );
    
  FILE* input_file = fopen(m_filename.c_str(), "r");
  if (!input_file) vw_throw( IOErr() << "Failed to open \"" << m_filename << "\"." );
  fsetpos(input_file,&m_image_data_position);

  // Reading image data
  ImageBuffer src;
  int32 num_pixels = m_format.cols*m_format.rows;
  src.format = m_format;
  src.cstride = 1;
  src.rstride = m_format.cols;
  src.pstride = num_pixels;

  if ( m_magic == "P1" ) {
    // Bool ASCII
    bool* image_data = new bool[num_pixels];
    bool* point = image_data;
    char buffer;
    for ( int32 i = 0; i < num_pixels; i++ ) {
      fscanf( input_file, "%c", &buffer );
      if ( buffer == '1' )
	*point = true;
      else
	*point = false;
      point++;
    }
    src.data = image_data;
    convert( dest, src, m_rescale );
    delete[] image_data;
  } else if ( m_magic == "P2" ) {
    // Gray uint8 ASCII
    uint8* image_data = new uint8[num_pixels];
    uint8* point = image_data;
    for ( int32 i = 0; i < num_pixels; i++ ) {
      fscanf( input_file, "%hhu", point );
      point++;
    }
    normalize( image_data, num_pixels, m_max_value );
    src.data = image_data;
    convert( dest, src, m_rescale );
    delete[] image_data;
  } else if ( m_magic == "P3" ) {
    // RGB uint8 ASCII
    uint8* image_data = new uint8[num_pixels*3];
    uint8* point = image_data;
    for ( int32 i = 0; i < num_pixels; i++ ) {
      fscanf( input_file, "%hhu", point );
      point++;
      fscanf( input_file, "%hhu", point );
      point++;
      fscanf( input_file, "%hhu", point );
      point++;
    }
    normalize( image_data, num_pixels*3, m_max_value );
    src.data = image_data;
    src.cstride = 3;
    src.rstride = src.cstride*m_format.cols;
    src.pstride = src.rstride*m_format.rows;
    convert( dest, src, m_rescale );
    delete[] image_data;
  } else if ( m_magic == "P4" ) {
    // Bool Binary
    bool* image_data = new bool[num_pixels];
    fread( image_data, sizeof(bool), num_pixels, input_file );
    src.data = image_data;
    convert( dest, src, m_rescale );
    delete[] image_data;
  } else if ( m_magic == "P5" ) {
    // Gray uint8 Binary
    uint8* image_data = new uint8[num_pixels];
    fread( image_data, sizeof(uint8), num_pixels, input_file );
    normalize( image_data, num_pixels, m_max_value );
    src.data = image_data;
    convert( dest, src, m_rescale );
    delete[] image_data;
  } else if ( m_magic == "P6" ) {
    // RGB uint8 Binary
    uint8* image_data = new uint8[num_pixels*3];
    fread( image_data, sizeof(uint8), num_pixels*3, input_file );
    normalize( image_data, num_pixels*3, m_max_value );
    src.data = image_data;
    src.cstride = 3;
    src.rstride = src.cstride*m_format.cols;
    src.pstride = src.rstride*m_format.rows;
    convert( dest, src, m_rescale );
    delete[] image_data;
  } else 
    vw_throw( NoImplErr() << "Unknown input channel type." );

  fclose(input_file);
}

// Bind the resource to a file for writing.
void DiskImageResourcePBM::create( std::string const& filename,
				   ImageFormat const& format ) {

  if ( format.planes != 1 )
    vw_throw( NoImplErr() << "DiskImageResourcePBM doesn't support multi-plane images.");

  m_filename = filename;
  m_format = format;

  // Deciding output pixel format based on the file extension regardless of input type. If they choose an unknown extension, we'll guess the writing format.
  if ( boost::algorithm::iends_with( filename, ".pbm" ) || 
       boost::algorithm::iends_with( filename, ".PBM" ) ) {
    m_format.pixel_format = VW_PIXEL_UNKNOWN_MASKED;
    m_format.channel_type = VW_CHANNEL_BOOL;
    m_magic = "P4";
    m_max_value = 1;
  } else if ( boost::algorithm::iends_with( filename, ".pgm" ) ||
	      boost::algorithm::iends_with( filename, ".PGM" )  ) {
    m_format.pixel_format = VW_PIXEL_GRAY;
    m_format.channel_type = VW_CHANNEL_UINT8;
    m_magic = "P5";
    m_max_value = 255;
  } else if ( boost::algorithm::iends_with( filename, ".ppm" ) ||
	      boost::algorithm::iends_with( filename, ".PPM" ) ) {
    m_format.pixel_format = VW_PIXEL_RGB;
    m_format.channel_type = VW_CHANNEL_UINT8;
    m_magic = "P6";
    m_max_value = 255;
  } else {
    // Give up and guess based on image format

    // Deciding channel type (there are few options).
    if ( m_format.channel_type != VW_CHANNEL_BOOL )
      m_format.channel_type = VW_CHANNEL_UINT8;
    
    // Deciding pixel format
    switch ( format.pixel_format ) {
    case VW_PIXEL_UNKNOWN_MASKED:
    case VW_PIXEL_SCALAR_MASKED:
      m_format.pixel_format = VW_PIXEL_UNKNOWN_MASKED;
      m_format.channel_type = VW_CHANNEL_BOOL;
      m_magic = "P4";
      m_max_value = 1;
      break;
    case VW_PIXEL_UNKNOWN:
    case VW_PIXEL_SCALAR:
    case VW_PIXEL_GRAY:
    case VW_PIXEL_GRAYA:
    case VW_PIXEL_GRAY_MASKED:
    case VW_PIXEL_GRAYA_MASKED:
    case VW_PIXEL_GENERIC_1_CHANNEL:
      m_format.pixel_format = VW_PIXEL_GRAY;
      m_magic = "P5";
      m_max_value = 255;
      break;
    default:
      m_format.pixel_format = VW_PIXEL_RGB;
      m_magic = "P6";
      m_max_value = 255;
      break;
    }

  }

  FILE* output_file = fopen(filename.c_str(), "w");
  fprintf( output_file, "%s\n", m_magic.c_str() );
  fprintf( output_file, "%d\n", m_format.cols );
  fprintf( output_file, "%d\n", m_format.rows );
  if ( m_magic != "P1" && m_magic != "P4" ) 
    fprintf( output_file, "%d\n", m_max_value );
  fgetpos( output_file, &m_image_data_position );
  fclose( output_file );
}

// Write the given buffer into the disk image.
void DiskImageResourcePBM::write( ImageBuffer const& src,
				  BBox2i const& bbox ) {
  VW_ASSERT( bbox.width()==int(cols()) && bbox.height()==int(rows()),
	     NoImplErr() << "DiskImageResourcePBM does not support partial writes." );
  VW_ASSERT( src.format.cols==cols() && src.format.rows==rows(),
	     IOErr() << "Buffer has wrong dimensions in PBM write." );

  FILE* output_file = fopen(m_filename.c_str(), "a");
  fsetpos(output_file,&m_image_data_position);
  
  ImageBuffer dst;
  int32 num_pixels = m_format.cols*m_format.rows;
  dst.format = m_format;
  if ( m_magic == "P6" || m_magic == "P3" ) 
    dst.cstride = 3;
  else
    dst.cstride = 1;
  dst.rstride = dst.cstride*cols();
  dst.pstride = dst.rstride*rows();
  
  // Ready to start writing
  if ( m_magic == "P1" ) {
    // Bool ASCII
    bool* image_data = new bool[num_pixels];
    bool* point = image_data;
    dst.data = image_data;
    convert( dst, src, m_rescale );
    for ( int32 i = 0; i < num_pixels; i++ ) {
      if ( *point == true )
	fprintf( output_file, "1 " );
      else
	fprintf( output_file, "0 " );
      point++;
    }
    delete[] image_data;
  } else if ( m_magic == "P2" ) {
    // Gray uint8 ASCII
    uint8* image_data = new uint8[num_pixels];
    uint8* point = image_data;
    dst.data = image_data;
    convert( dst, src, m_rescale );
    for ( int32 i = 0; i < num_pixels; i++ ) {
      fprintf( output_file, "%hu ", *point );
      point++;
    }
    delete[] image_data;
  } else if ( m_magic == "P3" ) {
    // RGB uint8 ASCII
    uint8* image_data = new uint8[num_pixels*3];
    uint8* point = image_data;
    dst.data = image_data;
    convert( dst, src, m_rescale );
    for ( int32 i = 0; i < num_pixels*3; i++ ) {
      fprintf( output_file, "%hu ", *point );
      point++;
    }
    delete[] image_data;
  } else if ( m_magic == "P4" ) {
    // Bool Binary
    bool* image_data = new bool[num_pixels];
    dst.data = image_data;
    convert( dst, src, m_rescale );
    fwrite( image_data, sizeof(bool), num_pixels, output_file );
    delete[] image_data;
  } else if ( m_magic == "P5" ) {
    // Gray uint8 Binary
    uint8* image_data = new uint8[num_pixels];
    dst.data = image_data;
    convert( dst, src, m_rescale );
    fwrite( image_data, sizeof(uint8), num_pixels, output_file );
    delete[] image_data;
  } else if ( m_magic == "P6" ) {
    // RGB uint8 Binary
    uint8* image_data = new uint8[num_pixels*3];
    dst.data = image_data;
    convert( dst, src, m_rescale );
    fwrite( image_data, sizeof(uint8), num_pixels*3, output_file );
    delete[] image_data;
  }

  fclose( output_file );
}

// A FileIO hook to open a file for reading
DiskImageResource* DiskImageResourcePBM::construct_open( std::string const& filename ) {
  return new DiskImageResourcePBM( filename );
}

// A FileIO hook to open a file for writing
DiskImageResource* DiskImageResourcePBM::construct_create( std::string const& filename,
							   ImageFormat const& format ) {
  return new DiskImageResourcePBM( filename, format );
}
