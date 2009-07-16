// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file DiskImageResourcePGM.cc
///
/// Provides support for the PGM file format.
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
#include <vw/FileIO/DiskImageResourcePGM.h>

using namespace vw;

void skip_any_comments( FILE * stream ) {
  char temp;
  temp = fgetc(stream);
  while ( temp == '#' ) {
    do {
      temp = fgetc(stream);
    } while ( temp != '\n' );
    temp = fgetc(stream);
  }
  fseek(stream,-1,SEEK_CUR);
}

// Constructors
DiskImageResourcePGM::DiskImageResourcePGM( std::string const& filename ) : DiskImageResource( filename ) {
  open( filename );
}

DiskImageResourcePGM::DiskImageResourcePGM( std::string const& filename, ImageFormat const& format ) : DiskImageResource( filename ) {
  create( filename, format );
}

// Bind the resource to a file for reading. Confirm that we can open
// the file and that it has a sane pixel format.
void DiskImageResourcePGM::open( std::string const& filename ) { 

  FILE* input_file = fopen(filename.c_str(), "r");
  if (!input_file) vw_throw( vw::IOErr() << "Failed to open \"" << filename << "\"." );

  char c_line[2048];

  // Reading version info
  skip_any_comments(input_file);
  fscanf(input_file,"%s",c_line);
  std::string string_temp(c_line);
  if ( string_temp != "P5" && string_temp != "P2" )
    vw_throw( IOErr() << "DiskImageResourcePGM: unsupported / or incorrect magic number identifer \"" << string_temp << "\"." );

  // Getting image width, height, and max gray value.
  int32 iwidth, iheight, imax;
  skip_any_comments(input_file);
  fscanf(input_file,"%d",&iwidth);
  skip_any_comments(input_file);
  fscanf(input_file,"%d",&iheight);
  skip_any_comments(input_file);
  fscanf(input_file,"%d",&imax);
  fgetpos(input_file, &m_image_data_position);
  fclose(input_file);

  // Checking bit sanity
  if (imax <= 0 || imax >= 65536 )
    vw_throw( IOErr() << "DiskImageResourcePGM: invalid bit type, must be in range of an unsigned 16 bit integer or less." );

  m_format.cols = iwidth;
  m_format.rows = iheight;
  m_format.planes = 1;

  if (imax < 2)
    m_format.channel_type = VW_CHANNEL_BOOL;
  else if ( imax < 256 )
    m_format.channel_type = VW_CHANNEL_UINT8;
  else
    m_format.channel_type = VW_CHANNEL_UINT16;
  
  // There can be only one!
  m_format.pixel_format = VW_PIXEL_GRAY;

}

// Bind the resource to a file for writing.
void DiskImageResourcePGM::create( std::string const& filename,
				   ImageFormat const& format ) {
  vw_throw( NoImplErr() << "I need an example of how to write before this is accomplished." );

}

// Read the disk image into the given buffer.
void DiskImageResourcePGM::read( ImageBuffer const& dest, BBox2i const& bbox )  const {

  VW_ASSERT( bbox.width()==int(cols()) && bbox.height()==int(rows()),
	     NoImplErr() << "DiskImageResourcePGM does not support partial reads." );
  VW_ASSERT( dest.format.cols==cols() && dest.format.rows==rows(),
	     IOErr() << "Buffer has wrong dimensions in PGM read." );
    
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
  if ( m_format.channel_type == VW_CHANNEL_UINT8 ) {
    uint8* image_data = new uint8[num_pixels];
    fread( image_data, sizeof(uint8), num_pixels, input_file );
    src.data = image_data;
    convert( dest, src, m_rescale );
    delete[] image_data;
  } else if ( m_format.channel_type == VW_CHANNEL_UINT16 ) {
    uint16* image_data = new uint16[num_pixels];
    fread( image_data, sizeof(uint16), num_pixels, input_file );
    src.data = image_data;
    convert( dest, src, m_rescale );
    delete[] image_data;
  } else if ( m_format.channel_type == VW_CHANNEL_BOOL ) {
    bool* image_data = new bool[num_pixels];
    fread( image_data, sizeof(bool), num_pixels, input_file );
    src.data = image_data;
    convert( dest, src, m_rescale );
    delete[] image_data;
  } else {
    // Probably never thrown
    vw_throw( NoImplErr() << "Unknown input channel type.");
  }
  
  fclose(input_file);
}

// Write the given buffer into the disk image.
void DiskImageResourcePGM::write( ImageBuffer const& src,
				  BBox2i const& bbox ) {
  vw_throw( NoImplErr() << "The PGM driver does not yet support creation of PGM files." );
}

// A FileIO hook to open a file for reading
DiskImageResource* DiskImageResourcePGM::construct_open( std::string const& filename ) {
  return new DiskImageResourcePGM( filename );
}

// A FileIO hook to open a file for writing
DiskImageResource* DiskImageResourcePGM::construct_create( std::string const& filename,
							   ImageFormat const& format ) {
  return new DiskImageResourcePGM( filename, format );
}
