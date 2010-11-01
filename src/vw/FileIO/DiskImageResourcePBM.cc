// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
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
#include <cstdio>

#include <boost/scoped_array.hpp>
#include <boost/algorithm/string.hpp>
using namespace boost;

#include <vw/Core/Exception.h>
#include <vw/Core/Debugging.h>
#include <vw/FileIO/DiskImageResourcePBM.h>

using namespace vw;
using std::fstream;
using std::ifstream;
using std::ofstream;

// Default PBM file creation options, and related functions.
namespace {

// Used to skip comment lines found in file
// Will skip all whitespace, then skip blocks delimited by # and \n, then repeat
void skip_any_comments( std::istream& f ) {
  while (isspace(f.peek()))
    f.ignore();
  while (f.peek() == '#') {
    f.ignore(1024, '\n');
    while (isspace(f.peek()))
      f.ignore();
  }
}

// Used to normalize an array of uint8s
void normalize( scoped_array<uint8>& data, uint32 count, uint8 max_value ) {
  uint8* pointer = data.get();
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

static bool default_ascii = false;

} // end anonymous

void DiskImageResourcePBM::default_to_ascii(bool ascii) {
  default_ascii = ascii;
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

  ifstream input(filename.c_str(), fstream::in|fstream::binary);

  if (!input.is_open())
    vw_throw( vw::ArgumentErr() << "DiskImageResourcePBM: Failed to open \"" << filename << "\"." );

  // Reading version info
  input >> m_magic;
  if ( !(m_magic == "P6" || m_magic == "P5" || m_magic == "P4" ||
         m_magic == "P3" || m_magic == "P2" || m_magic == "P1" ) )
    vw_throw( ArgumentErr() << "DiskImageResourcePBM: unsupported / or incorrect magic number identifer \"" << m_magic << "\". Possibly not PBM image." );

  // Getting image width, height, and max gray value.
  int32 iwidth, iheight;
  skip_any_comments(input);
  input >> iwidth;
  skip_any_comments(input);
  input >> iheight;

  if ( m_magic != "P1" && m_magic != "P4" ) {
    skip_any_comments(input);
    input >> m_max_value;
  } else
    m_max_value = 1; // Binary image

  // skip the one required whitespace
  if (!isspace(input.get()))
    vw_throw( IOErr() << "DiskImageResourcePBM: badly-formed file: " << filename);

  m_image_data_position = input.tellg();

  input.close();

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
  VW_ASSERT( dest.format.cols==uint32(cols()) && dest.format.rows==uint32(rows()),
             IOErr() << "Buffer has wrong dimensions in PBM read." );

  ifstream input(m_filename.c_str(), fstream::in|fstream::binary);

  if (!input.is_open())
    vw_throw( IOErr() << "DiskImageResourcePBM: Failed to open \""
              << m_filename << "\"." );
  input.seekg(m_image_data_position);

  // Reading image data
  ImageBuffer src(m_format, 0);
  size_t data_size = src.pstride * src.format.planes;
  scoped_array<uint8> image_data(new uint8[data_size]);

  // TODO: P4 is broken; binary bool is packed, and we're not doing that yet
  if ( m_magic == "P4" )
    vw_throw( NoImplErr() << "P4 (PBM Binary) is not currently implemented" );

  if ( m_magic == "P1" || m_magic == "P2" || m_magic == "P3" ) {
    // Bool/Grey/RGB (respectively) stored as ASCII
    uint32 buf;
    for ( size_t i = 0; i < data_size; ++i ) {
      input >> buf;
      image_data[i] = buf;
    }
  } else if ( m_magic == "P4" || m_magic == "P5" || m_magic == "P6" ) {
    // Bool/Grey/RGB (respectively) stored as Binary
    input.read( reinterpret_cast<char*>(image_data.get()), data_size );
  } else {
    vw_throw( NoImplErr() << "Unknown input channel type." );
  }

  if ( m_magic != "P1" && m_magic != "P4" )
    normalize( image_data, data_size, m_max_value );

  src.data = image_data.get();
  convert( dest, src, m_rescale );
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
    m_magic = "P4";
  } else if ( boost::algorithm::iends_with( filename, ".pgm" ) ||
              boost::algorithm::iends_with( filename, ".PGM" )  ) {
    m_magic = "P5";
  } else if ( boost::algorithm::iends_with( filename, ".ppm" ) ||
              boost::algorithm::iends_with( filename, ".PPM" ) ) {
    m_magic = "P6";
  } else {
    // Give up and guess based on image format
    int channels = num_channels(m_format.pixel_format);

    if (channels == 1) {
      if (m_format.channel_type == VW_CHANNEL_BOOL)
        m_magic = "P4";
      else
        m_magic = "P5";
    }
    else if (channels == 3) {
      m_magic = "P6";
    } else {
      vw_throw( NoImplErr() << "Unsupported number of channels: " << channels );
    }
  }

  // TODO: P4 is broken; binary bool is packed, and we're not doing that yet
  if ( m_magic == "P4" )
    vw_throw( NoImplErr() << "P4 (PBM Binary) is not currently implemented" );

  if (m_magic == "P4") {
    m_format.pixel_format = VW_PIXEL_SCALAR;
    m_format.channel_type = VW_CHANNEL_BOOL;
    m_max_value = 1;
  } else if (m_magic == "P5") {
    m_format.pixel_format = VW_PIXEL_GRAY;
    m_format.channel_type = VW_CHANNEL_UINT8;
    m_max_value = 255;
  } else if (m_magic == "P6") {
    m_format.pixel_format = VW_PIXEL_RGB;
    m_format.channel_type = VW_CHANNEL_UINT8;
    m_max_value = 255;
  }

  if (default_ascii) {
    if (m_magic == "P4")
      m_magic = "P1";
    else if (m_magic == "P5")
      m_magic = "P2";
    else if (m_magic == "P6")
      m_magic = "P3";
  }

  fstream output(filename.c_str(), fstream::out|fstream::binary);
  output.exceptions(ofstream::failbit | ofstream::badbit);

  output << m_magic       << "\n"
         << m_format.cols << " "
         << m_format.rows << "\n";

  if ( m_magic != "P1" && m_magic != "P4" )
    output << m_max_value << "\n";

  m_image_data_position = output.tellp();
}

// Write the given buffer into the disk image.
void DiskImageResourcePBM::write( ImageBuffer const& src,
                                  BBox2i const& bbox ) {
  VW_ASSERT( bbox.width()==int(cols()) && bbox.height()==int(rows()),
             NoImplErr() << "DiskImageResourcePBM does not support partial writes." );
  VW_ASSERT( src.format.cols==uint32(cols()) && src.format.rows==uint32(rows()),
             IOErr() << "Buffer has wrong dimensions in PBM write." );

  fstream output(m_filename.c_str(), fstream::out|fstream::binary|fstream::app);
  output.exceptions(ofstream::failbit | ofstream::badbit);
  output.seekp(m_image_data_position);

  ImageBuffer dst(m_format, 0);
  size_t data_size = dst.pstride * dst.format.planes;

  // TODO: P4 is broken; binary bool is packed, and we're not doing that yet
  if ( m_magic == "P4" )
    vw_throw( NoImplErr() << "P4 (PBM Binary) is not currently implemented" );

  scoped_array<uint8> image_data(new uint8[data_size]);
  dst.data = image_data.get();
  convert( dst, src, m_rescale );

  if ( m_magic == "P1" || m_magic == "P2" || m_magic == "P3" ) {
    // Bool/Grey/RGB (respectively) stored as ASCII
    if (data_size > 0)
      output << static_cast<uint32>(image_data[0]);

    for ( size_t i = 1; i < data_size; i++ )
      output << " " << static_cast<uint32>(image_data[i]);
  } else if ( m_magic == "P4" || m_magic == "P5" || m_magic == "P6" ) {
    // Bool/Grey/RGB (respectively) stored as binary
    output.write( reinterpret_cast<const char*>(image_data.get()), data_size );
  } else {
    vw_throw( NoImplErr() << "Unknown input channel type." );
  }
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
