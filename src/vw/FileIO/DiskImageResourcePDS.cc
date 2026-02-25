// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


/// \file DiskImageResourcePDS.cc
///
#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <vw/Core/Exception.h>
#include <vw/Core/Log.h>
#include <vw/FileIO/DiskImageResourcePDS.h>

#include <vector>
#include <string>
#include <cstring> // For memset()

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
using namespace boost;

static bool cpu_is_big_endian() {
#if VW_BYTE_ORDER == VW_BIG_ENDIAN
  return true;
#else
  return false;
#endif
}

void vw::DiskImageResourcePDS::treat_invalid_data_as_alpha() {
  // We currently only support this feature under very specific circumstances
  std::string format_str, sample_bits_str, valid_minimum_str;
  if( !( query("SAMPLE_TYPE",format_str) && query("SAMPLE_BITS",sample_bits_str) && query("VALID_MINIMUM",valid_minimum_str) ) )
    vw_throw( vw::NoImplErr() << "Invalid data not supported for this PDS image." );
  if( format_str != "MSB_INTEGER" || sample_bits_str != "16" || m_format.pixel_format != VW_PIXEL_GRAY )
    vw_throw( vw::NoImplErr() << "Invalid data not supported for this PDS image format." );
  m_invalid_as_alpha = true;
}

void vw::DiskImageResourcePDS::parse_pds_header(std::vector<std::string> const& header) {

  for ( unsigned i=0; i<header.size(); ++i ) {
    std::string line = header[i];

    // Locate lines that have a key/value pair (with key = value syntax)
    std::vector<std::string> split_vector;
    split( split_vector, line, is_any_of("=") );
    if ( split_vector.size() == 2 ) {
      trim_left(split_vector[0]);
      trim_right(split_vector[0]);
      trim_left(split_vector[1]);
      trim_right(split_vector[1]);
      m_header_entries[split_vector[0]] = split_vector[1];
    }
  }
}

vw::PixelFormatEnum vw::DiskImageResourcePDS::planes_to_pixel_format(int32 planes) const {
  switch( planes ) {
  case 1:  return VW_PIXEL_GRAY;   break;
  case 2:  return VW_PIXEL_GRAYA;  break;
  case 3:  return VW_PIXEL_RGB;    break;
  case 4:  return VW_PIXEL_RGBA;   break;
  }
  return VW_PIXEL_SCALAR;
}

/// Bind the resource to a file for reading.  Confirm that we can open
/// the file and that it has a sane pixel format.
void vw::DiskImageResourcePDS::open( std::string const& filename ) {

  FILE* input_file = fopen(filename.c_str(), "r");
  if( !input_file )
    vw_throw( vw::ArgumentErr() << "DiskImageResourcePDS: Failed to open \""
              << filename << "\"." );

  char c_line[2048];
  int i = 0;

  std::vector<std::string> header;

  // Read the entire header section and place it into a vector of
  // strings where each string is one line of the file.
  const static int MAX_PDS_HEADER_SIZE = 1000;
  while ( fgets(c_line, 2048, input_file) ) {
    i++;
    if ((c_line[0] == 'E' && c_line[1] == 'N' && c_line[2] == 'D' && c_line[3] != '_') ||
        (i > MAX_PDS_HEADER_SIZE))
      break;
    header.push_back(std::string(c_line));
  }
  fclose(input_file);

  // The the data into an associative contain (std::map).  Key/value
  // pairs are located by searching for strings separated by the
  // equals sign "=".
  parse_pds_header(header);

  std::vector<std::string> keys;
  std::string value;
  bool valid = true;

  // Query for the number of columns
  keys.push_back("LINE_SAMPLES");
  keys.push_back("NS");
  keys.push_back("/IMAGE/LINE_SAMPLES");
  valid = valid && query( keys, value );
  m_format.cols = atol(value.c_str());

  // Query for the number of rows
  keys.clear();
  keys.push_back("NL");
  keys.push_back("IMAGE_LINES");
  keys.push_back("LINES");
  keys.push_back("/IMAGE/LINES");

  valid = valid && query( keys, value );
  m_format.rows = atol(value.c_str());

  // Image planes (if there is no tag in the PDS for bands, we
  // assume one image plane)
  keys.clear();
  keys.push_back("BANDS");
  keys.push_back("/IMAGE/BANDS");
  if( query( keys, value ) ) {
    m_format.planes = atol(value.c_str());
  } else {
    m_format.planes = 1;
  }

  // Band storage type (if not specified, we assume interleaved samples)
  keys.clear();
  keys.push_back("BAND_STORAGE_TYPE");
  if ( query( keys, value ) ) {
    if ( value == "SAMPLE_INTERLEAVED" ) {
      m_band_storage = SAMPLE_INTERLEAVED;
    } else if ( value == "BAND_SEQUENTIAL" ) {
      m_band_storage = BAND_SEQUENTIAL;
    } else {
      vw_throw( NoImplErr() << "DiskImageResourcePDS: unsupported band storage type " << value );
    }
  } else {
    m_band_storage = SAMPLE_INTERLEAVED;
  }

  // Channel type
  keys.clear();
  keys.push_back("SAMPLE_TYPE");
  std::string format_str;
  valid = valid && query( keys, format_str );
  keys.clear();
  keys.push_back("SAMPLE_BITS");
  std::string sample_bits_str;
  valid = valid && query( keys, sample_bits_str );

  // Number of bytes in the PDS header (essentially the offset
  // before the image data begins.
  int record_size = 1;
  keys.clear();
  keys.push_back("RECORD_BYTES");
  keys.push_back("/RECSIZE");
  keys.push_back("HEADER_RECORD_BYTES");
  if( query( keys, value ) ) {
    record_size = atol(value.c_str());
  }

  keys.clear();
  keys.push_back("^IMAGE");
  if ( query( keys, value ) ) {
    // Check to see if the label contained a quoted string.  If so,
    // then we take this to mean that the actual image data is in a
    // separate file, and the filename is contained in this key.
    // Otherwise, we set the data filename to be the same as the
    // filename that contains the image header.
    if (value[0] == '\"' && value[value.size()-1] == '\"') {
      VW_OUT(InfoMessage, "fileio") << "PDS header points to a separate data file: " << value << ".\n";
      m_pds_data_filename = value.substr(1,value.size()-2);
      m_image_data_offset = 0;
    } else {
      m_pds_data_filename = DiskImageResource::m_filename;
      // m_image_data_offset = record_size * atol(value.c_str()) - 1;  // ????
      m_image_data_offset = record_size * (atol(value.c_str()) - 1);
    }
  } else {
    keys.clear();
    keys.push_back("LABEL_RECORDS");
    if( query( keys, value ) ) {
      m_image_data_offset = record_size * atol(value.c_str());
    }
    else m_image_data_offset = record_size;
  }

  if( ! valid ) {
    vw_throw( ArgumentErr() << "DiskImageResourcePDS: Could not find complete image header. Possibly not PDS image." );
  }

  m_file_is_msb_first = true;
  if (format_str == "UNSIGNED_INTEGER" ||
      format_str == "MSB_UNSIGNED_INTEGER" ||
      format_str == "LSB_UNSIGNED_INTEGER") {

    if (sample_bits_str == "8")
      m_format.channel_type = VW_CHANNEL_UINT8;
    else if (sample_bits_str == "16")
      m_format.channel_type = VW_CHANNEL_UINT16;

    if (format_str == "LSB_UNSIGNED_INTEGER")
      m_file_is_msb_first = false;

  } else if (format_str == "INTEGER" ||
             format_str == "MSB_INTEGER" ||
             format_str == "LSB_INTEGER") {

    if (sample_bits_str == "8")
      m_format.channel_type = VW_CHANNEL_INT8;
    else if (sample_bits_str == "16")
      m_format.channel_type = VW_CHANNEL_INT16;

    if (format_str == "LSB_INTEGER")
      m_file_is_msb_first = false;

  } else {
    vw_throw( IOErr() << "DiskImageResourcePDS: Unsupported pixel type in \"" << filename << "\"." );
  }

  // Match buffer format to band storage type
  m_format.pixel_format = planes_to_pixel_format(m_format.planes);
  if (m_format.pixel_format != VW_PIXEL_SCALAR) m_format.planes = 1;

  VW_OUT(DebugMessage, "fileio")
    << "Opening PDS Image\n"
    << "\tImage Dimensions: " << m_format.cols << "x" << m_format.rows << "x" << m_format.planes << "\n"
    << "\tImage Format: " << m_format.channel_type << "   " << m_format.pixel_format << "\n";
}

/// Bind the resource to a file for writing.
void vw::DiskImageResourcePDS::create( std::string const& /*filename*/,
                                       ImageFormat const& /*format*/ )
{
  vw_throw( NoImplErr() << "The PDS driver does not yet support creation of PDS files." );
}

/// Read the disk image into the given buffer.
void vw::DiskImageResourcePDS::read( ImageBuffer const& dest, BBox2i const& bbox ) const
{
  VW_ASSERT( bbox.width()==int(cols()) && bbox.height()==int(rows()),
             NoImplErr() << "DiskImageResourcePDS does not support partial reads." );
  VW_ASSERT( dest.format.cols==uint32(cols()) && dest.format.rows==uint32(rows()),
             IOErr() << "Buffer has wrong dimensions in PDS read." );

  // Re-open the file, and shift the file offset to the position of
  // the first image byte (as indicated by the PDS header).  Some PDS
  // files will have the actual data in a separate file that is
  // pointed to by the PDS image header, so we may actually be opening
  // that file here instead of the original file that the user
  // specified.
  //
  // NOTE: The filename encoded in the ^IMAGE tag seems to sometimes
  // differ in case from the actual data file, so we try a few
  // different combinations here.
  std::ifstream image_file(m_pds_data_filename.c_str(), std::ios::in | std::ios::binary);
  if (image_file.bad()) {
    image_file.open(boost::to_lower_copy(m_pds_data_filename).c_str(), std::ios::in | std::ios::binary);
    if (image_file.bad()) {
      image_file.open(boost::to_upper_copy(m_pds_data_filename).c_str(), std::ios::in | std::ios::binary);
      if (image_file.bad()) {
        vw_throw( vw::ArgumentErr() << "DiskImageResourcePDS: Failed to open \""
                  << DiskImageResource::m_filename << "\"." );
      }
    }
  }
  image_file.seekg(m_image_data_offset, std::ios::beg);

  // Grab the pixel data from the file.
  unsigned total_pixels = (unsigned)( m_format.cols * m_format.rows * m_format.planes );
  unsigned bytes_per_pixel = 1;
  if ( m_format.channel_type == VW_CHANNEL_UINT16 ||
       m_format.channel_type == VW_CHANNEL_INT16 ) {
    bytes_per_pixel = 2;
  }
  else if ( ! ( m_format.channel_type == VW_CHANNEL_UINT8 ||
                m_format.channel_type == VW_CHANNEL_INT8 ) ) {
    vw_throw( IOErr() << "DiskImageResourcePDS: Unsupported channel type (" << m_format.channel_type << ")." );
  }
  bytes_per_pixel *= num_channels(m_format.pixel_format);
  uint8* image_data = new uint8[total_pixels * bytes_per_pixel];
  image_file.read((char*)image_data, bytes_per_pixel*total_pixels);

  if (image_file.bad())
    vw_throw(IOErr() << "DiskImageResourcePDS: an unrecoverable error occurred while reading the image data.");

  // Convert the endian-ness of the data if the architecture of the
  // machine and the endianness of the file do not match.
  if (m_format.channel_type == VW_CHANNEL_INT16 ||
      m_format.channel_type == VW_CHANNEL_UINT16) {
    if ((cpu_is_big_endian() && !m_file_is_msb_first) ||
        (!cpu_is_big_endian() && m_file_is_msb_first) ) {
      for ( unsigned i=0; i<total_pixels*bytes_per_pixel; i+=2 ) {
        uint8 temp = image_data[i+1];
        image_data[i+1] = image_data[i];
        image_data[i] = temp;
      }
    }
  }

  // For band sequential images, we must copy the data over into
  // interleaved format.
  if ( m_band_storage == BAND_SEQUENTIAL && m_format.pixel_format != VW_PIXEL_SCALAR) {
    uint8* intermediate_data = new uint8[total_pixels * bytes_per_pixel];
    int n_channels = num_channels(m_format.pixel_format);
    int n_pixels = m_format.cols * m_format.rows;
    for (int n = 0; n < n_channels; ++n) {
      for (int p = 0; p < n_pixels; ++p) {
        intermediate_data[n_channels*p+n] = image_data[n_pixels*n+p];
      }
    }
    // Swap over to the new image buffer
    uint8* temporary_ptr = image_data;
    image_data = intermediate_data;
    delete[] temporary_ptr;
  }

  // set up an image buffer around the PDS data.
  ImageBuffer src;
  src.data = image_data;
  src.format = m_format;
  src.cstride = bytes_per_pixel;
  src.rstride = bytes_per_pixel * m_format.cols;
  src.pstride = bytes_per_pixel * m_format.cols * m_format.rows;
  convert( dest, src, m_rescale );

  if ( m_invalid_as_alpha ) {
    // We checked earlier that the source format is as we
    // expect.  Now we sanity-check the destination.
    if( dest.format.planes == 1 &&
        ( dest.format.pixel_format == VW_PIXEL_GRAYA ||
          dest.format.pixel_format == VW_PIXEL_RGBA ) ) {
      int dst_bpp = num_channels(dest.format.pixel_format) * channel_size(dest.format.channel_type);
      std::string valid_minimum_str;
      if ( query( "VALID_MINIMUM", valid_minimum_str ) ) {
        int16 valid_minimum = atoi(valid_minimum_str.c_str());
        uint8* src_row = (uint8*)src.data;
        uint8* dst_row = (uint8*)dest.data;
        for( uint32 y=0; y<m_format.rows; ++y ) {
          uint8* src_data = src_row;
          uint8* dst_data = dst_row;
          for( uint32 x=0; x<m_format.cols; ++x ) {
            if( *((int16*)src_data) < valid_minimum ) {
              std::memset( dst_data, 0, dst_bpp );
            }
            src_data += src.cstride;
            dst_data += dest.cstride;
          }
          src_row += src.rstride;
          dst_row += dest.rstride;
        }
      }
    }
  }

  delete[] image_data;
  image_file.close();
}

// Write the given buffer into the disk image.
void vw::DiskImageResourcePDS::write( ImageBuffer const& /*src*/, BBox2i const& /*bbox*/ )
{
  vw_throw( NoImplErr() << "The PDS driver does not yet support creation of PDS files." );
}

// A FileIO hook to open a file for reading
vw::DiskImageResource* vw::DiskImageResourcePDS::construct_open( std::string const& filename ) {
  return new DiskImageResourcePDS( filename );
}

// A FileIO hook to open a file for writing
vw::DiskImageResource* vw::DiskImageResourcePDS::construct_create( std::string const& filename,
                                                                   ImageFormat const& format ) {
  return new DiskImageResourcePDS( filename, format );
}
