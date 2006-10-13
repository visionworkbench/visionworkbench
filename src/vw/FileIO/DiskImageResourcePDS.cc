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

/// \file DiskImageResourcePDS.cc
/// 
#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>   
using namespace std;

#include <boost/algorithm/string.hpp>
using namespace boost;

#include <vw/Core/Exception.h>
#include <vw/FileIO/DiskImageResourcePDS.h>


/// Close the file when the object is destroyed
vw::DiskImageResourcePDS::~DiskImageResourcePDS() {
  this->flush();
}
 
/// Flush the buffered data to disk
void vw::DiskImageResourcePDS::flush() {
}

void vw::DiskImageResourcePDS::parse_pds_header(std::vector<std::string> const& header) {

  for (int i = 0; i < header.size(); i++) {
    string line = header[i];

    // Locate lines that have a key/value pair (with key = value syntax)
    vector<string> split_vector;
    split( split_vector, line, is_any_of("=") );
    if (split_vector.size() == 2) {
      trim_left(split_vector[0]);
      trim_right(split_vector[0]);
      trim_left(split_vector[1]);
      trim_right(split_vector[1]);
      m_header_entries[split_vector[0]] = split_vector[1];
    } 
  }
}

/// Bind the resource to a file for reading.  Confirm that we can open
/// the file and that it has a sane pixel format.  
void vw::DiskImageResourcePDS::open( std::string const& filename ) {

  FILE* input_file = fopen(filename.c_str(), "r");
  if( ! input_file ) throw vw::IOErr() << "Failed to open \"" << filename << "\".";
  m_filename = filename;

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
  // pairs are located by searching for strings seperated by the
  // equals sign "=".
  parse_pds_header(header);

  try {
    vector<string> keys;

    // Query for the number of columns
    keys.push_back("LINE_SAMPLES");
    keys.push_back("NS");
    keys.push_back("/IMAGE/LINE_SAMPLES");
    m_format.cols = atol(query(keys).c_str());

    // Query for the number of rows
    keys.clear();
    keys.push_back("NL");
    keys.push_back("IMAGE_LINES");
    keys.push_back("LINES");
    keys.push_back("/IMAGE/LINES");
    m_format.rows = atol(query(keys).c_str());

    // Image planes (if there is no tag in the PDS for bands, we
    // assume one image plane)
    try {
      keys.clear();
      keys.push_back("BANDS");
      keys.push_back("/IMAGE/BANDS");
      m_format.planes = atol(query(keys).c_str());
    } catch (NotFoundErr &e) {
      m_format.planes = 1;
    }

    // pixel type
    keys.clear();
    keys.push_back("SAMPLE_TYPE");
    string format_str = query(keys);
    keys.clear();
    keys.push_back("SAMPLE_BITS");
    string sample_bits_str = query(keys);
    
    if (format_str == "UNSIGNED_INTEGER") {
     
      if (sample_bits_str == "8") 
        m_format.channel_type = VW_CHANNEL_UINT8; 
      else if (sample_bits_str == "16") 
        m_format.channel_type = VW_CHANNEL_UINT16; 

    } else if (format_str == "MSB_INTEGER") {

      if (sample_bits_str == "8") 
        m_format.channel_type = VW_CHANNEL_INT8; 
      else if (sample_bits_str == "16") 
        m_format.channel_type = VW_CHANNEL_INT16; 

    } else {
     throw IOErr() << "DiskImageResourcePDS: Unsupported pixel type in \"" << filename << "\".";
    }

    // Number of bytes in the PDS header (essentially the offset
    // before the image data begins.
    keys.clear();
    keys.push_back("RECORD_BYTES");
    keys.push_back("/RECSIZE");
    keys.push_back("HEADER_RECORD_BYTES");
    m_image_data_offset = atol(query(keys).c_str());
    
    try {
      keys.clear();
      keys.push_back("LABEL_RECORDS");
      m_image_data_offset *= atol(query(keys).c_str());
    } catch (NotFoundErr &e) {}
        
  } catch (vw::NotFoundErr &e) {
    throw IOErr() << "DiskImageResourcePDS: could not find critical information in the image header.\n";
  }
  
  switch( m_format.planes ) {
  case 1:  m_format.pixel_format = VW_PIXEL_GRAY;   break;
  case 2:  m_format.pixel_format = VW_PIXEL_GRAYA;  m_format.planes=1; break;
  case 3:  m_format.pixel_format = VW_PIXEL_RGB;    m_format.planes=1; break;
  case 4:  m_format.pixel_format = VW_PIXEL_RGBA;   m_format.planes=1; break;
  default: m_format.pixel_format = VW_PIXEL_SCALAR; break;
  }

  // For debugging:
  //   std::cout << "\tImage Dimensions: " << m_format.cols << "x" << m_format.rows << "x" << m_format.planes << "\n";
  //   std::cout << "\tImage Format: " << m_format.channel_type << "   " << m_format.pixel_format << "\n";
  fclose(input_file);
}

/// Bind the resource to a file for writing.
void vw::DiskImageResourcePDS::create( std::string const& filename, 
                                       GenericImageFormat const& format )
{
  throw NoImplErr() << "The PDS driver does not yet support creation of PDS files";
}

/// Read the disk image into the given buffer.
void vw::DiskImageResourcePDS::read( GenericImageBuffer const& dest ) const 
{
  VW_ASSERT( dest.format.cols==cols() && dest.format.rows==rows(),
             IOErr() << "Buffer has wrong dimensions in PDS read." );
  
  // Re-open the file, and shift the file offset to the position of
  // the first image byte (as indicated by the PDS header)
  FILE* input_file = fopen(m_filename.c_str(), "r");
  if( ! input_file ) throw vw::IOErr() << "Failed to open \"" << m_filename << "\".";
  fseek(input_file, m_image_data_offset, 0);

  // Grab the pixel data from the file.
  unsigned int total_pixels = m_format.cols * m_format.rows * m_format.planes;
  unsigned int bytes_per_pixel;
  if (m_format.channel_type == VW_CHANNEL_UINT8 ||
      m_format.channel_type == VW_CHANNEL_INT8) {
    bytes_per_pixel = 1;
  } else if (m_format.channel_type == VW_CHANNEL_UINT16 ||
             m_format.channel_type == VW_CHANNEL_INT16) {
    bytes_per_pixel = 2;
  }
  uint8* image_data = new uint8[total_pixels * bytes_per_pixel];
  unsigned bytes_read = fread(image_data, bytes_per_pixel, total_pixels, input_file);
  if (bytes_read != total_pixels) 
    throw IOErr() << "DiskImageResourcePDS: An error occured while reading the image data.";

  // Convert the endian-ness of the data
  if (m_format.channel_type == VW_CHANNEL_INT16 ||
      m_format.channel_type == VW_CHANNEL_UINT16) {
    for (int i = 0; i < total_pixels * bytes_per_pixel; i+=2) {
      uint8 temp = image_data[i+1];
      image_data[i+1] = image_data[i];
      image_data[i] = temp;
    }
  }
  
  uint8 max = 0, min = 255;
  for (int i = 0; i < total_pixels; i++) {
    max = image_data[i] > max ? image_data[i] : max;
    min = image_data[i] < min ? image_data[i] : min;
  }

  // set up a generic image buffer around the PDS data.
  GenericImageBuffer src;
  src.data = image_data;
  src.format = m_format;
  src.cstride = bytes_per_pixel;
  src.rstride = bytes_per_pixel * m_format.cols;
  src.pstride = bytes_per_pixel * m_format.cols * m_format.rows;
  
  convert( dest, src );

  delete[] image_data;
  fclose(input_file);
}

// Write the given buffer into the disk image.
void vw::DiskImageResourcePDS::write( GenericImageBuffer const& src ) 
{
  throw NoImplErr() << "The PDS driver does not yet support creation of PDS files";
}

// A FileIO hook to open a file for reading
vw::DiskImageResource* vw::DiskImageResourcePDS::construct_open( std::string const& filename ) {
  return new DiskImageResourcePDS( filename );
}

// A FileIO hook to open a file for writing
vw::DiskImageResource* vw::DiskImageResourcePDS::construct_create( std::string const& filename,
                                                                    GenericImageFormat const& format ) {
  return new DiskImageResourcePDS( filename, format );
}
