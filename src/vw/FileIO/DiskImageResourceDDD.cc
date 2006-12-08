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

/// \file DiskImageResourceDDD.cc
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
#include <sstream>   

#include <values.h>			   // for BITSPERBYTE

#include <boost/algorithm/string.hpp>

#include <vw/Core/Exception.h>
#include <vw/FileIO/DiskImageResourceDDD.h>

using namespace std;
using namespace boost;

namespace vw
{
  static const int IMAGE_HEADER_LENGTH = 1024;
  static const int IMAGE_LABEL_OFFSET = 24;
  static const int IMAGE_LABEL_LENGTH = IMAGE_HEADER_LENGTH-IMAGE_LABEL_OFFSET;
  static const int MAGIC = 1659;

  struct DDDHeader
  {
    long magic;				   // always set to MAGIC
    long numScanlines, bytesPerScanline;
    long bitsPerElement; // the data size unless there's padding
    long spare1, spare2;
    char label[IMAGE_LABEL_LENGTH];
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // The label contains linefeed delimited strings, e.g.:
  //
  // decompressed-from ./4A_04_1001004F00_01.DAT
  // id 79 time 844381966:250
  // start 0 cross 5056 down 7168
  // cam ctx
  // mode 0x0
  // dac 195
  // offset 234 232
  // sram_base 0
  // start_addr 0
  // exposure 1.87 msec (19)
  // FPA-temp 17.9 C (2395)
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /// Close the file when the object is destroyed
  DiskImageResourceDDD::~DiskImageResourceDDD()
  {
    this->flush();
  }
 
  /// Flush the buffered data to disk
  void
  DiskImageResourceDDD::flush()
  {
  }

  std::string
  DiskImageResourceDDD::query(std::string const& key) const
  {
    std::map<std::string, std::string>::const_iterator entry =
      m_header_entries.find(key);

    if (entry != m_header_entries.end()) 
      return (*entry).second;
    else 
      throw NotFoundErr()
	<< "DiskImageResourceDDD::query(): no matching value found for \""
	<< key << "\""; 
  }

  void
  DiskImageResourceDDD::parse_ddd_header(const char* header_label)
  {
    string header_label_string(header_label);
    istringstream label_stream(header_label_string);

    while (!label_stream.eof())
    {
      char line[IMAGE_LABEL_LENGTH+1];
      vector<string> tokens;

      label_stream.getline(line, sizeof(line));
      split(tokens, line, is_any_of(" "));

      // First handle any special cases... yes this is a bit hacky,
      // but their label format is a bit hacky...
      if (tokens[0] == string("id"))
      {
	m_header_entries[tokens[2]] = tokens[3];
      }
      else if (tokens[0] == string("start"))
      {
	m_header_entries[tokens[2]] = tokens[3];
	m_header_entries[tokens[4]] = tokens[5];
      }
      else if (tokens[0] == string("offset"))
      {
	tokens[1] = tokens[1] + tokens[2];
      }
      else if ((tokens[0] == string("exposure")) ||
	       (tokens[0] == string("FPA-temp")))
      {
	tokens[1] = tokens[1] + tokens[2] + tokens[3];
      }

      m_header_entries[tokens[0]] = tokens[1];
    }
  }

  inline long
  swap_long(long data)
  {
    unsigned char *bytes = reinterpret_cast<unsigned char *>(&data);
    unsigned char bytes2 = bytes[2];
    unsigned char bytes3 = bytes[3];

    bytes[3] = bytes[0];
    bytes[2] = bytes[1];
    bytes[1] = bytes2;
    bytes[0] = bytes3;

    return data;
  }

  inline void
  swap_longs(void *buffer, size_t numLongs)
  {
    unsigned long *longs = reinterpret_cast<unsigned long *>(buffer);
    unsigned long *endLong = longs + numLongs;

    while (longs != endLong)
    {
      unsigned char *bytes = (unsigned char *)(longs);
      unsigned char bytes2 = bytes[2];
      unsigned char bytes3 = bytes[3];

      bytes[3] = bytes[0];
      bytes[2] = bytes[1];
      bytes[1] = bytes2;
      bytes[0] = bytes3;

      longs++;
    }
  }

  inline void
  swap_shorts(void *buffer, size_t numShorts)
  {
    unsigned short *shorts = reinterpret_cast<unsigned short *>(buffer);
    unsigned short *endShort = shorts + numShorts;

    while (shorts != endShort)
    {
      unsigned char *bytes = (unsigned char *)(shorts);
      unsigned char bytes1 = bytes[1];
      bytes[1] = bytes[0];
      bytes[0] = bytes1;
      shorts++;
    }
  }

  /// Bind the resource to a file for reading.  Confirm that we can open
  /// the file and that it has a sane pixel format.  
  void
  DiskImageResourceDDD::open(std::string const& filename)
  {
    ifstream image_file(filename.c_str(), ios::in | ios::binary);

    if (!image_file)
      throw IOErr() << "DiskImageResourceDDD::open(): could not open \""
		    << filename << "\".";

    DDDHeader header;

    m_filename = filename;

    image_file.read((char *)(&header), sizeof(DDDHeader));

    if (image_file.bad())
      throw IOErr() << "DiskImageResourceDDD::open(): could not read "
		    << filename << " header.";
      
    if (header.magic == MAGIC)
      m_is_other_endian = false;
    else if (header.magic == swap_long(MAGIC))
      m_is_other_endian = true;
    else
      throw IOErr() << "DiskImageResourceDDD::open(): " << filename
		    << " has bad magic number (" << header.magic << " != "
		    << MAGIC << ").";

    if (m_is_other_endian)
      swap_longs(&header, 4);

    int bytes_per_pixel = header.bitsPerElement / BITSPERBYTE;
    m_format.planes = 1;
    m_format.pixel_format = VW_PIXEL_GRAY;
    m_format.cols = header.bytesPerScanline / bytes_per_pixel;
    m_format.rows = header.numScanlines;
 
   switch (header.bitsPerElement)
    {
    case BITSPERBYTE:
      m_format.channel_type = VW_CHANNEL_UINT8;
      break;
    case sizeof(short) * BITSPERBYTE:
      m_format.channel_type = VW_CHANNEL_INT16;
      break;
    default:
      throw IOErr() << "DiskImageResourceDDD::open(): unsupported pixel size ("
		    << header.bitsPerElement << " bits) in " << filename
		    << ".";
    }

    // Put the data into an associative contain (std::map).  Key/value
    // pairs are located by searching for strings seperated by the
    // equals sign "=".
    parse_ddd_header(header.label);

    image_file.close();
  }

  /// Bind the resource to a file for writing.
  void
  DiskImageResourceDDD::create(std::string const& filename,
			       GenericImageFormat const& format)
  {
    throw NoImplErr()
      << "The DDD driver does not yet support creation of DDD files";
  }

  /// Read the disk image into the given buffer.
  void
  DiskImageResourceDDD::read_generic(GenericImageBuffer const& dest) const
  {
    VW_ASSERT((dest.format.cols == cols()) && (dest.format.rows == rows()),
	      IOErr() << "Buffer has wrong dimensions in DDD read.");
  
    ifstream image_file(m_filename.c_str(), ios::in | ios::binary);

    if (!image_file)
      throw IOErr()
	<< "  DiskImageResourceDDD::read_generic(): failed to open \""
	<< m_filename << "\".";

    // Set the file offset to the position of the first image byte
    // (read in from the previously opened DDD header)
    image_file.seekg(m_image_data_offset, ios::beg);

    // Read the pixel data from the file.
    unsigned int total_pixels = (m_format.cols * m_format.rows *
				 m_format.planes);
    unsigned int bytes_per_pixel;
    switch (m_format.channel_type)
    {
    case VW_CHANNEL_UINT8: case VW_CHANNEL_INT8:
      bytes_per_pixel = sizeof(char);
      break;
    case VW_CHANNEL_UINT16: case VW_CHANNEL_INT16:
      bytes_per_pixel = sizeof(short);
      break;
    }
    uint8* image_data = new uint8[total_pixels * bytes_per_pixel];
    image_file.read((char *) image_data, bytes_per_pixel * total_pixels);

    if (image_file.bad())
      throw IOErr() << "DiskImageResourceDDD::read_generic():"
	" An unrecoverable error occured while reading the image data.";

    // DDD images are always big-endian, swab bytes if this is a
    // little-endian system
    if (m_is_other_endian)
    {
      switch (bytes_per_pixel)
      {
      case sizeof(short):
	std::cout << "DiskImageResourceDDD::read_generic(): "
		  << "swapping bytes (shorts)... " << std::flush;
	swap_shorts(image_data, m_format.cols * m_format.rows);
	std::cout << "DiskImageResourceDDD::read_generic(): done."
		  << std::endl;
	break;
      }
    }
  
    // Create generic image buffer from the DDD data.
    GenericImageBuffer src;
    src.data = image_data;
    src.format = m_format;
    src.cstride = bytes_per_pixel;
    src.rstride = bytes_per_pixel * m_format.cols;
    src.pstride = bytes_per_pixel * m_format.cols * m_format.rows;
  
    convert(dest, src);

    delete[] image_data;

    image_file.close();
  }

  // Write the given buffer into the disk image.
  void
  DiskImageResourceDDD::write_generic(GenericImageBuffer const& src)
  {
    throw NoImplErr()
      << "The DDD driver does not yet support creation of DDD files";
  }

  // A FileIO hook to open a file for reading
  DiskImageResource*
  DiskImageResourceDDD::construct_open(std::string const& filename)
  {
    return new DiskImageResourceDDD(filename);
  }

  // A FileIO hook to open a file for writing
  DiskImageResource*
  DiskImageResourceDDD::construct_create(std::string const& filename,
					 GenericImageFormat const& format)
  {
    return new DiskImageResourceDDD(filename, format);
  }
}
