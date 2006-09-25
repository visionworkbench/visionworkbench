// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
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

/// \file FileIO/DiskImageResourceTIFF.cc
/// 
/// Provides support for TIFF image files.
///

#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <vector>

#include <tiffio.h>

#include <vw/Core/Exception.h>
#include <vw/FileIO/DiskImageResourceTIFF.h>

// Define the following macro if you want libtiff to show warning messages
// #define VW_TIFF_SHOW_WARNINGS

void tiff_warning_handler(const char* module, const char* frmt, va_list ap) {
#ifdef VW_TIFF_SHOW_WARNINGS
  std::cout << "Warning  ";
  if (module) { printf(frmt, module); }
  else { printf(frmt, ""); }
  printf("\n");
#endif
}

void tiff_error_handler(const char* module, const char* frmt, va_list ap) {
  char error_msg[2048];
  if (module) { snprintf(error_msg, 2048, frmt, module); }
  else { snprintf(error_msg, 2048, frmt, ""); }
  throw vw::IOErr() << "DiskImageResourceTIFF: A libtiff error occured. " << std::string(error_msg);
}


/// Close the TIFF file when the object is destroyed
vw::DiskImageResourceTIFF::~DiskImageResourceTIFF() {
  this->flush();
}
  
/// Flush the buffered data to disk
void vw::DiskImageResourceTIFF::flush() {
  if (m_tif_ptr) {
    TIFFFlush((TIFF*) m_tif_ptr);
    m_tif_ptr = NULL;
  }
}

/// Bind the resource to a file for reading.  Confirm that we can open
/// the file and that it has a sane pixel format.  
void vw::DiskImageResourceTIFF::open( std::string const& filename )
{
  TIFFSetWarningHandler(&tiff_warning_handler);
  TIFFSetErrorHandler(&tiff_error_handler);

  TIFF* tif = TIFFOpen(filename.c_str(), "r");
  if( !tif  ) throw vw::IOErr() << "Failed to open \"" << filename << "\" using libTIFF.";
  m_tif_ptr = (void*) tif;
  m_filename = filename;

  TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &(m_format.cols));
  TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &(m_format.rows));
  TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &(m_format.planes));

  uint32 sample_format = 0, bits_per_sample = 0;
  TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bits_per_sample);
  TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT, &sample_format);
  // Some TIFF files don't actually specify the sample format, but in
  // this case the format is almost always UINT, so we force this
  // assumption here.
  if (sample_format == 0) {
    sample_format = SAMPLEFORMAT_UINT;
  }

  switch (sample_format) {
  case SAMPLEFORMAT_UINT:
    if (bits_per_sample == 8)
      m_format.channel_type = VW_CHANNEL_UINT8;
    else if (bits_per_sample == 16) 
      m_format.channel_type = VW_CHANNEL_UINT16;
    else if (bits_per_sample == 32) 
      m_format.channel_type = VW_CHANNEL_UINT32;
    else if (bits_per_sample == 64) 
      m_format.channel_type = VW_CHANNEL_UINT64;
    else 
      throw IOErr() << "DiskImageResourceTIFF: Unsupported pixel format.";
    break;
  case SAMPLEFORMAT_INT:
    if (bits_per_sample == 8)
      m_format.channel_type = VW_CHANNEL_INT8;
    else if (bits_per_sample == 16) 
      m_format.channel_type = VW_CHANNEL_INT16;
    else if (bits_per_sample == 32) 
      m_format.channel_type = VW_CHANNEL_INT32;
    else if (bits_per_sample == 64) 
      m_format.channel_type = VW_CHANNEL_INT64;
    else 
      throw IOErr() << "DiskImageResourceTIFF: Unsupported pixel format.";
    break;
  case SAMPLEFORMAT_IEEEFP:
    if (bits_per_sample == 16) 
      m_format.channel_type = VW_CHANNEL_FLOAT16;
    else if (bits_per_sample == 32) 
      m_format.channel_type = VW_CHANNEL_FLOAT32;
    else if (bits_per_sample == 64) 
      m_format.channel_type = VW_CHANNEL_FLOAT64;    
    else 
      throw IOErr() << "DiskImageResourceTIFF: Unsupported pixel format.";
    break;
  default:
    throw IOErr() << "DiskImageResourceTIFF: Unsupported pixel format.";
  }
  
  uint32 plane_configuration = 0;
  TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &plane_configuration);

  // XXX TODO: Tiff might actually provide us with some info on
  // colorimetric interpretation of the channels, so maybe we should
  // try to use that here as well?
  if (plane_configuration == PLANARCONFIG_CONTIG) {
    switch( m_format.planes ) {
    case 1:  m_format.pixel_format = VW_PIXEL_GRAY;   break;
    case 2:  m_format.pixel_format = VW_PIXEL_GRAYA;  m_format.planes=1; break;
    case 3:  m_format.pixel_format = VW_PIXEL_RGB;    m_format.planes=1; break;
    case 4:  m_format.pixel_format = VW_PIXEL_RGBA;   m_format.planes=1; break;
    default: m_format.pixel_format = VW_PIXEL_SCALAR; break;
    }
  } else { 
    m_format.pixel_format = VW_PIXEL_SCALAR;
  }
}

/// Bind the resource to a file for writing.
void vw::DiskImageResourceTIFF::create( std::string const& filename, 
                                        GenericImageFormat const& format )
{
  if( format.planes!=1 && format.pixel_format!=VW_PIXEL_SCALAR )
    throw NoImplErr() << "TIFF doesn't support multi-plane images with compound pixel types.";

  // Set the TIFF warning and error handlers to Vision Workbench
  // functions, so that we can handle them ourselves.
  TIFFSetWarningHandler(&tiff_warning_handler);
  TIFFSetErrorHandler(&tiff_error_handler);

  TIFF* tif = TIFFOpen(filename.c_str(), "w");
  if( !tif  ) throw vw::IOErr() << "Failed to create \"" << filename << "\" using libTIFF.";
  m_tif_ptr = (void*) tif;

  m_filename = filename;
  m_format = format;
}

/// Read the disk image into the given buffer.
void vw::DiskImageResourceTIFF::read( GenericImageBuffer const& dest ) const
{
  VW_ASSERT( dest.format.cols==cols() && dest.format.rows==rows(),
             IOErr() << "Buffer has wrong dimensions in TIFF read." );

  TIFF* tif = (TIFF*) m_tif_ptr;
  uint32 image_length = 0, image_width = 0, config = 0;

  TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &image_length);
  TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &image_width);
  TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &config);
  uint32 scanline_size = TIFFScanlineSize(tif);

  // Allocate a buffer for reading in the data and read it, scanline
  // by scanline from disk.
  tdata_t buf = _TIFFmalloc(scanline_size * image_length * m_format.planes);
  uint16 s, nsamples;
  
  for (uint32 p = 0; p < m_format.planes; p++) {
    for (uint32 row = 0; row < image_length; row++) {
      TIFFReadScanline(tif, (uint8*)buf + row*scanline_size+p*scanline_size*image_length, row, p);
    }
  }

  // Set up a generic image buffer around the tdata_t buf object that
  // tiff used to copy it's data.
  GenericImageBuffer src;
  src.data = (uint8*)buf;
  src.format = m_format;
  src.cstride = scanline_size / image_width;
  src.rstride = scanline_size;
  src.pstride = scanline_size * image_length;
  
  convert( dest, src );
  _TIFFfree(buf);
}

// Write the given buffer into the disk image.
void vw::DiskImageResourceTIFF::write( GenericImageBuffer const& src )
{
  VW_ASSERT( src.format.cols==cols() && src.format.rows==rows(),
             IOErr() << "Buffer has wrong dimensions in TIFF write." );

  TIFF* tif = (TIFF*) m_tif_ptr;

  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, m_format.cols);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH, m_format.rows);
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8*channel_size(m_format.channel_type));

  if (m_format.pixel_format == VW_PIXEL_RGB ||
      m_format.pixel_format == VW_PIXEL_RGBA) {
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  } else {
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
  }

  switch (m_format.channel_type) {
  case VW_CHANNEL_INT8:
  case VW_CHANNEL_INT16:
  case VW_CHANNEL_INT32:
  case VW_CHANNEL_INT64:
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_INT);
    break;
  case VW_CHANNEL_UINT8:
  case VW_CHANNEL_UINT16:
  case VW_CHANNEL_UINT32:
  case VW_CHANNEL_UINT64:
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
    break;
  case VW_CHANNEL_FLOAT16:
  case VW_CHANNEL_FLOAT32:
  case VW_CHANNEL_FLOAT64:
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    break;
  default:
    throw IOErr() << "DiskImageResourceTIFF: Unsupported VW channel type.";
  }

  // Allocate some buffer memory for the output data
  uint32 scanline_size = num_channels(m_format.pixel_format) * channel_size(m_format.channel_type) * m_format.cols;
  tdata_t buf = _TIFFmalloc(scanline_size * m_format.rows * m_format.planes);

  // Set up the generic image buffer and convert the data into this buffer.
  GenericImageBuffer dst;
  dst.data = (uint8*)buf;
  dst.format = m_format;
  dst.cstride = num_channels(m_format.pixel_format) * channel_size(m_format.channel_type);
  dst.rstride = dst.cstride * m_format.cols;
  dst.pstride = dst.rstride * m_format.rows;
  convert( dst, src );

  if (m_format.pixel_format == VW_PIXEL_SCALAR) {
    // Multi-plane images with simple pixel types are stored in seperate
    // planes in the TIFF image.
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_SEPARATE);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, m_format.planes);
  } else {  
    // Compound pixel types are stored contiguously in TIFF files
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, num_channels(m_format.pixel_format));
  }

  // Write the image data to disk.
  for (uint32 p = 0; p < m_format.planes; p++) {
    for (uint32 row = 0; row < m_format.rows; row++) {
      TIFFWriteScanline(tif, (uint8*)buf + row*scanline_size + p * scanline_size * m_format.rows, row);
    }
  }

  // Clean up
  _TIFFfree(buf);
}

// A FileIO hook to open a file for reading
vw::DiskImageResource* vw::DiskImageResourceTIFF::construct_open( std::string const& filename ) {
  return new DiskImageResourceTIFF( filename );
}

// A FileIO hook to open a file for writing
vw::DiskImageResource* vw::DiskImageResourceTIFF::construct_create( std::string const& filename,
                                                                    GenericImageFormat const& format ) {
  return new DiskImageResourceTIFF( filename, format );
}
