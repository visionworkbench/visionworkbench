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
#include <vw/Core/Debugging.h>
#include <vw/FileIO/DiskImageResourceTIFF.h>

#ifndef VW_ERROR_BUFFER_SIZE
#define VW_ERROR_BUFFER_SIZE 2048
#endif

/// Handle libTIFF warning conditions by outputting message text at the 
/// DebugMessage verbosity level.
static void tiff_warning_handler(const char* module, const char* frmt, va_list ap) {
  char msg[VW_ERROR_BUFFER_SIZE];
  vsnprintf( msg, VW_ERROR_BUFFER_SIZE, frmt, ap );
  vw::vw_out(vw::DebugMessage) << "DiskImageResourceTIFF (" << (module?module:"none") << ") Warning: " << msg << std::endl;
}

/// Handle libTIFF error conditions by throwing an IOErr with the 
/// message text.
static void tiff_error_handler(const char* module, const char* frmt, va_list ap) {
  char msg[VW_ERROR_BUFFER_SIZE];
  vsnprintf( msg, VW_ERROR_BUFFER_SIZE, frmt, ap );
  throw vw::IOErr() << "DiskImageResourceTIFF (" << (module?module:"none") << ") Error: " << msg;
}

/// Destructor: flush and close the file.
vw::DiskImageResourceTIFF::~DiskImageResourceTIFF() {
  flush();
  if( m_tif_ptr ) {
    TIFFClose( (TIFF*)m_tif_ptr );
  }
  m_tif_ptr = 0;
}
  
/// Flush any buffered data to disk.
void vw::DiskImageResourceTIFF::flush() {
  if( m_tif_ptr ) {
    TIFFFlush((TIFF*) m_tif_ptr);
  }
}

vw::Vector2i vw::DiskImageResourceTIFF::native_read_block_size() const {
  TIFF *tif = (TIFF*)m_tif_ptr;
  if( TIFFIsTiled(tif) ) {
    vw_out(WarningMessage) << "DiskImageResourceTIFF Warning: Tile-based partial image reading is untested!" << std::cout;
    uint32 tile_width, tile_length;
    TIFFGetField( tif, TIFFTAG_TILEWIDTH, &tile_width );
    TIFFGetField( tif, TIFFTAG_TILELENGTH, &tile_length );
    return Vector2i(tile_width,tile_length);
  }
  else {
    uint32 rows_per_strip;
    TIFFGetField( tif, TIFFTAG_ROWSPERSTRIP, &rows_per_strip );
    return Vector2i(cols(),rows_per_strip);
  }
}

/// Bind the resource to a file for reading.  Confirm that we can open
/// the file and that it has a sane pixel format.  
void vw::DiskImageResourceTIFF::open( std::string const& filename ) {
  TIFFSetWarningHandler( &tiff_warning_handler );
  TIFFSetErrorHandler( &tiff_error_handler );

  TIFF* tif = TIFFOpen( filename.c_str(), "r" );
  if( !tif ) throw vw::IOErr() << "DiskImageResourceTIFF: Failed to open \"" << filename << "\" for reading!";
  m_filename = filename;
  m_tif_ptr = (void*) tif;

  TIFFGetField( tif, TIFFTAG_IMAGEWIDTH, &(m_format.cols) );
  TIFFGetField( tif, TIFFTAG_IMAGELENGTH, &(m_format.rows) );
  TIFFGetField( tif, TIFFTAG_SAMPLESPERPIXEL, &(m_format.planes) );

  uint32 sample_format = 0, bits_per_sample = 0;
  TIFFGetField( tif, TIFFTAG_BITSPERSAMPLE, &bits_per_sample );
  TIFFGetField( tif, TIFFTAG_SAMPLEFORMAT, &sample_format );
  // Some TIFF files don't actually specify the sample format, but in
  // this case the format is almost always UINT, so we force this
  // assumption here.
  if( sample_format == 0 ) {
    vw_out(DebugMessage) << "DiskImageResourceTIFF: " << m_filename << " does not specify sample format; assuming UINT!" << std::endl;
    sample_format = SAMPLEFORMAT_UINT;
  }

  m_format.channel_type = VW_CHANNEL_UNKNOWN;
  switch( sample_format ) {
  case SAMPLEFORMAT_UINT:
    if (bits_per_sample == 8)
      m_format.channel_type = VW_CHANNEL_UINT8;
    else if (bits_per_sample == 16) 
      m_format.channel_type = VW_CHANNEL_UINT16;
    else if (bits_per_sample == 32) 
      m_format.channel_type = VW_CHANNEL_UINT32;
    else if (bits_per_sample == 64) 
      m_format.channel_type = VW_CHANNEL_UINT64;
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
    break;
  case SAMPLEFORMAT_IEEEFP:
    if (bits_per_sample == 16) 
      m_format.channel_type = VW_CHANNEL_FLOAT16;
    else if (bits_per_sample == 32) 
      m_format.channel_type = VW_CHANNEL_FLOAT32;
    else if (bits_per_sample == 64) 
      m_format.channel_type = VW_CHANNEL_FLOAT64;    
    break;
  }
  if( ! m_format.channel_type ) {
    throw IOErr() << "DiskImageResourceTIFF: " << m_filename << " has an unsupported channel type ("
                  << sample_format << "," << bits_per_sample << ")!";
  }

  uint32 plane_configuration = 0;
  TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &plane_configuration);

  // XXX TODO: Tiff might actually provide us with some info on
  // colorimetric interpretation of the channels, so maybe we should
  // try to use that here as well?
  if( plane_configuration == PLANARCONFIG_CONTIG ) {
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
void vw::DiskImageResourceTIFF::read( GenericImageBuffer const& dest, BBox2i bbox ) const
{
  VW_ASSERT( dest.format.cols==bbox.width() && dest.format.rows==bbox.height(),
             ArgumentErr() << "DiskImageResourceTIFF (read) Error: Destination buffer has wrong dimensions!" );

  TIFF* tif = (TIFF*) m_tif_ptr;

  if( TIFFIsTiled(tif) ) {
    throw NoImplErr() << "DiskImageResourceTIFF (read) Error: Reading from tile-based TIFF files is not yet supported!";
  }
  else {
    tdata_t buf = _TIFFmalloc( TIFFStripSize(tif) );
    
    uint32 config = 0;
    TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &config);
    if( config==PLANARCONFIG_SEPARATE )
      throw NoImplErr() << "DiskImageResourceTIFF (read) Error: Reading from separate-plane TIFF files is not yet supported!";

    uint32 rows_per_strip;
    uint32 scanline_size = TIFFScanlineSize(tif);    
    TIFFGetField( tif, TIFFTAG_ROWSPERSTRIP, &rows_per_strip );

    GenericImageBuffer strip_src, strip_dest=dest;
    strip_src.format = m_format;
    strip_src.format.cols = bbox.width();
    strip_src.cstride = scanline_size / m_format.cols;
    strip_src.rstride = scanline_size;
    strip_src.pstride = scanline_size * m_format.rows;

    int minstrip=bbox.min().y()/rows_per_strip, maxstrip=(bbox.max().y()-1)/rows_per_strip;
    for( int strip=minstrip; strip<=maxstrip; ++strip ) {
      int strip_top = std::max(strip*rows_per_strip,uint32(bbox.min().y()));
      int strip_rows = std::min((strip+1)*rows_per_strip,uint32(bbox.max().y()))-strip_top;

      VW_DEBUG( vw_out(DebugMessage) << "DiskImageResourceTIFF reading strip " << strip 
                << " (rows " << strip_top << "-" << strip_top+strip_rows-1 << ") from " << m_filename << std::endl; )
      TIFFReadEncodedStrip( tif, strip, buf, (tsize_t) -1 );

      strip_src.data = ((uint8*)buf) + bbox.min().x()*strip_src.cstride + (strip_top-strip*rows_per_strip)*strip_src.rstride;
      strip_src.format.rows = strip_rows;

      strip_dest.data = ((uint8*)dest.data) + (strip_top-bbox.min().y())*strip_dest.rstride;
      strip_dest.format.rows = strip_rows;

      convert( strip_dest, strip_src );
    }
  
    _TIFFfree(buf);
  }
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
