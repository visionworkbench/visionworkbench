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

namespace vw {
  class DiskImageResourceInfoTIFF {
  public:
    TIFF *tif;
    Vector2i block_size;

    DiskImageResourceInfoTIFF() : tif(0), block_size() {}
    ~DiskImageResourceInfoTIFF() { 
      if( tif ) TIFFClose(tif);
    }
  };
}

/// Handle libTIFF warning conditions by outputting message text at the 
/// DebugMessage verbosity level.
static void tiff_warning_handler(const char* module, const char* frmt, va_list ap) {
  char msg[VW_ERROR_BUFFER_SIZE];
  vsnprintf( msg, VW_ERROR_BUFFER_SIZE, frmt, ap );
  vw::vw_out(vw::VerboseDebugMessage+1) << "DiskImageResourceTIFF (" << (module?module:"none") << ") Warning: " << msg << std::endl;
}

/// Handle libTIFF error conditions by vw_throwing an IOErr with the 
/// message text.
static void tiff_error_handler(const char* module, const char* frmt, va_list ap) {
  char msg[VW_ERROR_BUFFER_SIZE];
  vsnprintf( msg, VW_ERROR_BUFFER_SIZE, frmt, ap );
  vw_throw( vw::IOErr() << "DiskImageResourceTIFF (" << (module?module:"none") << ") Error: " << msg );
}

vw::DiskImageResourceTIFF::DiskImageResourceTIFF( std::string const& filename )
  : DiskImageResource( filename ), m_info( new DiskImageResourceInfoTIFF() )
{
  open( filename );
}
    
vw::DiskImageResourceTIFF::DiskImageResourceTIFF( std::string const& filename, 
                                                  vw::ImageFormat const& format )
  : DiskImageResource( filename ), m_info( new DiskImageResourceInfoTIFF() )
{
  create( filename, format );
}

vw::Vector2i vw::DiskImageResourceTIFF::native_block_size() const {
  return m_info->block_size;
}

/// Bind the resource to a file for reading.  Confirm that we can open
/// the file and that it has a sane pixel format.  
void vw::DiskImageResourceTIFF::open( std::string const& filename ) {
  TIFFSetWarningHandler( &tiff_warning_handler );
  TIFFSetErrorHandler( &tiff_error_handler );

  TIFF* tif = TIFFOpen( filename.c_str(), "r" );
  if( !tif ) vw_throw( vw::IOErr() << "DiskImageResourceTIFF: Failed to open \"" << filename << "\" for reading!" );

  TIFFGetField( tif, TIFFTAG_IMAGEWIDTH, &(m_format.cols) );
  TIFFGetField( tif, TIFFTAG_IMAGELENGTH, &(m_format.rows) );
  TIFFGetFieldDefaulted( tif, TIFFTAG_SAMPLESPERPIXEL, &(m_format.planes) );

  uint32 sample_format = 0, bits_per_sample = 0, photometric = 0;
  TIFFGetFieldDefaulted( tif, TIFFTAG_BITSPERSAMPLE, &bits_per_sample );
  TIFFGetFieldDefaulted( tif, TIFFTAG_SAMPLEFORMAT, &sample_format );
  TIFFGetField( tif, TIFFTAG_PHOTOMETRIC, &photometric );

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
    vw_throw( IOErr() << "DiskImageResourceTIFF: " << m_filename << " has an unsupported channel type ("
              << sample_format << "," << bits_per_sample << ")!" );
  }

  uint32 plane_configuration = 0;
  TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &plane_configuration);

  // FIXME: Tiff might actually provide us with some info on
  // colorimetric interpretation of the channels, so maybe we should
  // try to use that here as well?
  if( photometric == PHOTOMETRIC_PALETTE ) {
    m_format.channel_type = VW_CHANNEL_UINT16;
    m_format.pixel_format = VW_PIXEL_RGB;
  }
  else if( m_format.planes == 1 ) {
    m_format.pixel_format = VW_PIXEL_GRAY;
  }
  else {
    switch( m_format.planes ) {
    case 2:  m_format.pixel_format = VW_PIXEL_GRAYA;  m_format.planes=1; break;
    case 3:  m_format.pixel_format = VW_PIXEL_RGB;    m_format.planes=1; break;
    case 4:  m_format.pixel_format = VW_PIXEL_RGBA;   m_format.planes=1; break;
    default: m_format.pixel_format = VW_PIXEL_SCALAR; break;
    }
  }
  
  if( TIFFIsTiled(tif) ) {
    uint32 tile_width, tile_length;
    TIFFGetField( tif, TIFFTAG_TILEWIDTH, &tile_width );
    TIFFGetField( tif, TIFFTAG_TILELENGTH, &tile_length );
    m_info->block_size = Vector2i(tile_width,tile_length);
  }
  else {
    uint32 rows_per_strip;
    TIFFGetField( tif, TIFFTAG_ROWSPERSTRIP, &rows_per_strip );
    m_info->block_size = Vector2i(cols(),rows_per_strip);
  }

  m_info->tif = tif;
}

/// Bind the resource to a file for writing.
void vw::DiskImageResourceTIFF::create( std::string const& filename, 
                                        ImageFormat const& format )
{
  if( format.planes!=1 && format.pixel_format!=VW_PIXEL_SCALAR )
    vw_throw( NoImplErr() << "TIFF doesn't support multi-plane images with compound pixel types." );

  // Set the TIFF warning and error handlers to Vision Workbench
  // functions, so that we can handle them ourselves.
  TIFFSetWarningHandler(&tiff_warning_handler);
  TIFFSetErrorHandler(&tiff_error_handler);

  m_format = format;

  TIFF* tif = TIFFOpen(m_filename.c_str(), "w");
  if( !tif  ) vw_throw( vw::IOErr() << "Failed to create \"" << m_filename << "\" using libTIFF." );

  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, m_format.cols);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH, m_format.rows);
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8*channel_size(m_format.channel_type));

  if (m_format.pixel_format == VW_PIXEL_RGB ||
      m_format.pixel_format == VW_PIXEL_RGBA) {
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  } else {
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
  }

  TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, 1);
  TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, RESUNIT_INCH);
  TIFFSetField(tif, TIFFTAG_XRESOLUTION, 70.0);
  TIFFSetField(tif, TIFFTAG_YRESOLUTION, 70.0);

  TIFFSetField(tif, TIFFTAG_COMPRESSION,COMPRESSION_LZW);

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
    vw_throw( IOErr() << "DiskImageResourceTIFF: Unsupported VW channel type." );
  }

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

  uint32 rows_per_strip = TIFFDefaultStripSize( tif, 0 );
  TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, rows_per_strip);
  m_info->block_size = Vector2i(cols(),rows_per_strip);
  
  m_info->tif = tif;
}

/// Read the disk image into the given buffer.
void vw::DiskImageResourceTIFF::read( ImageBuffer const& dest, BBox2i const& bbox ) const
{
  VW_ASSERT( int(dest.format.cols)==bbox.width() && int(dest.format.rows)==bbox.height(),
             ArgumentErr() << "DiskImageResourceTIFF (read) Error: Destination buffer has wrong dimensions!" );

  TIFF* tif = m_info->tif;
  if( !tif ) vw_throw( vw::IOErr() << "DiskImageResourceTIFF: Failed to open \"" << m_filename << "\" for reading!" );

  if( TIFFIsTiled(tif) ) {
    vw_throw( NoImplErr() << "DiskImageResourceTIFF (read) Error: Reading from tile-based TIFF files is not yet supported!" );
  }
  else {
    uint32 config = 0, nsamples = 0, photometric = 0;
    TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &config);
    TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &nsamples);
    TIFFGetField( tif, TIFFTAG_PHOTOMETRIC, &photometric );

    tdata_t buf=0, plane_buf=0;
    if( config==PLANARCONFIG_SEPARATE && m_format.pixel_format != VW_PIXEL_SCALAR ) {
      plane_buf = _TIFFmalloc( TIFFStripSize(tif) );
      buf = _TIFFmalloc( TIFFStripSize(tif) * nsamples );
    }
    else {
      buf = _TIFFmalloc( TIFFStripSize(tif) );
    }

    uint32 rows_per_strip;
    uint32 scanline_size = TIFFScanlineSize(tif);    
    TIFFGetField( tif, TIFFTAG_ROWSPERSTRIP, &rows_per_strip );

    ImageBuffer strip_src, strip_dest=dest;
    strip_src.format = m_format;
    strip_src.format.cols = bbox.width();

    uint16 *red_table, *green_table, *blue_table;
    uint16 *palette_buf = 0;
    if( photometric == PHOTOMETRIC_PALETTE ) {
      TIFFGetField( tif, TIFFTAG_COLORMAP, &red_table, &green_table, &blue_table );
      palette_buf = new uint16[rows_per_strip*m_format.cols*3];
      strip_src.cstride = 6;
      strip_src.rstride = m_format.cols*6;
      strip_src.pstride = rows_per_strip*m_format.cols*6;
    }
    else if( config==PLANARCONFIG_SEPARATE && m_format.pixel_format != VW_PIXEL_SCALAR ) {
      strip_src.cstride = nsamples * scanline_size / m_format.cols;
      strip_src.rstride = nsamples * scanline_size;
      strip_src.pstride = nsamples * scanline_size * m_format.rows;
    }
    else {
      strip_src.cstride = scanline_size / m_format.cols;
      strip_src.rstride = scanline_size;
      strip_src.pstride = scanline_size * m_format.rows;
    }

    int minstrip=bbox.min().y()/rows_per_strip, maxstrip=(bbox.max().y()-1)/rows_per_strip;
    int nstrips = (rows()-1)/rows_per_strip + 1;
    for( int strip=minstrip; strip<=maxstrip; ++strip ) {
      int strip_top = std::max(strip*rows_per_strip,uint32(bbox.min().y()));
      int strip_rows = std::min((strip+1)*rows_per_strip,uint32(bbox.max().y()))-strip_top;

      VW_DEBUG( vw_out(DebugMessage) << "DiskImageResourceTIFF reading strip " << strip 
                << " (rows " << strip_top << "-" << strip_top+strip_rows-1 << ") from " << m_filename << std::endl; );

      if( config==PLANARCONFIG_SEPARATE && m_format.pixel_format != VW_PIXEL_SCALAR ) {
        // At the moment we make an extra copy here to spoof plane contiguity
        for( int i=0; i<nsamples; ++i ) {
          TIFFReadEncodedStrip( tif, strip+i*nstrips, plane_buf, (tsize_t) -1 );
          for( int y=0; y<rows_per_strip; ++y ) {
            for( int x=bbox.min().x(); x<bbox.max().x(); ++x ) {
              ((uint8*)buf)[(y*m_format.cols+x)*nsamples+i] = ((uint8*)plane_buf)[y*m_format.cols+x];
            }
          }
        }
        strip_src.data = ((uint8*)buf) + bbox.min().x()*strip_src.cstride + (strip_top-strip*rows_per_strip)*strip_src.rstride;
      }
      else {

        TIFFReadEncodedStrip( tif, strip, buf, (tsize_t) -1 );

        // FIXME: This could stand some serious tidying up / optimizing
        if( photometric == PHOTOMETRIC_PALETTE ) {
          for( int y=0; y<strip_rows; ++y ) {
            for( int x=bbox.min().x(); x<bbox.max().x(); ++x ) {
              int p = ((uint8*)buf)[ (strip_top-strip*rows_per_strip)*scanline_size + x*(scanline_size/m_format.cols) ];
              palette_buf[3*(m_format.cols*y+x)+0] = red_table[p];
              palette_buf[3*(m_format.cols*y+x)+1] = green_table[p];
            palette_buf[3*(m_format.cols*y+x)+2] = blue_table[p];
            }
          }
          strip_src.data = ((uint8*)palette_buf) + bbox.min().x()*strip_src.cstride + (strip_top-strip*rows_per_strip)*strip_src.rstride;
        }
        else {
          strip_src.data = ((uint8*)buf) + bbox.min().x()*strip_src.cstride + (strip_top-strip*rows_per_strip)*strip_src.rstride;
        }

      }
      strip_src.format.rows = strip_rows;
        
      strip_dest.data = ((uint8*)dest.data) + (strip_top-bbox.min().y())*strip_dest.rstride;
      strip_dest.format.rows = strip_rows;
        
      convert( strip_dest, strip_src );
    }
  
    delete[] palette_buf;
    _TIFFfree(buf);
  }
}

// Write the given buffer into the disk image.
void vw::DiskImageResourceTIFF::write( ImageBuffer const& src, BBox2i const& bbox )
{
  VW_ASSERT(bbox.width() == m_format.cols, 
            ArgumentErr() << "DiskImageResourceTIFF: bounding box must be the same width as image.\n");

  // Allocate some buffer memory for the output data
  uint32 scanline_size = num_channels(m_format.pixel_format) * channel_size(m_format.channel_type) * m_format.cols;
  tdata_t buf = _TIFFmalloc(scanline_size);

  // Set up the image buffer and convert the data into this buffer.
  ImageBuffer dst;
  dst.data = (uint8*)buf;
  dst.format = m_format;
  dst.cstride = num_channels(m_format.pixel_format) * channel_size(m_format.channel_type);
  dst.rstride = dst.cstride * m_format.cols;
  dst.pstride = dst.rstride * m_format.rows;

  ImageBuffer src_plane = src;
  src_plane.format.rows = 1;
  src_plane.format.planes = 1;
  dst.format.rows = 1;
  dst.format.planes = 1;

  // Write the image data to disk.
  for (uint32 p = 0; p < m_format.planes; p++) {
    ImageBuffer src_row = src_plane;
    for (uint32 row = 0; row < bbox.height(); row++) {
      convert( dst, src_row );
      TIFFWriteScanline(m_info->tif, (uint8*)buf, bbox.min()[1]+row, p);
      src_row.data = (uint8*)src_row.data + src_row.rstride;
    }
    src_plane.data = (uint8*)src_plane.data + src_plane.pstride;
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
                                                                    ImageFormat const& format ) {
  return new DiskImageResourceTIFF( filename, format );
}
