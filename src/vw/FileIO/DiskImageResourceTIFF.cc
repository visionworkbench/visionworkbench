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
    std::string filename;

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
  m_use_compression = false;
  open( filename );
}
    
vw::DiskImageResourceTIFF::DiskImageResourceTIFF( std::string const& filename, 
                                                  vw::ImageFormat const& format )
  : DiskImageResource( filename ), m_info( new DiskImageResourceInfoTIFF() )
{
  m_use_compression = false;
  create( filename, format );
}

vw::Vector2i vw::DiskImageResourceTIFF::native_block_size() const {
  return m_info->block_size;
}

/// Bind the resource to a file for reading.  Confirm that we can open
/// the file and that it has a sane pixel format.  
void vw::DiskImageResourceTIFF::open( std::string const& filename ) {
  m_info->filename = filename;

  TIFFSetWarningHandler( &tiff_warning_handler );
  TIFFSetErrorHandler( &tiff_error_handler );

  TIFF* tif = TIFFOpen( filename.c_str(), "r" );
  if( !tif ) vw_throw( vw::IOErr() << "DiskImageResourceTIFF: Failed to open \"" << filename << "\" for reading!" );

  // Read into temp variables first to ensure we are using the right integer type.
  // Otherwise we can run into endianness problems.
  uint32 cols_tmp, rows_tmp;
  uint16 planes_tmp;
  TIFFGetField( tif, TIFFTAG_IMAGEWIDTH, &cols_tmp );
  TIFFGetField( tif, TIFFTAG_IMAGELENGTH, &rows_tmp );
  TIFFGetFieldDefaulted( tif, TIFFTAG_SAMPLESPERPIXEL, &planes_tmp );
  m_format.cols = cols_tmp;
  m_format.rows = rows_tmp;
  m_format.planes = planes_tmp;

  uint16 sample_format = 0, bits_per_sample = 0, photometric = 0;
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

  uint16 plane_configuration = 0;
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

  TIFFClose(tif);
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

  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, (uint32)m_format.cols);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH, (uint32)m_format.rows);
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, (uint16)(8*channel_size(m_format.channel_type)));

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

  if (m_use_compression)
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
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, (uint16)m_format.planes);
  } else {
    // Compound pixel types are stored contiguously in TIFF files
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, (uint16)num_channels(m_format.pixel_format));
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

  TIFF* tif = TIFFOpen( m_info->filename.c_str(), "r" );
  if( !tif ) vw_throw( vw::IOErr() << "DiskImageResourceTIFF: Failed to open \"" << m_filename << "\" for reading!" );

  uint16 config = 0, bpsample = 0, nsamples = 0, photometric = 0;
  TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &config);
  TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bpsample);
  TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &nsamples);
  TIFFGetField(tif, TIFFTAG_PHOTOMETRIC, &photometric);

  bool is_planar = (config == PLANARCONFIG_SEPARATE) && (m_format.pixel_format != VW_PIXEL_SCALAR);
  bool is_tiled = TIFFIsTiled(tif);

  // Compute the tile or strip geometry
  uint32 block_cols, block_rows, block_size, blocks_per_row, blocks_per_plane;
  if( is_tiled ) {
    TIFFGetField(tif, TIFFTAG_TILEWIDTH, &block_cols);
    TIFFGetField(tif, TIFFTAG_TILELENGTH, &block_rows);
    block_size = TIFFTileSize(tif);
    blocks_per_row = (cols()-1) / block_cols + 1;
    blocks_per_plane = blocks_per_row * ( (rows()-1) / block_rows + 1 );
  }
  else {
    block_cols = cols();
    TIFFGetField( tif, TIFFTAG_ROWSPERSTRIP, &block_rows );
    block_size = TIFFStripSize(tif);
    blocks_per_row = 1;
    blocks_per_plane = (rows()-1) / block_rows + 1;
  }

  tdata_t buf = _TIFFmalloc( block_size );

  // Allocate a buffer interleave planar data
  tdata_t plane_buf = 0;
  if( is_planar ) {
    plane_buf = buf;
    buf = _TIFFmalloc( block_size * nsamples );
  }

  // Palettized TIFFs are always uint16 RGB
  tdata_t palette_buf = 0;
  uint16 *red_table, *green_table, *blue_table;
  if( photometric == PHOTOMETRIC_PALETTE ) {
    palette_buf = buf;
    buf = _TIFFmalloc( block_cols*block_rows*6 );
    TIFFGetField( tif, TIFFTAG_COLORMAP, &red_table, &green_table, &blue_table );
  }

  // Set up the source and destination image buffers
  ImageBuffer src_buf, dest_buf=dest;
  src_buf.format = m_format;
  if( photometric == PHOTOMETRIC_PALETTE ) src_buf.cstride = 6;
  else src_buf.cstride = bpsample * nsamples / 8;
  src_buf.rstride = block_cols*src_buf.cstride;
  src_buf.pstride = block_rows*src_buf.rstride;
    
  for( int block_y = bbox.min().y()/block_rows; block_y <= (bbox.max().y()-1)/block_rows; ++block_y ) {
    for( int block_x = bbox.min().x()/block_cols; block_x <= (bbox.max().x()-1)/block_cols; ++block_x ) {
      int block_id = block_y * blocks_per_row + block_x;
      int data_left = (std::max)(block_x*block_cols,uint32(bbox.min().x()))-block_x*block_cols;
      int data_top  = (std::max)(block_y*block_rows,uint32(bbox.min().y()))-block_y*block_rows;
      int data_right = (std::min)((block_x+1)*block_cols,uint32(bbox.max().x()))-block_x*block_cols;
      int data_bottom = (std::min)((block_y+1)*block_rows,uint32(bbox.max().y()))-block_y*block_rows;

      // Read the block into the buffer, converting planar or palettized data as needed.
      if( is_planar ) {
        // At the moment we make an extra copy here to spoof plane contiguity
        for( int i=0; i<nsamples; ++i ) {
          if( is_tiled ) TIFFReadEncodedTile( tif, block_id+nsamples*blocks_per_plane, plane_buf, (tsize_t) -1 );
          else TIFFReadEncodedStrip( tif, block_id+nsamples*blocks_per_plane, plane_buf, (tsize_t) -1 );
          // Oh man, this is horrible!
          switch(bpsample/8) {
          case 1:
            for( int y=data_top; y<data_bottom; ++y ) {
              for( int x=data_left; x<data_right; ++x ) {
                ((uint8*)buf)[(y*block_cols+x)*nsamples+i] = ((uint8*)plane_buf)[y*block_cols+x];
              }
            }
            break;
          case 2:
            for( int y=data_top; y<data_bottom; ++y ) {
              for( int x=data_left; x<data_right; ++x ) {
                ((uint16*)buf)[(y*block_cols+x)*nsamples+i] = ((uint16*)plane_buf)[y*block_cols+x];
              }
            }
            break;
          case 4:
            for( int y=data_top; y<data_bottom; ++y ) {
              for( int x=data_left; x<data_right; ++x ) {
                ((uint32*)buf)[(y*block_cols+x)*nsamples+i] = ((uint32*)plane_buf)[y*block_cols+x];
              }
            }
            break;
          default:
            vw_throw( NoImplErr() << "Unsupported bit depth in separate-plane TIFF!" );
          }
        }
      }
      else if( photometric == PHOTOMETRIC_PALETTE ) {
        if( is_tiled ) TIFFReadEncodedTile( tif, block_id, palette_buf, (tsize_t) -1 );
        else TIFFReadEncodedStrip( tif, block_id, palette_buf, (tsize_t) -1 );
        if( photometric == PHOTOMETRIC_PALETTE ) {
          for( int y=data_top; y<data_bottom; ++y ) {
            for( int x=data_left; x<data_right; ++x ) {
              int i = y*block_cols + x;
              int p = ((uint8*)palette_buf)[i];
              ((uint16*)buf)[3*i+0] = red_table[p];
              ((uint16*)buf)[3*i+1] = green_table[p];
              ((uint16*)buf)[3*i+2] = blue_table[p];
            }
          }
        }
      }
      else {
        if( is_tiled ) TIFFReadEncodedTile( tif, block_id, buf, (tsize_t) -1 );
        else TIFFReadEncodedStrip( tif, block_id, buf, (tsize_t) -1 );
      }
      
      src_buf.data = (uint8*)buf + data_left*src_buf.cstride + data_top*src_buf.rstride;
      dest_buf.data = (uint8*)dest.data + (data_left+block_x*block_cols-bbox.min().x())*dest.cstride + (data_top+block_y*block_rows-bbox.min().y())*dest.rstride;
      src_buf.format.cols = dest_buf.format.cols = data_right-data_left;
      src_buf.format.rows = dest_buf.format.rows = data_bottom-data_top;
      
      convert( dest_buf, src_buf );
    }
  
  }

  _TIFFfree(buf);
  if( plane_buf ) _TIFFfree(plane_buf);
  if( palette_buf ) _TIFFfree(palette_buf);
  TIFFClose(tif);
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
