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

#ifdef WIN32
#define snprintf _snprintf
#endif

namespace vw {
  class DiskImageResourceInfoTIFF {
  public:
    TIFF *tif;
    Vector2i block_size;
    std::string filename;
    int current_line;
    bool striped;

    DiskImageResourceInfoTIFF() : tif(0), block_size(), current_line(0) {}
    ~DiskImageResourceInfoTIFF() {
      close();
    }

    void reopen_read() {
      close();
      tif = TIFFOpen(filename.c_str(), "r");
      if( !tif ) vw_throw( vw::IOErr() << "DiskImageResourceTIFF: Failed to open \"" << filename << "\" for reading!" );
      current_line = 0;
    }

    void close() {
      if( tif ) {
        TIFFClose(tif);
        tif=NULL;
      }
    }
  };
}

/// Handle libTIFF warning conditions by outputting message text at the
/// DebugMessage verbosity level.
static void tiff_warning_handler(const char* module, const char* frmt, va_list ap) {
  char msg[VW_ERROR_BUFFER_SIZE];
  vsnprintf( msg, VW_ERROR_BUFFER_SIZE, frmt, ap );
  VW_OUT(vw::WarningMessage, "fileio") << "DiskImageResourceTIFF (" << (module?module:"none") << ") Warning: " << msg << std::endl;
}

/// Handle libTIFF error conditions by vw_throwing an IOErr with the
/// message text.
/*static void tiff_error_handler(const char* module, const char* frmt, va_list ap) {
  char msg[VW_ERROR_BUFFER_SIZE];
  vsnprintf( msg, VW_ERROR_BUFFER_SIZE, frmt, ap );
  vw_throw( vw::IOErr() << "DiskImageResourceTIFF (" << (module?module:"none") << ") Error: " << msg );
}
*/

// Handle libTIFF error conditions by writing the error and hope the calling
// program checks the return value for the function
static char tiff_error_msg[VW_ERROR_BUFFER_SIZE];
static void tiff_error_handler(const char* module, const char* frmt, va_list ap) {
  char msg[VW_ERROR_BUFFER_SIZE];
  vsnprintf( msg, VW_ERROR_BUFFER_SIZE, frmt, ap );
  snprintf( tiff_error_msg, VW_ERROR_BUFFER_SIZE,
    "DiskImageResourceTIFF (%s) Error: %s",
    (module?module:"none"), msg );
}


vw::DiskImageResourceTIFF::DiskImageResourceTIFF( std::string const& filename )
  : DiskImageResource( filename ), m_info( new DiskImageResourceInfoTIFF() )
{
  m_use_compression = false;
  open( filename );
}

vw::DiskImageResourceTIFF::DiskImageResourceTIFF( std::string const& filename,
                                                  vw::ImageFormat const& format,
                                                  bool use_compression )
  : DiskImageResource( filename ), m_info( new DiskImageResourceInfoTIFF() ),
    m_use_compression( use_compression )
{
  create( filename, format );
}

vw::Vector2i vw::DiskImageResourceTIFF::block_read_size() const {
  return m_info->block_size;
}

/// Bind the resource to a file for reading.  Confirm that we can open
/// the file and that it has a sane pixel format.
void vw::DiskImageResourceTIFF::open( std::string const& filename ) {
  m_info->filename = filename;

  TIFFSetWarningHandler( &tiff_warning_handler );
  TIFFSetErrorHandler( &tiff_error_handler );

  TIFF* tif = TIFFOpen( filename.c_str(), "r" );
  if( !tif ) vw_throw( vw::ArgumentErr() << "DiskImageResourceTIFF: Failed to open \"" << filename << "\" for reading!" );

  // Read into temp variables first to ensure we are using the right integer type.
  // Otherwise we can run into endianness problems.
  uint32 cols_tmp, rows_tmp;
  uint16 planes_tmp;
  check_retval(TIFFGetField( tif, TIFFTAG_IMAGEWIDTH, &cols_tmp ), 0);
  check_retval(TIFFGetField( tif, TIFFTAG_IMAGELENGTH, &rows_tmp ), 0);
  check_retval(
    TIFFGetFieldDefaulted( tif, TIFFTAG_SAMPLESPERPIXEL, &planes_tmp ), 0);
  m_format.cols = cols_tmp;
  m_format.rows = rows_tmp;
  m_format.planes = planes_tmp;

  uint16 sample_format = 0, bits_per_sample = 0, photometric = 0;
  check_retval(
    TIFFGetFieldDefaulted( tif, TIFFTAG_BITSPERSAMPLE, &bits_per_sample ), 0);
  check_retval(
    TIFFGetFieldDefaulted( tif, TIFFTAG_SAMPLEFORMAT, &sample_format ), 0);
  check_retval(TIFFGetField( tif, TIFFTAG_PHOTOMETRIC, &photometric ), 0);

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
  check_retval(TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &plane_configuration), 0);

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
    check_retval(TIFFGetField( tif, TIFFTAG_TILEWIDTH, &tile_width ), 0);
    check_retval(TIFFGetField( tif, TIFFTAG_TILELENGTH, &tile_length ), 0);
    m_info->block_size = Vector2i(tile_width,tile_length);
  }
  else {
    uint32 rows_per_strip;
    check_retval(TIFFGetField( tif, TIFFTAG_ROWSPERSTRIP, &rows_per_strip ), 0);
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
  if( !tif  ) vw_throw( vw::ArgumentErr() << "Failed to create \"" << m_filename << "\" using libTIFF." );

  check_retval(TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, (uint32)m_format.cols), 0);
  check_retval(TIFFSetField(tif, TIFFTAG_IMAGELENGTH, (uint32)m_format.rows), 0);
  check_retval(TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, (uint16)(8*channel_size(m_format.channel_type))), 0);

  if (m_format.pixel_format == VW_PIXEL_RGB ||
      m_format.pixel_format == VW_PIXEL_RGBA) {
    check_retval(TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB), 0);
  } else {
    check_retval(TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK), 0);
  }

  check_retval(TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, 1), 0);
  check_retval(TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, RESUNIT_INCH), 0);
  check_retval(TIFFSetField(tif, TIFFTAG_XRESOLUTION, 70.0), 0);
  check_retval(TIFFSetField(tif, TIFFTAG_YRESOLUTION, 70.0), 0);

  if (m_use_compression) {
    check_retval(TIFFSetField(tif, TIFFTAG_COMPRESSION,COMPRESSION_LZW), 0);
  }

  switch (m_format.channel_type) {
  case VW_CHANNEL_INT8:
  case VW_CHANNEL_INT16:
  case VW_CHANNEL_INT32:
  case VW_CHANNEL_INT64:
    check_retval(TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_INT), 0);
    break;
  case VW_CHANNEL_UINT8:
  case VW_CHANNEL_UINT16:
  case VW_CHANNEL_UINT32:
  case VW_CHANNEL_UINT64:
    check_retval(TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT), 0);
    break;
  case VW_CHANNEL_FLOAT16:
  case VW_CHANNEL_FLOAT32:
  case VW_CHANNEL_FLOAT64:
    check_retval(TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP), 0);
    break;
  default:
    vw_throw( IOErr() << "DiskImageResourceTIFF: Unsupported VW channel type." );
  }

  if (m_format.pixel_format == VW_PIXEL_SCALAR) {
    // Multi-plane images with simple pixel types are stored in separate
    // planes in the TIFF image.
    check_retval(TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_SEPARATE), 0);
    check_retval(
      TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, (uint16)m_format.planes), 0);
  } else {
    // Compound pixel types are stored contiguously in TIFF files
    check_retval(TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG), 0);
    check_retval(TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, (uint16)num_channels(m_format.pixel_format)), 0);
  }

  uint32 rows_per_strip = TIFFDefaultStripSize( tif, 0 );
  check_retval(TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, rows_per_strip), 0);
  m_info->block_size = Vector2i(cols(),rows_per_strip);

  m_info->tif = tif;
}

/// Read the disk image into the given buffer.
void vw::DiskImageResourceTIFF::read( ImageBuffer const& dest, BBox2i const& bbox ) const
{
  VW_ASSERT( int(dest.format.cols)==bbox.width() && int(dest.format.rows)==bbox.height(),
             ArgumentErr() << "DiskImageResourceTIFF (read) Error: Destination buffer has wrong dimensions!" );

  // Only support sequential reading on striped TIFFs right now.
  if( !m_info || !(m_info->tif) || !(m_info->striped) || (m_info->striped && m_info->current_line > bbox.min().y()) )
    m_info->reopen_read();

  uint16 config = 0, bpsample = 0, nsamples = 0, photometric = 0;
  check_retval(TIFFGetField(m_info->tif, TIFFTAG_PLANARCONFIG, &config), 0);
  check_retval(TIFFGetField(m_info->tif, TIFFTAG_BITSPERSAMPLE, &bpsample), 0);
  check_retval(TIFFGetField(m_info->tif, TIFFTAG_SAMPLESPERPIXEL, &nsamples), 0);
  check_retval(TIFFGetField(m_info->tif, TIFFTAG_PHOTOMETRIC, &photometric), 0);

  bool is_planar = (config == PLANARCONFIG_SEPARATE) && (m_format.pixel_format != VW_PIXEL_SCALAR);
  bool is_tiled = TIFFIsTiled(m_info->tif);
  m_info->striped = !is_tiled;

  // Compute the tile or strip geometry
  uint32 block_cols, block_rows, block_size, blocks_per_row, blocks_per_plane;
  if( is_tiled ) {
    check_retval(TIFFGetField(m_info->tif, TIFFTAG_TILEWIDTH, &block_cols), 0);
    check_retval(TIFFGetField(m_info->tif, TIFFTAG_TILELENGTH, &block_rows), 0);
    block_size = TIFFTileSize(m_info->tif);
    blocks_per_row = (cols()-1) / block_cols + 1;
    blocks_per_plane = blocks_per_row * ( (rows()-1) / block_rows + 1 );
  }
  else {
    block_cols = cols();
    check_retval(TIFFGetField( m_info->tif, TIFFTAG_ROWSPERSTRIP, &block_rows ), 0);
    block_size = TIFFStripSize(m_info->tif);
    blocks_per_row = 1;
    blocks_per_plane = (rows()-1) / block_rows + 1;
  }

  tdata_t buf = _TIFFmalloc( block_size );
  if( !buf ) vw_throw( vw::IOErr() << "DiskImageResourceTIFF: Failed to malloc!" );

  // Allocate a buffer interleave planar data
  tdata_t plane_buf = 0;
  if( is_planar ) {
    plane_buf = buf;
    buf = _TIFFmalloc( block_size * nsamples );
    if( !buf ) vw_throw( vw::IOErr() << "DiskImageResourceTIFF: Failed to malloc!" );
  }

  // Palettized TIFFs are always uint16 RGB
  tdata_t palette_buf = 0;
  uint16 *red_table, *green_table, *blue_table;
  if( photometric == PHOTOMETRIC_PALETTE ) {
    palette_buf = buf;
    buf = _TIFFmalloc( block_cols*block_rows*6 );
    if( !buf ) vw_throw( vw::IOErr() << "DiskImageResourceTIFF: Failed to malloc!" );

    check_retval(TIFFGetField( m_info->tif, TIFFTAG_COLORMAP, &red_table, &green_table, &blue_table ), 0);
  }

  // Set up the source and destination image buffers
  ImageBuffer src_buf, dest_buf=dest;
  src_buf.format = m_format;
  if( photometric == PHOTOMETRIC_PALETTE ) src_buf.cstride = 6;
  else src_buf.cstride = bpsample * nsamples / 8;
  src_buf.rstride = block_cols*src_buf.cstride;
  src_buf.pstride = block_rows*src_buf.rstride;

  for( int block_y = bbox.min().y()/block_rows; block_y <= int((bbox.max().y()-1)/block_rows); ++block_y ) {
    for( int block_x = bbox.min().x()/block_cols; block_x <= int((bbox.max().x()-1)/block_cols); ++block_x ) {
      int block_id = block_y * blocks_per_row + block_x;
      int data_left = (std::max)(block_x*block_cols,uint32(bbox.min().x()))-block_x*block_cols;
      int data_top  = (std::max)(block_y*block_rows,uint32(bbox.min().y()))-block_y*block_rows;
      int data_right = (std::min)((block_x+1)*block_cols,uint32(bbox.max().x()))-block_x*block_cols;
      int data_bottom = (std::min)((block_y+1)*block_rows,uint32(bbox.max().y()))-block_y*block_rows;

      // Read the block into the buffer, converting planar or palettized data as needed.
      if( is_planar ) {
        // At the moment we make an extra copy here to spoof plane contiguity
        for( int i=0; i<nsamples; ++i ) {
          if( is_tiled ) check_retval(TIFFReadEncodedTile( m_info->tif, block_id+i*blocks_per_plane, plane_buf, (tsize_t) -1 ), -1);
          else check_retval(TIFFReadEncodedStrip( m_info->tif, block_id+i*blocks_per_plane, plane_buf, (tsize_t) -1 ), 0);
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
        if( is_tiled ) check_retval(TIFFReadEncodedTile( m_info->tif, block_id, palette_buf, (tsize_t) -1 ), -1);
        else check_retval(TIFFReadEncodedStrip( m_info->tif, block_id, palette_buf, (tsize_t) -1 ), 0);
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
        if( is_tiled )  {
          check_retval(TIFFReadEncodedTile( m_info->tif, block_id, buf, (tsize_t) -1 ), -1);
        } else {
          check_retval(TIFFReadEncodedStrip( m_info->tif, block_id, buf, (tsize_t) -1 ), -1);
          m_info->current_line++;
        }
      }

      src_buf.data = (uint8*)buf + data_left*src_buf.cstride + data_top*src_buf.rstride;
      dest_buf.data = (uint8*)dest.data + (data_left+block_x*block_cols-bbox.min().x())*dest.cstride + (data_top+block_y*block_rows-bbox.min().y())*dest.rstride;
      src_buf.format.cols = dest_buf.format.cols = data_right-data_left;
      src_buf.format.rows = dest_buf.format.rows = data_bottom-data_top;

      convert( dest_buf, src_buf, m_rescale );
    }

  }

  _TIFFfree(buf);
  if( plane_buf ) _TIFFfree(plane_buf);
  if( palette_buf ) _TIFFfree(palette_buf);
  // Sorry .. this keeps us from incrementally reading
  m_info->close();
}

// Write the given buffer into the disk image.
void vw::DiskImageResourceTIFF::write( ImageBuffer const& src, BBox2i const& bbox )
{
  VW_ASSERT(bbox.width() == m_format.cols,
            ArgumentErr() << "DiskImageResourceTIFF: bounding box must be the same width as image.\n");

  // Allocate some buffer memory for the output data
  uint32 scanline_size = num_channels(m_format.pixel_format) * channel_size(m_format.channel_type) * m_format.cols;
  tdata_t buf = _TIFFmalloc(scanline_size);
  if( !buf ) vw_throw( vw::IOErr() << "DiskImageResourceTIFF: Failed to malloc!" );

  // Set up the image buffer and convert the data into this buffer.
  ImageBuffer dst(m_format, buf);

  ImageBuffer src_plane = src;
  src_plane.format.rows = 1;
  src_plane.format.planes = 1;
  dst.format.rows = 1;
  dst.format.planes = 1;

  // Write the image data to disk.
  for (int32 p = 0; p < m_format.planes; p++) {
    ImageBuffer src_row = src_plane;
    for (int32 row = 0; row < bbox.height(); row++) {
      convert( dst, src_row, m_rescale );
      check_retval(TIFFWriteScanline(m_info->tif, (uint8*)buf, bbox.min()[1]+row, p), -1);
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

// Helper routine for TIFF functions to check their return values and throw
// if there was an error.
void vw::DiskImageResourceTIFF::check_retval(const int retval, const int error_val) const {
  if (retval == error_val) {
    vw_throw( vw::IOErr() << "check_retval: " << tiff_error_msg );
  }
}

