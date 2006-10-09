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

/// \file DiskImageResourceJPEG.cc
/// 
/// Provides support for the JPEG file format.
///

#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <vector>
#include <vw/FileIO/DiskImageResourceJPEG.h>
#include <vw/Core/Exception.h>

extern "C" {
#include <jpeglib.h>
}


// This class and the following method are part of the structure for
// taking control of error handling from the JPEG library.  Here we
// are overriding the default error handling method with one of our
// own that throws an exception.
METHODDEF(void) vw_jpeg_error_handler (j_common_ptr cinfo) {

  // cinfo->err really points to a my_error_mgr struct, so coerce
  // pointer and return to the last location where setjmp() was
  // calledq.
  char error_msg[2048];
  (*cinfo->err->format_message) (cinfo, error_msg);
  throw vw::IOErr() << "DiskImageResourceJPEG: A libjpeg error occurred. " << error_msg;
}


/// Close the JPEG file when the object is destroyed
vw::DiskImageResourceJPEG::~DiskImageResourceJPEG() {
  this->flush();
}
  
/// Flush the buffered data to disk
void vw::DiskImageResourceJPEG::flush() {

  if (m_file_ptr) {
    fclose((FILE*)m_file_ptr);  
    m_file_ptr = NULL;
  }

}

/// Bind the resource to a file for reading.  Confirm that we can open
/// the file and that it has a sane pixel format.  
void vw::DiskImageResourceJPEG::open( std::string const& filename )
{
  if(m_file_ptr)
    throw IOErr() << "DiskImageResourceJPEG: A file is already open.";

  // Set up the JPEG data structures
  jpeg_decompress_struct cinfo;

  // Set up error handling.  If the JPEG library encounters an error,
  // it will return here and throw an exception.
  jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jerr.error_exit = vw_jpeg_error_handler;
  jpeg_create_decompress(&cinfo);

  // Open the file on disk
	FILE * infile;
	if ((infile = fopen(filename.c_str(), "rb")) == NULL) {
    throw vw::IOErr() << "Failed to open \"" << filename << "\" using libJPEG.";
	}
	jpeg_stdio_src(&cinfo, infile);
  m_filename = filename;
  m_file_ptr = infile;

  // Read the header information and parse it
	jpeg_read_header(&cinfo, TRUE);

  // If the user has requested a subsampled version of the original
  // image, we pass this request along to the JPEG driver here.
  cinfo.scale_num = 1; 
  cinfo.scale_denom = m_subsample_factor;

  // Start the decompression routine
  jpeg_start_decompress(&cinfo);
  
  m_format.cols = cinfo.output_width;
  m_format.rows = cinfo.output_height;
  m_format.channel_type = VW_CHANNEL_UINT8;
    
  switch (cinfo.output_components) {
  case 1:  m_format.pixel_format = VW_PIXEL_GRAY;   m_format.planes=1; break;
  case 2:  m_format.pixel_format = VW_PIXEL_GRAYA;  m_format.planes=1; break;
  case 3:  m_format.pixel_format = VW_PIXEL_RGB;    m_format.planes=1; break;
  case 4:  m_format.pixel_format = VW_PIXEL_RGBA;   m_format.planes=1; break;
  default: m_format.planes = cinfo.output_components; m_format.pixel_format = VW_PIXEL_SCALAR; break;
  }

  jpeg_abort_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);
  
  // Rewind the file so that we can restart the jpeg header reading
  // process when the user calls read().
  fseek(infile, 0, SEEK_SET);
}

/// Bind the resource to a file for writing.
void vw::DiskImageResourceJPEG::create( std::string const& filename, 
                                        GenericImageFormat const& format )
{
  if( format.planes!=1 && format.pixel_format!=VW_PIXEL_SCALAR )
    throw NoImplErr() << "JPEG doesn't support multi-plane images with compound pixel types.";
  if(m_file_ptr)
    throw IOErr() << "DiskImageResourceJPEG: A file is already open.";

  // Open the file on disk
	FILE * outfile;
	if ((outfile = fopen(filename.c_str(), "wb")) == NULL) {
    throw vw::IOErr() << "Failed to open \"" << filename << "\" using libJPEG.";
	}
  m_filename = filename;
  m_format = format;
  m_file_ptr = outfile;

  // The JPEG file format only supports 8-bit channel types, so we
  // force that setting here.  
  m_format.channel_type = VW_CHANNEL_UINT8;
  
  // The JPEG file format only supports color spaces without an alpha
  // channel, so we truncate the alpha channel here and print a
  // warning message.
  if (format.pixel_format == VW_PIXEL_GRAYA) {
    m_format.pixel_format = VW_PIXEL_GRAY;
    vw_out(WarningMessage) << "  Warning: alpha channel removed.  ";
  } else if (format.pixel_format == VW_PIXEL_RGBA) {
    m_format.pixel_format = VW_PIXEL_RGB;
    vw_out(WarningMessage) << "  Warning: alpha channel removed.  ";
  }
}

/// Read the disk image into the given buffer.
void vw::DiskImageResourceJPEG::read( GenericImageBuffer const& dest ) const
{
  VW_ASSERT( dest.format.cols==cols() && dest.format.rows==rows(),
             IOErr() << "Buffer has wrong dimensions in JPEG read." );

  // Set up the JPEG data structures
  jpeg_decompress_struct cinfo;

  // Set up error handling.  If the JPEG library encounters an error,
  // it will return here and throw an exception.
  jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jerr.error_exit = vw_jpeg_error_handler;

  jpeg_create_decompress(&cinfo);
	jpeg_stdio_src(&cinfo, (FILE*)m_file_ptr);

  // Read the header information and parse it
	jpeg_read_header(&cinfo, TRUE);

  // If the user has requested a grayscale image, we can speed up the
  // decoding process by forcing the jpeg driver to decode the image
  // in grayscale (since the color components need not be processed).
  if (dest.format.pixel_format == VW_PIXEL_GRAY) 
    cinfo.out_color_space = JCS_GRAYSCALE;

  // If the user has requested a subsampled version of the original
  // image, we pass this request along to the JPEG driver here.
  cinfo.scale_num = 1; 
  cinfo.scale_denom = m_subsample_factor;

  // Start the decompression
  jpeg_start_decompress(&cinfo);

  int scanline_size = cinfo.output_width * cinfo.output_components;
  JSAMPARRAY scanline = (*cinfo.mem->alloc_sarray)
		((j_common_ptr) &cinfo, JPOOL_IMAGE, scanline_size, 1);
  uint8* buf = new uint8[scanline_size * cinfo.output_height];
  while (cinfo.output_scanline < cinfo.output_height) {
    jpeg_read_scanlines(&cinfo, scanline, 1);
    
    // Copy the data over into a contiguous buffer, which is what
    // convert() expects when we call it below.  This extra copy of
    // the data is undesirable, but probably not a huge performance
    // hit compared to reading the data from disk.
    for (int i = 0; i < scanline_size; i++) 
      buf[scanline_size * (cinfo.output_scanline-1) + i] = (uint8)(scanline[0][i]);
  }

  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);
  
  // Set up a generic image buffer around the tdata_t buf object that
  // jpeg used to copy it's data.
  GenericImageBuffer src;
  src.data = (uint8*)buf;
  src.format = m_format;

  // This is part of the grayscale optimization above.
  if (dest.format.pixel_format == VW_PIXEL_GRAY) 
    src.format.pixel_format = VW_PIXEL_GRAY;

  src.cstride = scanline_size / cinfo.output_width;
  src.rstride = scanline_size;
  src.pstride = scanline_size * cinfo.output_height;
  convert( dest, src );
}

// Write the given buffer into the disk image.
void vw::DiskImageResourceJPEG::write( GenericImageBuffer const& src )
{
  VW_ASSERT( src.format.cols==cols() && src.format.rows==rows(),
             IOErr() << "Buffer has wrong dimensions in JPEG write." );

  // Set up the JPEG data structures
  jpeg_compress_struct cinfo;

  // Set up error handling.  If the JPEG library encounters an error,
  // it will return here and throw an exception.
  jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jerr.error_exit = vw_jpeg_error_handler;

  jpeg_create_compress(&cinfo);
	jpeg_stdio_dest(&cinfo, (FILE*)m_file_ptr);

  cinfo.image_width = m_format.cols;
  cinfo.image_height = m_format.rows;

  if (m_format.pixel_format == VW_PIXEL_SCALAR) {
    cinfo.input_components = m_format.planes;
    cinfo.in_color_space = JCS_UNKNOWN;
  } else if (m_format.pixel_format == VW_PIXEL_GRAY) {
    cinfo.input_components = 1;
    cinfo.in_color_space = JCS_GRAYSCALE;
  } else if (m_format.pixel_format == VW_PIXEL_RGB) {
    cinfo.input_components = 3;
    cinfo.in_color_space = JCS_RGB;
  } else {
    throw IOErr() << "DiskImageResourceJPEG: Unsupported pixel type.";
  }

  // Set up the default values for the header and set the compression
  // quality
  jpeg_set_defaults(&cinfo);
  jpeg_set_quality(&cinfo, (int)(100*m_quality), TRUE); // limit to baseline-JPEG values 

  // Set up the generic image buffer and convert the data into this buffer
  GenericImageBuffer dst;
  dst.data = new uint8[cinfo.image_width*cinfo.input_components*cinfo.image_height];
  dst.format = m_format;
  dst.cstride = num_channels(m_format.pixel_format) * channel_size(m_format.channel_type);
  dst.rstride = dst.cstride * m_format.cols;
  dst.pstride = dst.rstride * m_format.rows;
  convert( dst, src );

  jpeg_start_compress(&cinfo, TRUE);
  int row_stride = cinfo.image_width*cinfo.input_components;	

  // Write the image data to disk.
  JSAMPROW row_pointer[1];	
  while (cinfo.next_scanline < cinfo.image_height) {
    row_pointer[0] = &(((uint8*)dst.data)[cinfo.next_scanline * row_stride]);
    jpeg_write_scanlines(&cinfo, row_pointer, 1);
  }  

  // Clean up
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);
}

// A FileIO hook to open a file for reading
vw::DiskImageResource* vw::DiskImageResourceJPEG::construct_open( std::string const& filename ) {
  return new DiskImageResourceJPEG( filename );
}

// A FileIO hook to open a file for writing
vw::DiskImageResource* vw::DiskImageResourceJPEG::construct_create( std::string const& filename,
                                                                    GenericImageFormat const& format ) {
  return new DiskImageResourceJPEG( filename, format );
}
