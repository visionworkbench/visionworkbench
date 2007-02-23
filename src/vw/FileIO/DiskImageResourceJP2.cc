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

/// \file DiskImageResourceJP2.cc
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
#include <strings.h>			   // for strncasecmp

#include <boost/algorithm/string.hpp>

#include <openjpeg.h>

#include <vw/Core/Exception.h>
#include <vw/FileIO/DiskImageResourceJP2.h>

using namespace std;
using namespace boost;

extern "C" {
int cio_numbytesleft(opj_cio_t *cio);
}

namespace vw
{

  /// Close the file when the object is destroyed
  DiskImageResourceJP2::~DiskImageResourceJP2()
  {
    this->flush();
  }
 
  /// Flush the buffered data to disk
  void
  DiskImageResourceJP2::flush()
  {
  }

#if 0
  static int
  get_file_format(char *filename)
  {
    static const char *extension[] = { "j2k", "jp2", "jpt" };
    static const int format[] = { J2K_CFMT, JP2_CFMT, JPT_CFMT };
    char * ext = strrchr(filename, '.') + 1;

    if (ext != 0)
    {
      for (int i = 0; i < sizeof(format); i++)
      {
	if (strncasecmp(ext, extension[i], 3) == 0)
	  return format[i];
      }
    }

    return -1;
  }
#endif
  // sample error callback expecting a FILE* client object
  static void
  error_callback(const char *msg, void *client_data)
  {
    FILE *stream = (FILE*)client_data;
    fprintf(stream, "[ERROR] %s", msg);
  }
  // sample warning callback expecting a FILE* client object
  void
  warning_callback(const char *msg, void *client_data)
  {
    FILE *stream = (FILE*)client_data;
    fprintf(stream, "[WARNING] %s", msg);
  }
  // sample debug callback expecting a FILE* client object
  void
  info_callback(const char *msg, void *client_data)
  {
    FILE *stream = (FILE*)client_data;
    fprintf(stream, "[INFO] %s", msg);
  }

  static void
  setup_event_handling(opj_common_ptr cinfo)
  {
    opj_event_mgr_t event_manager;

    // configure the event callbacks (not required)
    // setting of each callback is optionnal 
    memset(&event_manager, 0, sizeof(opj_event_mgr_t));
    event_manager.error_handler = error_callback;
    event_manager.warning_handler = warning_callback;
    event_manager.info_handler = info_callback;
    // catch events using our callbacks and give a local context
    opj_set_event_mgr((opj_common_ptr)cinfo, &event_manager, stderr);
  }

  static opj_image_t *
  init_jpeg2k_image(opj_cparameters_t &compress_parameters,
		    int width, int height, int num_components,
		    OPJ_COLOR_SPACE color_space,
		    ImageBuffer &dst)
  {
    cout << "-- init_jpeg2k_image(): start" << endl;
    opj_image_t *new_image = new opj_image_t;

    if (new_image != 0)
    {
      new_image->color_space = color_space;
      new_image->numcomps = num_components;

      // allocate memory for the per-component information
      new_image->comps = new opj_image_comp_t [new_image->numcomps];
	
      if (new_image->comps == 0)
      {
	delete new_image;
	return 0;
      }

      int dx = compress_parameters.subsampling_dx;
      int dy = compress_parameters.subsampling_dy;

      // init the individual image components
      for (int i = 0; i < num_components; i++)
      {
	opj_image_comp_t *component = &new_image->comps[i];

	component->dx = dx;
	component->dy = dy;
	component->w = width;
	component->h = height;
	component->x0 = 0;
	component->y0 = 0;
	component->prec = 8;
	component->bpp = 8;
	component->sgnd = 0;
	component->resno_decoded = 0;
	component->factor = 0;

	component->data = (int*) calloc(width * height, sizeof(int));
	if (component->data == 0)
	{
	  cout << "init_jpeg2k_image(): could not allocate memory for data"
	       << endl;
	  delete [] new_image->comps;
	  delete new_image;
	  return 0;
	}
	// maybe we'll eventaully be able to just point it at the data...
// 	component->data = (int *) dst(0, 0, i);
	cout << "init_jpeg2k_image(): initializing data" << endl;
	// Only have had luck with 8 bit images so far...
	uint8 *values = (uint8 *) dst(0, 0, i);
	int minval = values[0];
	int maxval = minval;
	cout << "values[0] = " << int(values[0]) << endl;
	for (int index = 0; index < width * height; index++)
	{
	  component->data[index] = (int) values[index];
	  if (minval > component->data[index])
	    minval = component->data[index];
	  if (maxval < component->data[index])
	    maxval = component->data[index];;
	}
	cout << "done" << endl;
	cout << "minval = " << minval << endl;
	cout << "maxval = " << maxval << endl;
      }

      new_image->x0 = compress_parameters.image_offset_x0;
      new_image->y0 = compress_parameters.image_offset_y0;
      new_image->x1 = new_image->x0 + ((width - 1) * dx) + 1;
      new_image->y1 = new_image->y0 + ((height - 1) * dy) + 1;
    }

    cout << "-- init_jpeg2k_image(): end" << endl;

    return new_image;
  }

  /// Bind the resource to a file for reading.  Confirm that we can open
  /// the file and that it has a sane pixel format.  
  void
  DiskImageResourceJP2::open(std::string const& filename)
  {
    vw_throw(NoImplErr() << "DiskImageResourceJP2 (read) Error: "
	     "JPEG2000 file formats are not supported yet!");
#if 0
    ifstream image_file(filename.c_str(), ios::in | ios::binary);

    if (!image_file)
      throw IOErr() << "DiskImageResourceJP2::open(): could not open \""
		    << filename << "\".";

    // - - - - - - - - - - - -
    // decode the code-stream
    //
    // Note: A JPEG 2000 file (extension .jp2 and magic number
    // 0000000C6A5020200D0A870A) consists of a number of "boxes" that
    // contain metadata or code streams (i.e., encoded image data).
    //
    // A file with a .j2k extension contains only a "raw" code stream,
    // which actually has a header:
    //
    // FF4F -- Start of Codestream (SOC) marker
    // FF51 -- Size (SIZ) marker
    // #### -- Length of SIZ segment (Lsiz)
    // ...  -- SIZ parameters
    // FF5C -- Quantization Default (QCD) marker
    // #### -- Length of QCD segment (Lqcd)
    // ...  -- QCD parameters
    // FF52 -- Coding style Default (COD) marker
    // #### -- Length of COD segment (Lcod)
    // ...  -- COD parameters
    // Optional marker segments COC (FF53), QCC (FF5D), RGN (FF5E),
    // POD (FF5F), PPM (FF60), PLM (FF57), TLM (FF55), CME (Comment
    // and Extension info, FF64)
    //
    // followed by encoded image data:
    //
    // FF90 -- Start of Tile (SOT) marker
    // ...
    // FF93 -- Start of data (SOD) marker
    // ...
    // FFD9 -- End of Codestream (EOC) marker
    //
    // - - - - - - - - - - - -

    opj_dinfo_t* dinfo = 0;		   // handle to a decompressor
    opj_cio_t *cio = 0;
    opj_image_t *image = 0;

    int file_format = get_file_format(file_format.c_str();

    // We can handle *.j2k, *.jp2, *.jpc or *.jpt files
    switch (file_format)
    {
    case J2K_CFMT:			   // JPEG-2000 codestream
      {
	// get a decoder handle
	dinfo = opj_create_decompress(CODEC_J2K);
	
	// catch events using our callbacks and give a local context
	opj_set_event_mgr((opj_common_ptr)dinfo, &event_mgr, stderr);

	// setup the decoder decoding parameters using user parameters
	opj_setup_decoder(dinfo, &parameters);

	// open a byte stream
	cio = opj_cio_open((opj_common_ptr)dinfo, src, file_length);

	// decode the stream and fill the image structure
	image = opj_decode(dinfo, cio);
	if (!image)
	{
	  cerr << "ERROR in DiskImageResourceJP2::open(): "
	       << "failed to decode image!" << endl;
	  opj_destroy_decompress(dinfo);
	  opj_cio_close(cio);
	  return 1;
	}
	
	// close the byte stream
	opj_cio_close(cio);
      }
      break;
    case JP2_CFMT:			   // JPEG 2000 compressed image data
      {
	// get a decoder handle
	dinfo = opj_create_decompress(CODEC_JP2);
	
	// catch events using our callbacks and give a local context
	opj_set_event_mgr((opj_common_ptr)dinfo, &event_mgr, stderr);

	// setup the decoder decoding parameters using the current
	// image and using user parameters
	opj_setup_decoder(dinfo, &parameters);

	// open a byte stream
	cio = opj_cio_open((opj_common_ptr)dinfo, src, file_length);

	// decode the stream and fill the image structure
	image = opj_decode(dinfo, cio);
	if(!image)
	{
	  cerr << "ERROR in DiskImageResourceJP2::open(): "
	       << "failed to decode image!" << endl;
	  opj_destroy_decompress(dinfo);
	  opj_cio_close(cio);
	  return 1;
	}

	// close the byte stream
	opj_cio_close(cio);
      }
      break;
    case JPT_CFMT:			   // JPEG 2000, JPIP
      {
	// get a decoder handle
	dinfo = opj_create_decompress(CODEC_JPT);
	
	// catch events using our callbacks and give a local context
	opj_set_event_mgr((opj_common_ptr)dinfo, &event_mgr, stderr);

	// setup the decoder decoding parameters using user parameters
	opj_setup_decoder(dinfo, &parameters);

	// open a byte stream
	cio = opj_cio_open((opj_common_ptr)dinfo, src, file_length);

	// decode the stream and fill the image structure
	image = opj_decode(dinfo, cio);
	if (!image)
	{
	  cerr << "ERROR in DiskImageResourceJP2::open(): "
	       << "failed to decode image!" << endl;
	  opj_destroy_decompress(dinfo);
	  opj_cio_close(cio);
	  return 1;
	}

	// close the byte stream
	opj_cio_close(cio);
      }
      break;
    default:
      cerr << "ERROR in DiskImageResourceJP2::open(): "
	   << "unknown JPEG2000 image format variant" << endl;
      return 1;
      break;
    }

    m_format.planes = 1;
    m_format.pixel_format = VW_PIXEL_GRAY;
    m_format.cols = header.bytesPerScanline / bytes_per_pixel;
    m_format.rows = header.numScanlines;
 
    image_file.close();
#endif
  }

  /// Bind the resource to a file for writing.
  void
  DiskImageResourceJP2::create(std::string const& filename,
			       ImageFormat const& format)
  {
    cout << "DiskImageResourceJP2::create(): enter" << endl;

    VW_ASSERT(format.planes == 1 || format.pixel_format == VW_PIXEL_SCALAR,
	      NoImplErr() << "DiskImageResourceOpenJP2: Cannot create "
	      << filename << "\n\tThe image cannot have both multiple "
	      << "channels and multiple planes.\n");

    m_filename = filename;
    m_format = format;

    std::cout << "DiskImageResourceJP2::create(): channel type = " << std::flush;

    switch (format.channel_type)
    {
    case VW_CHANNEL_UNKNOWN:
      cout << "VW_CHANNEL_UNKNOWN" << endl;
      break;
    case VW_CHANNEL_INT8:
      cout << "VW_CHANNEL_INT8" << endl;
      break;
    case VW_CHANNEL_UINT8:
      cout << "VW_CHANNEL_UINT8" << endl;
      break;
    case VW_CHANNEL_INT16:
      cout << "VW_CHANNEL_INT16" << endl;
      break;
    case VW_CHANNEL_UINT16:
      cout << "VW_CHANNEL_UINT16" << endl;
      break;
    case VW_CHANNEL_INT32:
      cout << "VW_CHANNEL_INT32" << endl;
      break;
    case VW_CHANNEL_UINT32:
      cout << "VW_CHANNEL_UINT32" << endl;
      break;
    case VW_CHANNEL_INT64:
      cout << "VW_CHANNEL_INT64" << endl;
      break;
    case VW_CHANNEL_UINT64:
      cout << "VW_CHANNEL_UINT64" << endl;
      break;
    case VW_CHANNEL_FLOAT16:
      cout << "VW_CHANNEL_FLOAT16" << endl;
      break;
    case VW_CHANNEL_FLOAT32:
      cout << "VW_CHANNEL_FLOAT32" << endl;
      break;
    case VW_CHANNEL_FLOAT64:
      cout << "VW_CHANNEL_FLOAT64" << endl;
      break;
    case VW_CHANNEL_BOOL:
      cout << "VW_CHANNEL_BOOL" << endl;
      break;
    case VW_CHANNEL_USER:
      cout << "VW_CHANNEL_USER" << endl;
      break;
    }

    // The open JPEG code only seems to work with 8 bit images right now...
    m_format.channel_type = VW_CHANNEL_UINT8;
    m_format.planes = std::max(format.planes,
			       num_channels(format.pixel_format));

    cout << "DiskImageResourceJP2::create(): return" << endl;
  }

  /// Read the disk image into the given buffer.
  void
    DiskImageResourceJP2::read(ImageBuffer const& dest,
			       BBox2i const& bbox) const
  {
    VW_ASSERT((bbox.width() == int(cols())) && (bbox.height()==int(rows())),
	      NoImplErr() << "DiskImageResourceOpenJP2 does not support"
	      " partial reads." );
#if 0
    VW_ASSERT((dest.format.cols == cols()) && (dest.format.rows == rows()),
	      IOErr() << "Buffer has wrong dimensions in JP2 read.");
  
    ifstream image_file(m_filename.c_str(), ios::in | ios::binary);

    if (!image_file)
      vw_throw(IOErr() << "DiskImageResourceJP2::read(): "
	       "failed to open \"" << m_filename << "\" for reading!");

    // Set the file offset to the position of the first image byte
    // (read in from the previously opened JP2 header)
    image_file.seekg(m_image_data_offset, ios::beg);

    image_file.read((char *) image_data, bytes_per_pixel * total_pixels);

    if (image_file.bad())
      vw_throw(IOErr() << "DiskImageResourceJP2::read(): "
	       "an unrecoverable error occured while reading the image data.");

    // Create image buffer from the JP2 data.
    ImageBuffer src;
    src.data = image_data;
    src.format = m_format;
    src.cstride = bytes_per_pixel;
    src.rstride = bytes_per_pixel * m_format.cols;
    src.pstride = bytes_per_pixel * m_format.cols * m_format.rows;
  
    convert(dest, src);

    delete[] image_data;

    image_file.close();
#endif
  }

  // Write the given buffer into the disk image.
  void
    DiskImageResourceJP2::write(ImageBuffer const& src, BBox2i const& bbox)
  {
    VW_ASSERT((bbox.width()==int(cols())) && (bbox.height()==int(rows())),
	      NoImplErr() << "DiskImageResourceOpenJP2 does not support"
	      " partial writes." );

    VW_ASSERT((src.format.cols == cols()) && (src.format.rows == rows()),
	      IOErr() << "Buffer has wrong dimensions in "
	      "DiskImageResourceJP2::write().");

    VW_ASSERT(src.format.planes == 1, IOErr() << "Buffer has multiple planes "
	      "in DiskImageResourceJP2::write(). Only gray scale "
	      "images currently handled.");
  
    // This is pretty simple since we always write 32 bit integer
    // files.  Note that we handle multi-channel images with
    // interleaved planes.
    // We've already ensured that either planes == 1 or channels == 1.
    // Can't seem to get anything other than 8 bit JP2s working right now...
    ImageView<uint8> jp2_image_view(m_format.cols, m_format.rows,
				    m_format.planes);
    ImageBuffer dst = jp2_image_view.buffer();

    convert(dst, src);

    cout << "DiskImageResourceJP2::write(): "
	 << "initializing encoder and parameters" << endl;


    // For the moment we only handle grayscale images
    OPJ_COLOR_SPACE color_space = CLRSPC_GRAY;
    int num_components = 1;
    opj_cparameters_t coding_parameters;

    opj_set_default_encoder_parameters(&coding_parameters);
    coding_parameters.cp_comment = "Created by OpenJPEG version 0.9";
    // if no rate entered, lossless by default
    if (coding_parameters.tcp_numlayers == 0)
    {
      cout << "DiskImageResourceJP2::write(): encoding lossless!" << endl;
      coding_parameters.tcp_rates[0] = 0;	// MOD antonin : losslessbug
      coding_parameters.tcp_numlayers++;
      coding_parameters.cp_disto_alloc = 1;
    }
    // experiment with tiles -- LJE
    coding_parameters.cp_tdx = 128;
    coding_parameters.cp_tdy = 128;
    coding_parameters.cp_tx0 = 0;
    coding_parameters.cp_ty0 = 0;
    coding_parameters.tile_size_on = true;

    opj_image_t *jpeg2k_image = init_jpeg2k_image(coding_parameters,
						  cols(), rows(),
						  num_components,
						  color_space, dst);

    int *values = jpeg2k_image->comps[0].data;
    int minval = values[0];
    int maxval = minval;
    cout << "DiskImageResourceJP2::write(): values[0] = " << values[0] << endl;
    cout << "DiskImageResourceJP2::write(): width = " << jpeg2k_image->comps[0].w << endl;
    cout << "DiskImageResourceJP2::write(): height = " << jpeg2k_image->comps[0].h
	 << endl;
    for (int index = 0;
	 index < jpeg2k_image->comps[0].w * jpeg2k_image->comps[0].h;
	 index++)
    {
      if (minval > values[index])
	minval = values[index];
      if (maxval < values[index])
	maxval = values[index];;
    }
    cout << "DiskImageResourceJP2::write(): minval = " << minval << endl;
    cout << "DiskImageResourceJP2::write(): maxval = " << maxval << endl;

    cout << "DiskImageResourceJP2::write(): "
	 << "creating a JP2 encoder" << endl;

    // get a JP2 compression context pointer
    opj_cinfo_t* cinfo = opj_create_compress(CODEC_JP2);

    // Set up event manager (not necessary)
    setup_event_handling(opj_common_ptr(cinfo));

    cout << "DiskImageResourceJP2::write(): "
	 << "setting up JP2 encoder" << endl;


    cout << "num resolutions = " << coding_parameters.numresolution << endl;
    // setup the encoder parameters using the current image
    opj_setup_encoder(cinfo, &coding_parameters, jpeg2k_image);

    cout << "DiskImageResourceJP2::write(): "
	 << "opening codec i/o stream" << endl;

    // open a codec input/output byte stream for writing
    opj_cio_t *cio = opj_cio_open(opj_common_ptr(cinfo), 0, 0);
    if (cio == 0)
    {
      cout << "DiskImageResourceJP2::write(): "
	   << "failed to open a codec input output stream." << endl;

      vw_throw(IOErr() << "DiskImageResourceJP2::write(): "
	       "failed to open a codec input output stream.");
    }

    cout << "DiskImageResourceJP2::write(): "
	 << "encoding image" << endl;

    cout << "jpeg2k_image x0 = " << jpeg2k_image->x0 << endl;
    cout << "jpeg2k_image y0 = " << jpeg2k_image->y0 << endl;
    cout << "jpeg2k_image width = " << jpeg2k_image->x1 << endl;
    cout << "jpeg2k_image height = " << jpeg2k_image->y1 << endl;
    cout << "jpeg2k_image->comps = " << jpeg2k_image->comps << endl;
    cout << "jpeg2k_image->numcomps = " << jpeg2k_image->numcomps << endl;

    bool success = opj_encode(cinfo, cio, jpeg2k_image,
			      coding_parameters.index);
    if (!success)
    {
      cout << "DiskImageResourceJP2::write(): "
	   << "failed to encode image." << endl;

      opj_cio_close(cio);
      vw_throw(IOErr() << "DiskImageResourceJP2::write(): "
	       "failed to encode image.");
    }

    int codestream_length = cio_tell(cio);
    int bytes_to_end =  cio_numbytesleft(cio);

    cout << "DiskImageResourceJP2::write(): "
	 << "writing encoded image: " << m_filename << endl;

    cout << "DiskImageResourceJP2::write(): "
	 << "encoded buffer length = " << codestream_length << endl;
    cout << "DiskImageResourceJP2::write(): "
	 << "bytes to end of encode stream = " << bytes_to_end << endl;

    // write the buffer to disk
    FILE *output_fp = fopen(m_filename.c_str(), "wb");
    if (output_fp == 0)
      vw_throw(IOErr() << "DiskImageResourceJP2::write(): "
	       "failed to open " << m_filename);

    fwrite(cio->buffer, 1, codestream_length, output_fp);

    cout << "DiskImageResourceJP2::write(): "
	 << "closing streams and files" << endl;

    fclose(output_fp);

    // close and free the byte stream
    opj_cio_close(cio);

    // free remaining compression structures
    opj_destroy_compress(cinfo);

    cout << "DiskImageResourceJP2::write(): "
	 << "cleaning up" << endl;

    // delete image structures
//     delete [] jpeg2k_image->comps;
//     delete jpeg2k_image;
  }

  // A FileIO hook to open a file for reading
  DiskImageResource*
  DiskImageResourceJP2::construct_open(std::string const& filename)
  {
    return new DiskImageResourceJP2(filename);
  }

  // A FileIO hook to open a file for writing
  DiskImageResource*
  DiskImageResourceJP2::construct_create(std::string const& filename,
					 ImageFormat const& format)
  {
    return new DiskImageResourceJP2(filename, format);
  }
}
