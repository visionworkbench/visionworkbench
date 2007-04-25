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

/// \file FileIO/DiskImageResourceOpenEXR.cc
/// 
/// Provides support for the OpenEXR file format.
///

#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <vector>

#include <ImfInputFile.h>
#include <ImfOutputFile.h>
#include <ImfChannelList.h>
#include <ImfStringAttribute.h>
#include <ImfMatrixAttribute.h>
#include <ImfArray.h>
#include <ImfChannelList.h>

#include <vw/Core/Exception.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/Statistics.h>
#include <vw/FileIO/DiskImageResourceOpenEXR.h>

namespace {

  // This little type computation routine helps us to determine what
  // to label the channels in the OpenEXR file given the pixel type of
  // the source image.
  static std::string openexr_channel_string_of_pixel_type(const int pixel_format, 
                                                          const int channel) {
    if (pixel_format == vw::VW_PIXEL_RGB) {
      switch (channel) {
      case 0 : return "R"; break;
      case 1 : return "G"; break;
      case 2 : return "B"; break;
      default: vw_throw( vw::ArgumentErr() << "ChannelStringOfPixelType: Invalid channel number (" << channel << ")" );
      }

    } else if (pixel_format == vw::VW_PIXEL_RGBA) {
      switch (channel) {
      case 0 : return "R"; break;
      case 1 : return "G"; break;
      case 2 : return "B"; break;
      case 3 : return "A"; break;
      default: vw_throw( vw::ArgumentErr() << "ChannelStringOfPixelType: Invalid channel number (" << channel << ")" );
      }
    
    }
    // Default case:
    std::ostringstream m_stream;
    m_stream << "Channel" << channel; 
    return m_stream.str();

  }

}


// The destructor is here, despite being so brief, because deleting 
// an object safely requires knowing its full type.
vw::DiskImageResourceOpenEXR::~DiskImageResourceOpenEXR() {
  if (m_input_file_ptr) 
    delete m_input_file_ptr;
  if (m_output_file_ptr) 
    delete m_output_file_ptr;
}

vw::Vector2i vw::DiskImageResourceOpenEXR::native_block_size() const {
  return m_block_size;
}

// Bind the resource to a file for reading.  Confirm that we can open
// the file and that it has a sane pixel format.  In general VIL does 
// not give us any useful information about the pixel format, so we 
// make some guesses here based on the channel count.
void vw::DiskImageResourceOpenEXR::open( std::string const& filename )
{
  try {
    // Open Image file and read the header
    m_filename = filename;
    
    // Check to make sure that the file_ptr is not already in use.
    if (m_input_file_ptr) 
      vw_throw( IOErr() << "Disk image resources do not yet support reuse." );
    
    m_input_file_ptr = new Imf::InputFile(filename.c_str());
    Imf::FrameBuffer frameBuffer;
    
    // Find the width and height of the image 
    Imath::Box2i dw = m_input_file_ptr->header().dataWindow();
    m_format.cols  = int(dw.max.x - dw.min.x + 1);
    m_format.rows  = int(dw.max.y - dw.min.y + 1);
    
    // Determine the number of image channels 
    Imf::ChannelList::ConstIterator iter = m_input_file_ptr->header().channels().begin();
    int num_channels = 0;
    while( iter++ != m_input_file_ptr->header().channels().end() )
      num_channels++;
    m_format.planes = num_channels;
    
    // For now, we only support reading in multi-plane, single channel
    // images.
    m_format.pixel_format = VW_PIXEL_SCALAR;
    m_block_size = Vector2i(m_format.cols,m_openexr_rows_per_block);
    
  } catch (Iex::BaseExc e) {
    vw_throw( vw::IOErr() << "DiskImageResourceOpenEXR: could not open " << filename << ":\n\t" << e.what() ); 
  } 
}

// Bind the resource to a file for writing.  
void vw::DiskImageResourceOpenEXR::create( std::string const& filename, 
                                           ImageFormat const& format )
{
  VW_ASSERT(format.planes == 1 || format.pixel_format==VW_PIXEL_SCALAR,
            NoImplErr() << "DiskImageResourceOpenEXR: Cannot create " << filename << "\n\t"
            << "The image cannot have both multiple channels and multiple planes.\n");
  
  m_filename = filename;
  m_format = format;
  m_format.channel_type = VW_CHANNEL_FLOAT32;
  m_format.planes = std::max( format.planes, num_channels( format.pixel_format ) );
  m_block_size = Vector2i(m_format.cols,m_openexr_rows_per_block);

  // Open the EXR file and set up the header information
  m_labels.resize(m_format.planes);
  
  try {      
    // Create the file header with the appropriate number of
    // channels.  Label the channels in order starting with "Channel 0".
    Imf::Header header (m_format.cols,m_format.rows);
    for ( int32 nn = 0; nn < m_format.planes; nn++) {
      m_labels[nn] = openexr_channel_string_of_pixel_type(m_format.pixel_format, nn);
      //      std::cout << "Writing channel " << nn << ": " << labels[nn] << "\n";
      header.channels().insert (m_labels[nn].c_str(), Imf::Channel (Imf::FLOAT));
    }
    
    // Open the file handle and create an empty framebuffer object. 
    m_output_file_ptr = new Imf::OutputFile(m_filename.c_str(), header);

  } catch (Iex::BaseExc e) {
    vw_throw( vw::IOErr() << "DiskImageResourceOpenEXR: Failed to create " << m_filename << ".\n\t" << e.what() );
  }
}

// Read the disk image into the given buffer.
void vw::DiskImageResourceOpenEXR::read( ImageBuffer const& dest, BBox2i const& bbox ) const
{
  
  if (!m_input_file_ptr) 
    vw_throw( LogicErr() << "DiskImageResourceOpenEXR: Could not read file. No file has been opened." );
  
  try {
    // Find the width and height of the image and set the data window
    // to the beginning of the requesed block.
    Imath::Box2i dw = m_input_file_ptr->header().dataWindow();
    dw.min.y += bbox.min().y();
    int32 height = bbox.height();
    int32 width = bbox.width();

    // Set up the OpenEXR data structures necessary to read all of
    // the channels out of the file and execute the call to readPixels().
    Imf::Array2D<float> *inputArrays[m_format.planes];
    Imf::ChannelList::ConstIterator channel = m_input_file_ptr->header().channels().begin();
    std::vector<std::string> channel_names(m_format.planes);
    for (int i=0; channel != m_input_file_ptr->header().channels().end(); ++channel, ++i) {
      channel_names[i] = channel.name();
    }
    
    // OpenEXR seems to order channels in the file alphabetically
    // (dumb!), rather than in the order in which they were saved.
    // This means that we need to reorder the channel names when
    // they are labelled as RGB or RGBA.  For other channel naming
    // schemes, we just go with alphabetical, since that's all we've
    // got.
    if ( m_format.planes == 3 ) {
      if (find(channel_names.begin(), channel_names.end(), "R") != channel_names.end() &&
          find(channel_names.begin(), channel_names.end(), "G") != channel_names.end() &&
          find(channel_names.begin(), channel_names.end(), "B") != channel_names.end()) {
        channel_names[0] = "R";
        channel_names[1] = "G";
        channel_names[2] = "B";
      }
    } else if ( m_format.planes == 4 ) {
      if (find(channel_names.begin(), channel_names.end(), "R") != channel_names.end() &&
          find(channel_names.begin(), channel_names.end(), "G") != channel_names.end() &&
          find(channel_names.begin(), channel_names.end(), "B") != channel_names.end() &&
          find(channel_names.begin(), channel_names.end(), "A") != channel_names.end()) {
        channel_names[0] = "R";
        channel_names[1] = "G";
        channel_names[2] = "B";
        channel_names[3] = "A";
      }
    }
    
    Imf::FrameBuffer frameBuffer;
    for ( int32 nn = 0; nn < m_format.planes; ++nn ) {
      inputArrays[nn] = new Imf::Array2D<float>(height,width);
      //        std::cout << "Reading channel " << channel_names[nn] << "\n";
      frameBuffer.insert (channel_names[nn].c_str(), Imf::Slice (Imf::FLOAT, (char *) (&(*inputArrays[nn])[-dw.min.y][-dw.min.x]),
                                                                 sizeof ((*inputArrays[nn])[0][0]) * 1, 
                                                                 sizeof ((*inputArrays[nn])[0][0]) * width, 1, 1, 0.0)); 
    }
    m_input_file_ptr->setFrameBuffer (frameBuffer);
    m_input_file_ptr->readPixels (dw.min.y, std::min(int(dw.min.y + (height-1)), dw.max.y));
    
    // Copy the pixels over into a ImageView object.
    // 
    // Recast to the templatized pixel type in the process.
    ImageView<float> src_image(width, height, m_format.planes);
    for ( int32 nn=0; nn<m_format.planes; ++nn ) {
      for ( int32 i=0; i<width; ++i ) {
        for ( int32 j=0; j<height; ++j ) {
          src_image(i,j,nn) = (*inputArrays[nn])[j][i];
        }
      } 
    }

    ImageBuffer src = src_image.buffer();
    convert( dest, src );
        
    // Clean up
    for ( int32 nn = 0; nn < m_format.planes; nn++) {
      delete inputArrays[nn];
    }
    
  } catch (Iex::BaseExc e) {
    vw_throw( vw::IOErr() << "Failed to open " << m_filename << " using the OpenEXR image reader.\n\t" << e.what() );
  } 
}

// Write the given buffer into the disk image.
void vw::DiskImageResourceOpenEXR::write( ImageBuffer const& src, BBox2i const& bbox )
{
  if (!m_output_file_ptr) 
    vw_throw( LogicErr() << "DiskImageResourceOpenEXR: Could not write file. No file has been opened." );

  // This is pretty simple since we always write 32-bit floating point
  // files.  Note that we handle multi-channel images with interleaved
  // planes.  We've already ensured that either planes==1 or
  // channels==1.
  ImageView<float> openexr_image_block( bbox.width(), bbox.height(), m_format.planes );
  ImageBuffer dst = openexr_image_block.buffer();
  convert( dst, src );
  
  try {      
    Imf::FrameBuffer frameBuffer;
        
    // Build the framebuffer out of the various image channels 
    for (int32 nn = 0; nn < dst.format.planes; nn++) {
      frameBuffer.insert (m_labels[nn].c_str(), 
                          Imf::Slice (Imf::FLOAT, (char*) (&(openexr_image_block(-bbox.min()[0],-bbox.min()[1],nn))), 
                                      sizeof (openexr_image_block(0,0,nn)) * 1, 
                                      sizeof (openexr_image_block(0,0,nn)) * dst.format.cols));             
    } 
    
    // Write the data to disk 
    m_output_file_ptr->setFrameBuffer (frameBuffer);
    //    std::cout << "Writing --> " << m_output_file_ptr->currentScanLine() << "\n";
    m_output_file_ptr->writePixels (bbox.height());
    
  } catch (Iex::BaseExc e) {
    vw_throw( vw::IOErr() << "DiskImageResourceOpenEXR: Failed to write " << m_filename << ".\n\t" << e.what() );
  }
  
}

// A FileIO hook to open a file for reading
vw::DiskImageResource* vw::DiskImageResourceOpenEXR::construct_open( std::string const& filename ) {
  return new DiskImageResourceOpenEXR( filename );
}

// A FileIO hook to open a file for writing
vw::DiskImageResource* vw::DiskImageResourceOpenEXR::construct_create( std::string const& filename,
                                                                       ImageFormat const& format ) {
  return new DiskImageResourceOpenEXR( filename, format );
}
