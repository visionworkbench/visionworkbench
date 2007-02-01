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
  if (m_file_ptr) 
    delete m_file_ptr;
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
    if (m_file_ptr) 
      vw_throw( IOErr() << "Disk image resources do not yet support reuse." );
    
    m_file_ptr = new Imf::InputFile(filename.c_str());
    Imf::FrameBuffer frameBuffer;
    
    // Find the width and height of the image 
    Imath::Box2i dw = m_file_ptr->header().dataWindow();
    m_format.cols  = int(dw.max.x - dw.min.x + 1);
    m_format.rows  = int(dw.max.y - dw.min.y + 1);
    
    // Determine the number of image channels 
    Imf::ChannelList::ConstIterator iter = m_file_ptr->header().channels().begin();
    int num_channels = 0;
    while( iter++ != m_file_ptr->header().channels().end() )
      num_channels++;
    m_format.planes = num_channels;
    
    // For now, we only support reading in multi-plane, single channel
    // images.
    m_format.pixel_format = VW_PIXEL_SCALAR;
    
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
}

// Read the disk image into the given buffer.
void vw::DiskImageResourceOpenEXR::read( ImageBuffer const& dest, BBox2i const& bbox ) const
{
  VW_ASSERT( bbox.width()==int(cols()) && bbox.height()==int(rows()),
             NoImplErr() << "DiskImageResourceOpenEXR does not support partial reads." );
  VW_ASSERT( dest.format.cols==cols() && dest.format.rows==rows(),
             IOErr() << "Buffer has wrong dimensions in OpenEXR read." );
  
  if (!m_file_ptr) 
    vw_throw( LogicErr() << "DiskImageResourceOpenEXR: Could not read file. No file has been opened." );
  
  try {
    // Find the width and height of the image 
    Imath::Box2i dw = m_file_ptr->header().dataWindow();
    unsigned height = m_format.rows;
    unsigned width = m_format.cols;
    
    // Set up the OpenEXR data structures necessary to read all of
    // the channels out of the file and execute the call to readPixels().
    Imf::Array2D<float> *inputArrays[m_format.planes];
    Imf::ChannelList::ConstIterator channel = m_file_ptr->header().channels().begin();
    std::vector<std::string> channel_names(m_format.planes);
    for (int i=0; channel != m_file_ptr->header().channels().end(); ++channel, ++i) {
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
        channel_names[2] = "A";
      }
    }
    
    Imf::FrameBuffer frameBuffer;
    for ( unsigned nn = 0; nn < m_format.planes; ++nn ) {
      inputArrays[nn] = new Imf::Array2D<float>(height,width);
      //        std::cout << "Reading channel " << channel_names[nn] << "\n";
      frameBuffer.insert (channel_names[nn].c_str(), Imf::Slice (Imf::FLOAT, (char *) (&(*inputArrays[nn])[0][0]), 
                                                                 sizeof ((*inputArrays[nn])[0][0]) * 1, 
                                                                 sizeof ((*inputArrays[nn])[0][0]) * width, 1, 1, 0.0)); 
    }
    m_file_ptr->setFrameBuffer (frameBuffer);
    m_file_ptr->readPixels (dw.min.y, dw.max.y);
    
    // Copy the pixels over into a ImageView object.
    // 
    // Recast to the templatized pixel type in the process.
    ImageView<float> src_image(m_format.cols, m_format.rows, m_format.planes);
    for ( unsigned nn=0; nn<m_format.planes; ++nn ) {
      for ( unsigned i=0; i<width; ++i ) {
        for ( unsigned j=0; j<height; ++j ) {
          src_image(i,j,nn) = (*inputArrays[nn])[j][i];
        }
      } 
    }
    ImageBuffer src = src_image.buffer();
    convert( dest, src );
    
    // Print out the image size and number of channels
    std::cout << "OpenEXR file " << m_filename
              << "\t" << m_format.planes << " x " << width << " x " << height << "\n";
    
    // Clean up
    for ( unsigned nn = 0; nn < m_format.planes; nn++) {
      delete inputArrays[nn];
    }
    
  } catch (Iex::BaseExc e) {
    vw_throw( vw::IOErr() << "Failed to open " << m_filename << " using the OpenEXR image reader.\n\t" << e.what() );
  } 
}

// Write the given buffer into the disk image.
void vw::DiskImageResourceOpenEXR::write( ImageBuffer const& src, BBox2i const& bbox )
{
  VW_ASSERT( bbox.width()==int(cols()) && bbox.height()==int(rows()),
             NoImplErr() << "DiskImageResourceOpenEXR does not support partial writes." );
  VW_ASSERT( src.format.cols==cols() && src.format.rows==rows(),
             IOErr() << "Buffer has wrong dimensions in OpenEXR write." );
  
  // This is pretty simple since we always write 8-bit integer files.
  // Note that we handle multi-channel images with interleaved planes. 
  // We've already ensured that either planes==1 or channels==1.
  ImageView<float> openexr_image( m_format.cols, m_format.rows, m_format.planes );
  ImageBuffer dst = openexr_image.buffer();
  convert( dst, src );
  
  float* pixels[dst.format.planes];
  Imf::Array2D<float> *floatArrays[dst.format.planes];
  std::string labels[dst.format.planes];
  
  try {      
    // Create the file header with the appropriate number of
    // channels.  Label the channels in order starting with "Channel 0".
    Imf::Header header (dst.format.cols,dst.format.rows);
    for ( unsigned nn = 0; nn < dst.format.planes; nn++) {
      labels[nn] = openexr_channel_string_of_pixel_type(m_format.pixel_format, nn);
      //        std::cout << "Writing channel " << nn << ": " << labels[nn] << "\n";
      header.channels().insert (labels[nn].c_str(), Imf::Channel (Imf::FLOAT));
      floatArrays[nn] = new Imf::Array2D<float>(dst.format.rows,dst.format.cols);
    }
    
    // Open the file handle and create an empty framebuffer object. 
    Imf::OutputFile file (m_filename.c_str(), header);
    Imf::FrameBuffer frameBuffer;
    
    // Copy the actual data into temporary memory, which will
    // ultimately be written to the file.
    for ( unsigned nn = 0; nn < dst.format.planes; nn++) 
      for ( unsigned i = 0; i < dst.format.cols; i++) 
        for ( unsigned j = 0; j < dst.format.rows; j++) 
          (*floatArrays[nn])[j][i] = openexr_image(i,j,nn);
    
    // Build the framebuffer out of the various image channels 
    for (unsigned int nn = 0; nn < dst.format.planes; nn++) {
      pixels[nn] = &((*floatArrays[nn])[0][0]);
      frameBuffer.insert (labels[nn].c_str(), 
                          Imf::Slice (Imf::FLOAT, (char*) pixels[nn], 
                                      sizeof (*(pixels[nn])) * 1, 
                                      sizeof (*(pixels[nn])) * dst.format.cols));             
    } 
    
    // Write the data to disk 
    file.setFrameBuffer (frameBuffer);
    file.writePixels (dst.format.rows);
    
    // Clean up 
    for ( unsigned nn = 0; nn < dst.format.planes; ++nn )
      delete floatArrays[nn];
    
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
