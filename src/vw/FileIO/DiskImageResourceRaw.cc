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


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <vw/Core/Exception.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Math/BBox.h>
#include <vw/FileIO/DiskImageResourceRaw.h>

#include <fstream>

#include <boost/smart_ptr/scoped_array.hpp>
#include <boost/algorithm/string/predicate.hpp>

using std::fstream;
using std::ifstream;
using std::ofstream;


namespace vw {

void DiskImageResourceRaw::close() {
  m_stream.close();
  m_format.cols = 0;
  m_format.rows = 0;
}

void DiskImageResourceRaw::check_format() const {
  if (m_format.planes != 1)
    vw_throw( NoImplErr() << "DiskImageResourceRaw does not support multi-plane images.");
  if (num_channels(m_format.pixel_format) != 1)
    vw_throw( NoImplErr() << "DiskImageResourceRaw does not yet support multi-channel images.");
  if ((m_format.cols < 1) || (m_format.rows < 1))
    vw_throw( ArgumentErr() << "Input image is size zero!");    
}

void DiskImageResourceRaw::set_block_write_size(const Vector2i& block_size) {

  if ( (block_size[0] > cols()) || (block_size[1] > rows()) )
    vw_throw( ArgumentErr() << "Requested block size is too big!");

  if ( (block_size[0] < 0) || (block_size[1] < 0) ){ // Use a default size
    if (rows() < 1024)
      m_block_size = Vector2i(cols(), rows());
    else
      m_block_size = Vector2i(cols(), 1024);
  }
  else
    m_block_size = block_size;
}

void DiskImageResourceRaw::open(std::string const& filename, ImageFormat const& format,
                                bool read_only, Vector2i const& block_size) {
  // Update internal variables
  close();
  m_format   = format;
  m_filename = filename;
  check_format();
  
  set_block_write_size(block_size);

  // Open the input stream
  if (read_only)
    m_stream.open(filename.c_str(), fstream::in|fstream::binary);
  else
    m_stream.open(filename.c_str(), fstream::in|fstream::out|fstream::binary);
  if (!m_stream.is_open())
    vw_throw( vw::ArgumentErr() << "DiskImageResourceRaw: Failed to open \"" << filename << "\"." );
}

void DiskImageResourceRaw::read( ImageBuffer const& dest, BBox2i const& bbox )  const {

  // Check that the input parameters are valid
  VW_ASSERT( (bbox.width ()<=cols()) &&
             (bbox.height()<=rows()),
             IOErr() << "Requested read bbox is out of bounds." );
  VW_ASSERT( (static_cast<int>(dest.format.cols)>=bbox.width ()) &&
             (static_cast<int>(dest.format.rows)>=bbox.height()),
             IOErr() << "Buffer is too small for requested read bbox." );

  // Compute the raw data positions (in bytes) in the file we need to read.
  // - For now we only support a single channel so it is pretty simple.
  std::streamsize read_width = m_format.cstride() * bbox.width();
  std::streamsize stride     = m_format.rstride();
  std::streampos  offset     = bbox.min().y()*stride + bbox.min().x()*m_format.cstride();
  std::streamsize total_size = read_width * bbox.height();
  
  // TODO: Optimize speed by allowing a read directly into a compatible buffer
  //       instead of going through a temporary buffer.
  //// Read all the data directly into the output buffer
  //if ((dest.pixel_format() == m_format.pixel_format) &&
  //    (dest.channel_type() == m_format.channel_type) &&
  //    (dest.format.
  //   )
  

  // Create a temporary image buffer just big enough to contain the input data.  
  boost::scoped_array<uint8> image_data(new uint8[total_size]);
  ImageFormat buffer_format = m_format;
  buffer_format.cols = bbox.width();
  buffer_format.rows = bbox.height();
  ImageBuffer source(buffer_format, image_data.get());
  
  // Read all the data into a temporary buffer
  char* write_ptr = reinterpret_cast<char*>(image_data.get());
  for (int i=0; i<bbox.height(); ++i) { // For each line
    m_stream.seekg(offset + i*stride);     // Move to start of read line
    m_stream.read(write_ptr, read_width);  // Read the line
    write_ptr += read_width;               // Move to the next write line
  }

  // Copy the data from the temporary buffer to the output buffer,
  //  performing necessary format conversions.
  convert(dest, source, false); 
}

void DiskImageResourceRaw::write( ImageBuffer const& source, BBox2i const& bbox ) {
  // Check that the input parameters are valid
  VW_ASSERT( (bbox.width ()<=cols()) &&
             (bbox.height()<=rows()),
             IOErr() << "Requested read bbox is out of bounds." );
  VW_ASSERT( (static_cast<int>(source.format.cols)>=bbox.width ()) &&
             (static_cast<int>(source.format.rows)>=bbox.height()),
             IOErr() << "Buffer is too small for requested read bbox." );

  // Compute the raw data positions (in bytes) in the file we need to read.
  // - For now we only support a single channel so it is pretty simple.
  std::streamsize read_width = m_format.cstride() * bbox.width();
  std::streamsize stride     = m_format.rstride();
  std::streampos  offset     = bbox.min().y()*stride + bbox.min().x()*m_format.cstride();
  std::streamsize total_size = read_width * bbox.height();
  
  // TODO: Optimize speed by allowing a read directly into a compatible buffer
  //       instead of going through a temporary buffer.
  //// Read all the data directly into the output buffer
  //if ((dest.pixel_format() == m_format.pixel_format) &&
  //    (dest.channel_type() == m_format.channel_type) &&
  //    (dest.format.
  //   )

  // Create a temporary image buffer just big enough to contain the input data.  
  boost::scoped_array<uint8> image_data(new uint8[total_size]);
  ImageFormat buffer_format = m_format;
  buffer_format.cols = bbox.width();
  buffer_format.rows = bbox.height();
  ImageBuffer dest(buffer_format, image_data.get());
  
  // Copy the data from the input buffer to the temporary buffer,
  //  performing necessary format conversions.
  convert(dest, source, false);
  
  // Write all the data from the temporary buffer
  char* read_ptr = reinterpret_cast<char*>(image_data.get());
  for (int i=0; i<bbox.height(); ++i) { // For each line
    m_stream.seekg(offset + i*stride);     // Move to start of write line
    m_stream.write(read_ptr, read_width);  // write the line
    read_ptr += read_width;                // Move to the next read line
  }
}

} // end namespace vw



