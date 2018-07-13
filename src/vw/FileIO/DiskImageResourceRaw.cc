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
#include <boost/filesystem.hpp>

using std::fstream;
using std::ifstream;
using std::ofstream;


namespace fs = boost::filesystem;

namespace vw {

void DiskImageResourceRaw::close() {
  m_stream.close();
  m_format.cols = 0;
  m_format.rows = 0;
}

// Factory functions required by DiskImageResource.cc
DiskImageResource* DiskImageResourceRaw::construct_open( std::string const& filename ) {
  std::string dim_file = find_associated_spot5_dim_file(filename);
  if (dim_file.empty()) 
    vw_throw( ArgumentErr() << "Could not find .DIM file for: " << filename);
  
  return DiskImageResourceRaw::construct(filename, image_format_from_spot5_DIM(dim_file));
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
    std::streampos row_offset = i*stride;
    m_stream.seekg(offset + row_offset);     // Move to start of read line
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
    std::streampos row_offset = i*stride;
    m_stream.seekg(offset + row_offset);     // Move to start of write line
    m_stream.write(read_ptr, read_width);  // write the line
    read_ptr += read_width;                // Move to the next read line
  }
}


std::string DiskImageResourceRaw::find_associated_spot5_dim_file(std::string const& image_file){

  // Try replacing the extension with .DIM
  std::string cam_file = fs::path(image_file).replace_extension(".DIM").string();
  if (fs::exists(cam_file))
    return cam_file;

  // Now try .dim
  cam_file = fs::path(image_file).replace_extension(".dim").string();
  if (fs::exists(cam_file))
    return cam_file;

  // Now try something more complicated.  From
  // back/SEGMT01/IMAGERY.BIL go to back/SEGMT01/METADATA_BACK.DIM or
  // back/SEGMT01/METADATA.DIM
  
  std::string line = image_file;
  std::transform(line.begin(), line.end(), line.begin(), ::tolower); // lowercase
  
  // Try looking in some hard coded folders that the header file is often found in.
  std::size_t found = line.rfind("/imagery.b");
  if (found == std::string::npos)
    return "";
  std::string prefix = image_file.substr(0, found);

  found = line.rfind("front/");
  if (found != std::string::npos) {
    cam_file = prefix + "/METADATA_FRONT.DIM"; if (fs::exists(cam_file)) return cam_file;
    cam_file = prefix + "/metadata_front.dim"; if (fs::exists(cam_file)) return cam_file;
    cam_file = prefix + "/METADATA.DIM";       if (fs::exists(cam_file)) return cam_file;
    cam_file = prefix + "/metadata.dim";       if (fs::exists(cam_file)) return cam_file;
    return "";
  }

  found = line.rfind("back/");
  if (found != std::string::npos) {
    cam_file = prefix + "/METADATA_BACK.DIM";  if (fs::exists(cam_file)) return cam_file;
    cam_file = prefix + "/metadata_back.dim";  if (fs::exists(cam_file)) return cam_file;
    cam_file = prefix + "/METADATA.DIM";       if (fs::exists(cam_file)) return cam_file;
    cam_file = prefix + "/metadata.dim";       if (fs::exists(cam_file)) return cam_file;
    return "";
  }
  
  return "";
}
  
bool DiskImageResourceRaw::parse_int_between_tags(std::string const& line, std::string const& tag, int & val){

  // Keep the output uninitialized, initialize it only on success
  // val = 0;

  // advance to the beginning of <tag>
  std::size_t found = line.find(tag);
  if (found == std::string::npos)
    return false;

  // Advance to the end of <tag>
  std::size_t beg = line.find(">", found);
  if (beg == std::string::npos)
    return false;
  
  // Move past "<"
  beg++;

  // Advance to the beginning of </tag>
  std::size_t end = line.find("<", found);
  if (end == std::string::npos)
    return false;

  std::string val_str = line.substr(beg, end - beg);

  val = atoi(val_str.c_str());
  
  return true;
}

vw::ImageFormat DiskImageResourceRaw::image_format_from_spot5_DIM(std::string const& camera_file) {

  // Find the following text, and read ncols and nrows values
  // <Raster_Dimensions>
  //<NCOLS>12000</NCOLS>
  //<NROWS>46408</NROWS>
  //<NBANDS>1</NBANDS>
  //</Raster_Dimensions>

  vw::Vector2i image_size;
  std::ifstream ifs(camera_file.c_str());
  std::string line;
  while (std::getline(ifs, line)){

    std::transform(line.begin(), line.end(), line.begin(), ::tolower); // lowercase
    std::size_t found = line.find("<raster_dimensions");
    if (found == std::string::npos)
      continue;

    std::string cols_tag = "ncols";
    std::getline(ifs, line);
    std::transform(line.begin(), line.end(), line.begin(), ::tolower); // lowercase
    if (!parse_int_between_tags(line, cols_tag, image_size[0]))
      continue;
    
    std::string rows_tag = "nrows";
    std::getline(ifs, line);
    std::transform(line.begin(), line.end(), line.begin(), ::tolower); // lowercase
    if (!parse_int_between_tags(line, rows_tag, image_size[1]))
      continue;
    
    break;
  }

  if (image_size[0] <= 0 || image_size[1] <= 0)
    vw_throw( ArgumentErr() << "Could not parse correctly the image size from: " << camera_file);
    
  vw::ImageFormat format;
  format.cols          = image_size[0];
  format.rows          = image_size[1];
  format.planes        = 1;
  format.pixel_format  = vw::VW_PIXEL_GRAY; // This should be constant
  format.channel_type  = vw::VW_CHANNEL_UINT8;
  format.premultiplied = true; // Don't do anything funny to the data
  return format;
}


} // end namespace vw



