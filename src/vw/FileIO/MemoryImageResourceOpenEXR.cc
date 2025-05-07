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

#include <vw/config.h>
#if defined(VW_HAVE_PKG_OPENEXR) && VW_HAVE_PKG_OPENEXR==1

#include <vw/FileIO/MemoryImageResourceOpenEXR.h>
#include <vw/Core/Debugging.h>

#include <iostream>
#include <string.h>
#include <algorithm>
#include <string>
#include <vector>

#include <ImathBox.h>
#include <ImathInt64.h>
#include <ImathVec.h>
#include <ImfFrameBuffer.h>
#include <ImfHeader.h>
#include <ImfPixelType.h>
#include <ImfInputFile.h>
#include <ImfOutputFile.h>
#include <ImfIO.h>
#include <ImfChannelList.h>


using namespace vw;

Imf::PixelType vw_channel_id_to_imf_pix_fmt( ChannelTypeEnum vw_type) {
  switch (vw_type) {
  case VW_CHANNEL_UINT32:
    return Imf::UINT;
  case VW_CHANNEL_FLOAT32:
  default:
    return Imf::FLOAT;
  }
}

class SrcMemoryImageResourceOpenEXR::Data : public Imf::IStream, public std::istringstream {
public:
  Data( std::string const& str ) : IStream(""), std::istringstream(str) {}

  virtual bool isMemoryMapped() const { return false; }
  virtual bool read( char c[], int n ) {
    std::istringstream::read( c, n );
    return std::istringstream::good();
  }
  virtual char* readMemoryMapped( int /*n*/ ) {
    vw_throw( NoImplErr() << "readMemoryMapped should never be called." );
    return NULL;
  }
  virtual Imf::Int64 tellg() {
    return std::istringstream::tellg();
  }
  virtual void seekg( Imf::Int64 pos ) {
    std::istringstream::seekg( pos );
  }
  virtual void clear() {
    std::istringstream::clear();
  }

  ImageFormat format;
};

SrcMemoryImageResourceOpenEXR::SrcMemoryImageResourceOpenEXR( boost::shared_array<const uint8> buffer, size_t len ) {
  VW_ASSERT(buffer, ArgumentErr() << VW_CURRENT_FUNCTION << ": buffer must be non-null");
  VW_ASSERT(len,    ArgumentErr() << VW_CURRENT_FUNCTION << ": len must be non-zero");

  m_data.reset( new Data(std::string((const char*)buffer.get(), len)) );

  Imf::InputFile file( *m_data.get() );
  Imath::Box2i dw = file.header().dataWindow();
  m_data->format.cols = dw.max.x - dw.min.x + 1;
  m_data->format.rows = dw.max.y - dw.min.y + 1;

  // We only support multi-plane images. Everything will have to be
  // converted to the right type by vw::Convert.
  m_data->format.planes = 0;
  for ( Imf::ChannelList::ConstIterator iter = file.header().channels().begin();
        iter != file.header().channels().end(); iter++ )
    m_data->format.planes++;

  m_data->format.pixel_format = VW_PIXEL_SCALAR;
  switch ( file.header().channels().begin().channel().type ) {
  case Imf::UINT:
    m_data->format.channel_type = VW_CHANNEL_UINT32;
  case Imf::HALF:
  case Imf::FLOAT:
  default:
    m_data->format.channel_type = VW_CHANNEL_FLOAT32;
  }
}

void SrcMemoryImageResourceOpenEXR::read( ImageBuffer const& dst, BBox2i const& bbox ) const {
  VW_ASSERT( dst.format.cols == size_t(bbox.width()) && dst.format.rows == size_t(bbox.height()),
             ArgumentErr() << VW_CURRENT_FUNCTION << ": Destination buffer has wrong dimensions!" );
  VW_ASSERT( dst.format.cols == size_t(cols()) && dst.format.rows == size_t(rows()),
             ArgumentErr() << VW_CURRENT_FUNCTION << ": Partial reads are not supported");

  m_data->seekg( 0 );
  Imf::InputFile file( *m_data.get() );

  // Set up the OpenEXR data structures necessary to read all of
  // the channels out of the file and execute the call to readPixels().
  std::vector<std::string> channel_names;
  typedef Imf::ChannelList::ConstIterator ChIterT;
  for ( ChIterT i = file.header().channels().begin();
        i != file.header().channels().end(); i++ )
    channel_names.push_back( i.name() );

  // OpenEXR seems to order channels in the file alphabetically,
  // rather than in the order in which they were saved.  This means
  // that we need to reorder the channel names when they are
  // labelled as RGB or RGBA.  For other channel naming schemes, we
  // just go with alphabetical, since that's all we've got.
  if ( m_data->format.planes == 3 ) {
    if (find(channel_names.begin(), channel_names.end(), "R") != channel_names.end() &&
        find(channel_names.begin(), channel_names.end(), "G") != channel_names.end() &&
        find(channel_names.begin(), channel_names.end(), "B") != channel_names.end()) {
      channel_names[0] = "R";
      channel_names[1] = "G";
      channel_names[2] = "B";
    }
  } else if ( m_data->format.planes == 4 ) {
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

  // Allocate the memory we need
  boost::scoped_array<char> src_decoded_data(new char[m_data->format.byte_size()]);
  ImageBuffer src(m_data->format, src_decoded_data.get());

  // Work out offsets in the source data
  Imath::Box2i dw = file.header().dataWindow();

  // Attach memory to OpenEXR device
  Imf::FrameBuffer frameBuffer;
  for ( size_t nn = 0; nn < m_data->format.planes; ++nn ) {
    frameBuffer.insert( channel_names[nn].c_str(),
                        Imf::Slice( vw_channel_id_to_imf_pix_fmt( m_data->format.channel_type ),
                                    src_decoded_data.get() + nn * src.pstride
                                    - dw.min.y * src.rstride - dw.min.x * src.cstride,
                                    src.cstride, src.rstride ) );
  }

  // Attach buffer device to file and read
  file.setFrameBuffer( frameBuffer );
  file.readPixels( dw.min.y, dw.max.y );

  convert( dst, src, true );
}

ImageFormat SrcMemoryImageResourceOpenEXR::format() const {
  return m_data->format;
}

// This could be std::ostringstream, but I couldn't figure out a way
// to get a pointer to the persistant memory underneath the
// structure. str() provides a temporary copy.
class DstMemoryImageResourceOpenEXR::Data : public Imf::OStream {
  boost::shared_array<char> m_array;
  size_t m_perceived_size, m_actual_size, m_index;

  // Attempts to double current memory usage or meet the minimum
  // demand.
  void increase_size( size_t min_size ) {
    size_t new_size = m_actual_size * 2 < min_size ? min_size : m_actual_size * 2;
    boost::shared_array<char> new_array( new char[new_size] );
    m_actual_size = new_size;
    std::memcpy( new_array.get(), m_array.get(), m_perceived_size );
    m_array.swap( new_array );
  }

public:
  Data(ImageFormat const& fmt) : OStream(""), m_array( new char[512] ),
                                 m_perceived_size(0), m_actual_size(512),
                                 m_index(0), format(fmt) {}

  char* data() const { return m_array.get(); }

  virtual void write( const char c[], int n ) {
    if ( m_index + n >= m_actual_size )
      increase_size(m_index+n);
    std::memcpy( m_array.get() + m_index, c, n );
    m_index += n;
    m_perceived_size += n;
  }

  virtual Imf::Int64 tellp() {
    return m_index;
  }

  virtual void seekp( Imf::Int64 pos ) {
    if ( pos >= m_actual_size )
      increase_size(pos + 1);
    m_index = size_t( pos );
  }

  size_t size() const { return m_perceived_size; }
  size_t actual_size() const { return m_actual_size; }

  ImageFormat format;
  std::vector<std::string> labels;
};

DstMemoryImageResourceOpenEXR::DstMemoryImageResourceOpenEXR( const ImageFormat& fmt ) : m_data( new Data(fmt) ) {
  VW_ASSERT( fmt.planes == 1 || fmt.pixel_format == VW_PIXEL_SCALAR,
             NoImplErr() << VW_CURRENT_FUNCTION << ": Cannot create image with both multiple planes and multiple channels.\n");

  m_data->format.planes =
    std::max( fmt.planes, num_channels(fmt.pixel_format) );
  m_data->labels.resize( m_data->format.planes );

  if ( fmt.pixel_format == VW_PIXEL_RGB ) {
    m_data->labels[0] = "R";
    m_data->labels[1] = "G";
    m_data->labels[2] = "B";
  } else if ( fmt.pixel_format == VW_PIXEL_RGBA ) {
    m_data->labels[0] = "R";
    m_data->labels[1] = "G";
    m_data->labels[2] = "B";
    m_data->labels[3] = "A";
  } else {
    for ( size_t i = 0; i < m_data->format.planes; i++ ) {
      std::ostringstream stream;
      stream << "Channel" << i;
      m_data->labels[i] = stream.str();
    }
  }

  switch ( fmt.channel_type ) {
  case VW_CHANNEL_UINT32:
    break;
  default:
    m_data->format.channel_type = VW_CHANNEL_FLOAT32;
  }
  m_data->format.pixel_format = VW_PIXEL_SCALAR;
}

void DstMemoryImageResourceOpenEXR::write( ImageBuffer const& src, BBox2i const& bbox ) {
  VW_ASSERT( src.format.cols == size_t(bbox.width()) && src.format.rows == size_t(bbox.height()),
             ArgumentErr() << VW_CURRENT_FUNCTION << ": partial writes not supported." );
  VW_ASSERT( src.format.cols == m_data->format.cols && src.format.rows == m_data->format.rows,
             ArgumentErr() << VW_CURRENT_FUNCTION << ": Source buffer should equal internal format size." );

  // Converting USER's data to be plane format
  boost::scoped_array<char> src_converted_data(new char[m_data->format.byte_size()]);
  ImageBuffer src_converted_buffer(m_data->format, src_converted_data.get());
  convert( src_converted_buffer, src, true );

  // Attaching data to Imf's framebuffer
  Imf::FrameBuffer frameBuffer;
  Imf::Header header( m_data->format.cols, m_data->format.rows );
  for ( size_t nn = 0; nn < m_data->format.planes; nn++ ) {
    header.channels().insert( m_data->labels[nn].c_str(),
                              Imf::Channel( vw_channel_id_to_imf_pix_fmt( m_data->format.channel_type ) ) );
    frameBuffer.insert( m_data->labels[nn].c_str(),
                        Imf::Slice( vw_channel_id_to_imf_pix_fmt( m_data->format.channel_type ),
                                    src_converted_data.get() + nn * src_converted_buffer.pstride,
                                    src_converted_buffer.cstride,
                                    src_converted_buffer.rstride ));
  }

  Imf::OutputFile out( *m_data.get(), header );
  out.setFrameBuffer( frameBuffer );
  out.writePixels( src_converted_buffer.format.rows );
}

const uint8* DstMemoryImageResourceOpenEXR::data() const {
  return (const uint8*)m_data->data();
}

size_t DstMemoryImageResourceOpenEXR::size() const {
  return m_data->size();
}

#endif // VW_HAVE_PKG_OPENEXR
