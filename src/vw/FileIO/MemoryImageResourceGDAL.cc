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


#include <vw/FileIO/MemoryImageResourceGDAL.h>
#include <vw/FileIO/GdalIO.h>
#include <vw/Core/Debugging.h>

#include <boost/format.hpp>
#include <gdal.h>
#include <gdal_priv.h>

namespace {
  std::string make_fn(const char* name, const void* key) {
    static const boost::format GDAL_MEM_FILENAME("/vsimem/vw_%s_%u_%p");
    return std::string(boost::str(boost::format(GDAL_MEM_FILENAME) % name % vw::Thread::id() % key));
  }
  // GDAL warns when GDALClose is called on a NULL, even though it's documented
  // as being okay...
  void GDALCloseNullOk( GDALDatasetH x) {
    if (x)
      ::GDALClose(x);
  }
}

namespace vw {

class SrcMemoryImageResourceGDAL::Data : public fileio::detail::GdalIODecompress {
    boost::shared_array<const uint8> m_data;
    size_t m_len;
  protected:
    virtual void bind() {
      const std::string src_fn(make_fn("src", m_data.get()));
      VSIFCloseL(VSIFileFromMemBuffer(src_fn.c_str(), const_cast<uint8*>(m_data.get()), m_len, false));
      m_dataset.reset(reinterpret_cast<GDALDataset*>(GDALOpen(src_fn.c_str(), GA_ReadOnly)), GDALCloseNullOk);
      if (!m_dataset) {
        VSIUnlink(src_fn.c_str());
        vw_throw(IOErr() << "Unable to open memory dataset.");
      }
    }
    Data* rewind() const { vw_throw(NoImplErr() << VW_CURRENT_FUNCTION << ": not supported"); return NULL; }
  public:
    Data(boost::shared_array<const uint8> buffer, size_t len) : m_data(buffer), m_len(len) {
      VW_ASSERT(buffer, ArgumentErr() << VW_CURRENT_FUNCTION << ": buffer must be non-null");
      VW_ASSERT(len,    ArgumentErr() << VW_CURRENT_FUNCTION << ": len must be non-zero");
    }
};

SrcMemoryImageResourceGDAL::SrcMemoryImageResourceGDAL(boost::shared_array<const uint8> buffer, size_t len)
  : m_data(new Data(buffer, len)) {
    m_data->open();
}

void SrcMemoryImageResourceGDAL::read( ImageBuffer const& dst, BBox2i const& bbox_ ) const {
  size_t width = bbox_.width(), height = bbox_.height(), planes = dst.format.planes;
  VW_ASSERT( dst.format.cols == width && dst.format.rows == height,
             ArgumentErr() << VW_CURRENT_FUNCTION << ": Destination buffer has wrong dimensions!" );
  VW_ASSERT( dst.format.cols == size_t(cols()) && dst.format.rows == size_t(rows()),
             ArgumentErr() << VW_CURRENT_FUNCTION << ": Partial reads are not supported");

  // shared rather than scoped so we get a deleter.
  boost::shared_array<uint8> buf;

  // no conversion necessary?
  bool simple = m_data->fmt().simple_convert(dst.format);

  // If we don't need to convert, we read directly into the dst buffer (using a
  // noop_deleter, so the destructor doesn't try to delete it)
  if (simple) {
    m_data->read(dst.format, reinterpret_cast<uint8*>(dst.data), dst.format.byte_size());
    return;
  } else {
    size_t bufsize = m_data->line_bytes() * height * planes;
    buf.reset( new uint8[bufsize] );
    m_data->read(m_data->fmt(), buf.get(), bufsize);
  }

  ImageFormat src_fmt(m_data->fmt());
  src_fmt.rows = height;
  src_fmt.cols = width;

  ImageBuffer src(src_fmt, buf.get());
  convert(dst, src, true);
}

ImageFormat SrcMemoryImageResourceGDAL::format() const {
  return m_data->fmt();
}

bool SrcMemoryImageResourceGDAL::has_nodata_read() const {
  double value;
  return m_data->nodata_read_ok(value);
}

double SrcMemoryImageResourceGDAL::nodata_read() const {
  double value;
  bool ok = m_data->nodata_read_ok(value);
  VW_ASSERT(ok, IOErr() << VW_CURRENT_FUNCTION << ": This dataset does not have a nodata value.");
  return value;
}

class DstMemoryImageResourceGDAL::Data : public fileio::detail::GdalIOCompress {
  protected:
    virtual void bind() { m_fn = make_fn("dst", this); }
  public:
    Data(const ImageFormat &fmt) : GdalIOCompress(fmt) {}
    ~Data() {
      if (!m_fn.empty())
          VSIUnlink(m_fn.c_str()); // ignore return code, we can't do anything if it failed
    }
    const uint8* data() const { return VSIGetMemFileBuffer( m_fn.c_str(), 0, FALSE); }
    size_t size() const {vsi_l_offset s; VSIGetMemFileBuffer( m_fn.c_str(), &s, FALSE); return s; }
};


DstMemoryImageResourceGDAL::DstMemoryImageResourceGDAL(const ImageFormat& fmt)
  : m_data(new Data(fmt))
{
  m_data->open();
}

void DstMemoryImageResourceGDAL::write( ImageBuffer const& src, BBox2i const& bbox_ ) {
  size_t width = bbox_.width(), height = bbox_.height(), planes = src.format.planes;
  VW_ASSERT( src.format.cols == width && src.format.rows == height,
             ArgumentErr() << VW_CURRENT_FUNCTION << ": partial writes not supported." );
  VW_ASSERT(m_data->ready(), LogicErr() << "Multiple writes to one DstMemoryImageResourceGDAL. Probably a bug?");

  // shared rather than scoped so we get a deleter.
  boost::shared_array<uint8> buf;

  // no conversion necessary?
  bool simple = src.format.simple_convert(m_data->fmt());

  size_t bufsize = m_data->chan_bytes() * width * height * planes;

  // If we don't need to convert, we write directly from the src buffer (using a
  // noop_deleter, so the destructor doesn't try to delete it)
  if (simple)
    buf.reset( reinterpret_cast<uint8*>(const_cast<void*>(src.data)), NOP() );
  else {
    buf.reset( new uint8[bufsize] );

    ImageFormat dst_fmt(m_data->fmt());
    dst_fmt.rows = height;
    dst_fmt.cols = width;

    ImageBuffer dst(dst_fmt, buf.get());
    convert(dst, src, true);
  }

  m_data->write(buf.get(), bufsize, width, height, planes);
}

const uint8* DstMemoryImageResourceGDAL::data() const {
  return m_data->data();
}

size_t DstMemoryImageResourceGDAL::size() const {
  return m_data->size();
}

void DstMemoryImageResourceGDAL::set_nodata_write(double value) {
  m_data->set_nodata(value);
}

} // namespace vw
