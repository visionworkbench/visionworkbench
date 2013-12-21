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


#ifndef __VW_FILEIO_GDALIO_H__
#define __VW_FILEIO_GDALIO_H__

#include <vw/Core/Features.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/FileIO/ScanlineIO.h>

/// \file GdalIO.h Shares code between the on-disk and in-memory GDAL code.
class GDALDataset;
class GDALDriver;

namespace vw {
  class Mutex;

namespace fileio {
namespace detail {

Mutex& gdal() VW_WARN_UNUSED;

// These classes exist to share code between the on-disk and in-memory versions
// of the relevant image resources. They are not intended for use by users
// (thus the detail namespace).

class GdalIODecompress : public ScanlineReadBackend {
  protected:
    boost::shared_ptr<GDALDataset> m_dataset;
  public:
    GdalIODecompress();
    virtual ~GdalIODecompress();

    void open();
    bool ready() const;
    void read(uint8* data, size_t bufsize);
    void read(const ImageFormat& fmt, uint8* data, size_t bufsize);

    bool nodata_read_ok(double& value) const;
};

class GdalIOCompress : public ScanlineWriteBackend {
  protected:
    std::string m_fn;
    GDALDriver *m_driver;
    boost::shared_ptr<GDALDataset> m_dataset;
    double m_nodata;
    bool   m_has_nodata;

  public:
    // cols/rows/planes ignored in imageformat
    GdalIOCompress(const ImageFormat& fmt);
    virtual ~GdalIOCompress();

    void open();
    bool ready() const;
    void write(const uint8* buffer, size_t bufsize, size_t rows, size_t cols, size_t planes);
    void set_nodata(double value) {m_nodata = value; m_has_nodata = true;}
};

}}} // namespace vw::fileio::detail

#endif
