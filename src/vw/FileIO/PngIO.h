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


#ifndef __VW_FILEIO_PNGIO_H__
#define __VW_FILEIO_PNGIO_H__

#include <vw/FileIO/ScanlineIO.h>

extern "C" {
#include <png.h>
}

/// \file PngIO.h Shares code between the on-disk and in-memory PNG code.


namespace vw {
namespace fileio {
namespace detail {

// These classes exist to share code between the on-disk and in-memory versions
// of the relevant image resources. They are not intended for use by users
// (thus the detail namespace).

class PngIO {
  protected:
    png_structp m_ctx;
    png_infop m_info;
    PngIO();
};

class PngIODecompress : public PngIO, public ScanlineReadBackend {
  private:
    bool m_read;
  public:
    PngIODecompress();
    virtual ~PngIODecompress();

    void open();
    bool ready() const;
    void read(uint8* data, size_t bufsize);
};

class PngIOCompress : public PngIO, public ScanlineWriteBackend {
  private:
    bool m_written;
  public:
    // cols/rows/planes ignored in imageformat
    PngIOCompress(const ImageFormat& fmt);
    virtual ~PngIOCompress();

    void open();
    bool ready() const;
    void write(const uint8* buffer, size_t bufsize, size_t rows, size_t cols, size_t planes);
};

}}} // namespace vw::fileio::detail

#endif
