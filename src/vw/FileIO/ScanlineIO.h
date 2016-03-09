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


#ifndef __VW_FILEIO_SCANLINEIO_H__
#define __VW_FILEIO_SCANLINEIO_H__

#include <vw/Image/ImageResource.h>

/// \file ScanlineIO.h Abstract scanline reading interface.

namespace vw {
namespace fileio {
namespace detail {

/// Abstract scanline reading interface.
/// - This interface is common to several file types such as PNG and JPEG.
class ScanlineBackend {
  protected:
    // The subclass constructor must initialize these
    ImageFormat m_fmt;
    size_t      m_cstride;

    // attach a data source/sink
    virtual void bind() = 0;

  public:
    virtual ~ScanlineBackend() {}

    size_t chan_bytes() const;
    const ImageFormat& fmt() const;

    // After this, protected members should be initialized,  as well as the
    // data source (i.e., your implementation should call bind())
    virtual void open() = 0;

    // Best guess: Can we perform a read or write? If the return value is
    // false, the object should be destroyed.
    virtual bool ready() const = 0;
};

class ScanlineReadBackend : public ScanlineBackend {
  protected:
    // The subclass constructor must initialize these
    size_t m_rstride;
  public:
    size_t line_bytes() const;
    // decode data into buffer. buffer must be able to hold line_bytes() *
    // lines * sizeof(uint8) bytes
    virtual void read(uint8* buffer, size_t bufsize) = 0;

    // Rewind to allow another read. Might be a noop if it's not necessary.
    virtual ScanlineReadBackend* rewind() const = 0;
};

class ScanlineWriteBackend : public ScanlineBackend {
  public:
    // encode data from buffer. buffer must hold rows * cols * chan_bytes() bytes
    virtual void write(const uint8* buffer, size_t bufsize, size_t rows, size_t cols, size_t planes) = 0;
    virtual const uint8* data() const = 0;
    virtual size_t size() const = 0;
};

}}} // namespace vw::fileio::detail

#endif
