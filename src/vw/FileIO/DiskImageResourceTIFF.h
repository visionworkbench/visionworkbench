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


/// \file FileIO/DiskImageResourceTIFF.h
///
/// Provides support for TIFF image files.
///
#ifndef __VW_FILEIO_DISKIMAGERESOUCETIFF_H__
#define __VW_FILEIO_DISKIMAGERESOUCETIFF_H__

#include <string>

#include <vw/FileIO/DiskImageResource.h>

namespace vw {

  class DiskImageResourceInfoTIFF;

  class VW_API DiskImageResourceTIFF : public DiskImageResource {
  public:

    DiskImageResourceTIFF( std::string const& filename );

    DiskImageResourceTIFF( std::string const& filename,
                           ImageFormat const& format,
                           bool use_compression = false );

    virtual ~DiskImageResourceTIFF() {}

    /// Returns the type of disk image resource.
    static std::string type_static() { return "TIFF"; }

    /// Returns the type of disk image resource.
    virtual std::string type() { return type_static(); }

    virtual bool has_block_write()  const {return false;}
    virtual bool has_nodata_write() const {return false;}
    virtual bool has_block_read()   const {return true;}
    virtual bool has_nodata_read()  const {return false;}

    virtual Vector2i block_read_size() const;

    virtual void read( ImageBuffer const& buf, BBox2i const& bbox ) const;

    virtual void write( ImageBuffer const& dest, BBox2i const& bbox );

    void open( std::string const& filename );

    void create( std::string const& filename,
                 ImageFormat const& format );

    static DiskImageResource* construct_open( std::string const& filename );

    static DiskImageResource* construct_create( std::string const& filename,
                                                ImageFormat const& format );

    void use_lzw_compression(bool state) { m_use_compression = state; }

  protected:
    void check_retval(const int retval, const int error_val) const;

  private:
    boost::shared_ptr<DiskImageResourceInfoTIFF> m_info;
    bool m_use_compression;
  };

} // namespace vw

#endif // __VW_FILEIO_DISK_IMAGE_RESOUCE_TIFF_H__
