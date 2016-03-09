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


/// \file DiskImageResourcePBM.h
///
/// Provides support for the Netpbm format.
///
/// Supported File Types
/// PBM - Monochrome - P1 (means in ASCII) - P4 (means in Binary)
/// PGM - Grayscale  - P2 (ASCII) - P5 (Binary)
/// PPM - RGB Color  - P3 (ASCII) - P5 (Binary)
#ifndef __VW_FILEIO_DISKIMAGERESOURCEPBM_H__
#define __VW_FILEIO_DISKIMAGERESOURCEPBM_H__

#include <string>
#include <fstream>
#include <boost/shared_ptr.hpp>

#include <vw/FileIO/DiskImageResource.h>

namespace vw {

  /// Provides support for the Netpbm format.
  class DiskImageResourcePBM : public DiskImageResource {
  public:

    // Standard DiskImageResource interface:

    DiskImageResourcePBM( std::string const& filename ); // Reading

    DiskImageResourcePBM( std::string const& filename,
                          ImageFormat const& format ); // Writing

    virtual ~DiskImageResourcePBM() {}

    // Returns the type of disk image resource.
    static std::string type_static(){ return "PBM"; }

    // Returns the type of disk image resource.
    virtual std::string type() { return type_static(); }

    virtual void read (ImageBuffer const& buf,  BBox2i const& bbox) const;
    virtual void write(ImageBuffer const& dest, BBox2i const& bbox);
    virtual void flush() {}

    /// Bind the resource to a file for reading.
    void open( std::string const& filename );

    /// Bind the resource to a file for writing.
    void create( std::string const& filename,
                 ImageFormat const& format );

    // Factory functions
    static DiskImageResource* construct_open( std::string const& filename ) {
      return new DiskImageResourcePBM( filename );
    }
    static DiskImageResource* construct_create( std::string const& filename,
                                                ImageFormat const& format ) {
      return new DiskImageResourcePBM( filename, format );
    }

    // Use ascii modes instead of binary ones? (default is binary)
    static void default_to_ascii(bool ascii);

    virtual bool has_block_write () const {return false;}
    virtual bool has_nodata_write() const {return false;}
    virtual bool has_block_read  () const {return false;}
    virtual bool has_nodata_read () const {return false;}

  private:
    std::streampos m_image_data_position;
    std::string m_magic;
    int32 m_max_value;
  };

} // namespace VW

#endif//__VW_FILEIO_DISKIMAGERESOURCEPBM_H__
