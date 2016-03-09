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


/// \file DiskImageResourceRaw.h
///
/// Provides support for the any sort of raw image data on disk.
#ifndef __VW_FILEIO_DiskImageResourceRaw_H__
#define __VW_FILEIO_DiskImageResourceRaw_H__

#include <string>
#include <fstream>
#include <boost/shared_ptr.hpp>

#include <vw/FileIO/DiskImageResource.h>

namespace vw {

  /// Provides support for the any sort of raw image data on disk.
  /// - Currently this class only supports single channel images.
  /// - Currently the DiskImageResource internal factory function 
  ///   will FAIL to load this file type because the interface does
  ///   not allow the the header info to be passed in.  Until there is
  ///   an interface change this file type needs to be specifically loaded!
  class DiskImageResourceRaw : public DiskImageResource {
  public:

    /// The file contains no header info so the user must provide it.
    DiskImageResourceRaw( std::string const& filename,
                          ImageFormat const& format,
                          bool read_only = true,
                          Vector2i const& block_size=Vector2i(-1,-1))
      : DiskImageResource(filename) { 
      open(filename, format, read_only, block_size); 
    }

    virtual ~DiskImageResourceRaw() {close();}

    // Returns the type of disk image resource.
    static std::string type_static(){ return "Raw"; }

    // Returns the type of disk image resource.
    virtual std::string type() { return type_static(); }

    /// Read the image resource at the given location into the given buffer.
    /// - Though the ImageBuffer object is const, the contents of the buffer will change!
    virtual void read (ImageBuffer const& dest,   BBox2i const& bbox) const;

    /// Write the given buffer to the image resource at the given location.
    virtual void write(ImageBuffer const& source, BBox2i const& bbox);
    
    
    virtual void flush() {m_stream.flush();}

    /// Bind the resource to a file for reading and/or writing.
    void open( std::string const& filename,
               ImageFormat const& format,
               bool read_only = true,
               Vector2i const& block_size=Vector2i(-1,-1));

                 
    /// Closes any open file.
    void close();

    /// Factory function
    static DiskImageResource* construct( std::string const& filename,
                                         ImageFormat const& format,
                                         bool read_only = true,
                                         Vector2i const& block_size=Vector2i(-1,-1) ){
      return new DiskImageResourceRaw(filename, format, read_only, block_size);
    }
    
    // Factory functions required by DiskImageResource.cc
    static DiskImageResource* construct_open( std::string const& filename ) {
      vw_throw( NoImplErr() << "Cannot construct DiskImageResourceRaw without header information!");
      return 0;
    }
    static DiskImageResource* construct_create( std::string const& filename,
                                                ImageFormat const& format ) {
      return new DiskImageResourceRaw(filename, format);
    }


    // This format supports block read/write, but not nodata.

    virtual bool has_block_write () const {return true; }
    virtual bool has_nodata_write() const {return false;}
    virtual bool has_block_read  () const {return true; }
    virtual bool has_nodata_read () const {return false;}

    /// Returns the preferred block size/alignment for partial reads.
    virtual Vector2i block_read_size() const { return m_block_size; }

    /// Gets the preferred block size/alignment for partial writes.
    virtual Vector2i block_write_size() const { return m_block_size; }

    /// Sets the preferred block size/alignment for partial writes.
    virtual void set_block_write_size(const Vector2i& block_size);


  private:
  
    /// Throws an exception if the format is bad
    void check_format() const;
  
    mutable std::fstream m_stream;
    Vector2i m_block_size;
  };

} // namespace VW

#endif//__VW_FILEIO_DiskImageResourceRaw_H__
