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


/// \file FileIO/DiskImageResourceOpenEXR.h
///
/// Provides support for the OpenEXR file format.
///
#ifndef __VW_FILEIO_DISKIMAGERESOUCEOPENEXR_H__
#define __VW_FILEIO_DISKIMAGERESOUCEOPENEXR_H__

#include <string>

#include <vw/FileIO/DiskImageResource.h>

namespace vw {

  /// DiskImageResource implementation for the OpenEXR file format.
  /// - OpenEXR is a high dynamic range image file format developed by
  ///   Industrial Light and Magic (ILM).
  class DiskImageResourceOpenEXR : public DiskImageResource {
  protected:

    void set_tiled_write(int32 tile_width, int32 tile_height, bool random_tile_order = false);
    void set_scanline_write(int32 scanlines_per_block);

  public:

    DiskImageResourceOpenEXR( std::string const& filename )
      : DiskImageResource( filename )
    {
      m_input_file_ptr = 0;
      m_output_file_ptr = 0;
      open( filename );
    }

    DiskImageResourceOpenEXR( std::string const& filename,
                              ImageFormat const& format )
      : DiskImageResource( filename )
    {
      m_input_file_ptr = 0;
      m_output_file_ptr = 0;
      create( filename, format );
    }

    virtual ~DiskImageResourceOpenEXR();

    /// Returns the type of disk image resource.
    static std::string type_static() { return "OpenEXR"; }

    /// Returns the type of disk image resource.
    virtual std::string type() { return type_static(); }

    virtual void read( ImageBuffer const& dest, BBox2i const& bbox ) const;

    virtual void write( ImageBuffer const& dest, BBox2i const& bbox );

    void open( std::string const& filename );

    void create( std::string const& filename,
                 ImageFormat const& format );

    static DiskImageResource* construct_open( std::string const& filename );

    static DiskImageResource* construct_create( std::string const& filename,
                                                ImageFormat const& format );

    virtual bool has_block_write () const {return true; }
    virtual bool has_nodata_write() const {return false;}
    virtual bool has_block_read  () const {return true; }
    virtual bool has_nodata_read () const {return false;}

    virtual Vector2i block_read_size () const;
    virtual Vector2i block_write_size() const;
    virtual void set_block_write_size(const Vector2i&);

  private:
    const static int m_openexr_rows_per_block = 10;

    std::string m_filename;
    Vector2i    m_block_size;
    std::vector<std::string> m_labels;
    void* m_input_file_ptr;
    void* m_output_file_ptr;
    bool  m_tiled;
  };

} // namespace vw

#endif // __VW_FILEIO_DISKIMAGERESOUCEOPENEXR_H__
