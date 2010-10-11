// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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

  class DiskImageResourceOpenEXR : public DiskImageResource {
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

    virtual Vector2i block_size() const;
    virtual void set_block_size(Vector2i const&);

    void set_tiled_write(int32 tile_width, int32 tile_height, bool random_tile_order = false);
    void set_scanline_write(int32 scanlines_per_block);

    virtual void read( ImageBuffer const& dest, BBox2i const& bbox ) const;

    virtual void write( ImageBuffer const& dest, BBox2i const& bbox );

    void open( std::string const& filename );

    void create( std::string const& filename,
                 ImageFormat const& format );

    static DiskImageResource* construct_open( std::string const& filename );

    static DiskImageResource* construct_create( std::string const& filename,
                                                ImageFormat const& format );

  private:
    const static int m_openexr_rows_per_block = 10;

    std::string m_filename;
    Vector2i m_block_size;
    std::vector<std::string> m_labels;
    void* m_input_file_ptr;
    void* m_output_file_ptr;
    bool m_tiled;
  };

} // namespace vw

#endif // __VW_FILEIO_DISKIMAGERESOUCEOPENEXR_H__
