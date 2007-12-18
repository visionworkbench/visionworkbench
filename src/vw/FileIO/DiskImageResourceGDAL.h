// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
//
// __END_LICENSE__

/// \file DiskImageResourceGDAL.h
/// 
/// Provides support for georeferenced files via the GDAL library.
///
#ifndef __VW_FILEIO_DISKIMAGERESOUCEGDAL_H__
#define __VW_FILEIO_DISKIMAGERESOUCEGDAL_H__

#include <vw/config.h>
#include <string>

// VW Headers
#include <vw/Image/PixelTypes.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Math/Matrix.h>
#include <vw/Core/Cache.h>

namespace vw {

  class GdalDatasetHandle {
    std::string m_filename;
    void* m_dataset_ptr;
    
  public:
    GdalDatasetHandle(std::string filename);
    void* get_dataset() const { return m_dataset_ptr; }
    ~GdalDatasetHandle();
  };

  // GdalDatasetGenerator
  class GdalDatasetGenerator {
    std::string m_filename;
  public:
    typedef GdalDatasetHandle value_type;

    GdalDatasetGenerator( std::string filename ) : m_filename( filename ) {}
    
    size_t size() const {
      return 1;
    }
    
    boost::shared_ptr<GdalDatasetHandle> generate() const {
      return boost::shared_ptr<GdalDatasetHandle> ( new GdalDatasetHandle(m_filename) );
    }
  };

  class DiskImageResourceGDAL : public DiskImageResource {
  public:

    DiskImageResourceGDAL( std::string const& filename )
      : DiskImageResource( filename )
    {
      m_write_dataset_ptr = NULL;
      m_convert_jp2 = false;
      open( filename );
    }

    DiskImageResourceGDAL( std::string const& filename, 
                           ImageFormat const& format,
                           Vector2i block_size = Vector2i(-1,-1) )
      : DiskImageResource( filename )
    {
      m_write_dataset_ptr = NULL;
      m_convert_jp2 = false;
      create( filename, format, block_size );
    }
    
    virtual ~DiskImageResourceGDAL() {
      flush();
    }

    /// Returns the type of disk image resource.
    static std::string type_static() { return "GDAL"; }
    
    /// Returns the type of disk image resource.
    virtual std::string type() { return type_static(); }
    
    virtual void read( ImageBuffer const& dest, BBox2i const& bbox ) const;
    virtual void write( ImageBuffer const& dest, BBox2i const& bbox );

    /// Set the native block size
    ///
    /// Be careful here -- you can set any block size here, but you
    /// choice may lead to extremely inefficient FileIO operations.  You
    /// can choose to pass in -1 as the width and/or the height of the
    /// block, in which case the width and/or height is chosen by GDAL.
    /// 
    /// For example, if you pass in a vector of (-1,-1), the block size will be
    /// assigned based on GDAL's best guess of the best block or strip
    /// size. However, GDAL assumes a single-row stripsize even for file
    /// formats like PNG for which it does not support true strip access.
    /// Thus, we check the file driver type before accepting GDAL's block
    /// size as our own.
    void set_native_block_size(Vector2i block_size);

    virtual Vector2i native_block_size() const;
    virtual void flush();


    void open( std::string const& filename );    
    void create( std::string const& filename,
                 ImageFormat const& format,
                 Vector2i m_block_size = Vector2i(-1,-1) );
    
    static DiskImageResource* construct_open( std::string const& filename );

    static DiskImageResource* construct_create( std::string const& filename,
                                                ImageFormat const& format );

    void* get_write_dataset_ptr() const { return m_write_dataset_ptr; }
    void* get_read_dataset_ptr() const { return (*(m_dataset_cache_handle)).get_dataset(); }

  private:
    static vw::Cache& gdal_cache();

    std::string m_filename;
    void* m_write_dataset_ptr;
    Matrix<double,3,3> m_geo_transform;
    bool m_convert_jp2;
    std::vector<PixelRGBA<uint8> > m_palette;
    Vector2i m_native_blocksize;
    Cache::Handle<GdalDatasetGenerator> m_dataset_cache_handle;
  };

} // namespace vw

#endif // __VW_FILEIO_DISKIMAGERESOUCEGDAL_H__
