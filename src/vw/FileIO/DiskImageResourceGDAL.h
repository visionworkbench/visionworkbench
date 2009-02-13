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
/// Advanced users can pass custom options to GDAL when creating a 
/// resource.  Here is an example showing how to make a tiled, 
/// compressed BigTIFF (assuming libtiff 4.0 or greater):
///
///   DiskImageResourceGDAL::Options options;
///   options["COMPRESS"] = "LZW";
///   options["BIGTIFF"] = "YES";
///   DiskImageResourceGDAL resource( "filename.tif",
///                                   image.format(),
///                                   Vector2i(256,256),
///                                   options );
///   write_image( resource, image );
///
#ifndef __VW_FILEIO_DISKIMAGERESOUCEGDAL_H__
#define __VW_FILEIO_DISKIMAGERESOUCEGDAL_H__

#include <vw/config.h>
#include <string>
#include <map>

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

    typedef std::map<std::string,std::string> Options;

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
    
    DiskImageResourceGDAL( std::string const& filename, 
                           ImageFormat const& format,
                           Vector2i block_size,
                           Options const& options )
      : DiskImageResource( filename )
    {
      m_write_dataset_ptr = NULL;
      m_convert_jp2 = false;
      create( filename, format, block_size, options );
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

    virtual Vector2i block_size() const;
    virtual void set_block_size(Vector2i const&);

    virtual void flush();

    // Ask GDAL if it's compiled with support for this file
    static bool gdal_has_support(std::string const& filename);

    /// Some GDAL files store a "NoData" value.  You can use these
    /// methods to access that information. 
    double get_no_data_value(int band);

    void open( std::string const& filename );    
    void create( std::string const& filename,
                 ImageFormat const& format,
                 Vector2i block_size,
                 Options const& options );

    void create( std::string const& filename,
                 ImageFormat const& format,
                 Vector2i block_size = Vector2i(-1,-1) )
    {
      create( filename, format, block_size, Options() );
    }
    
    static DiskImageResource* construct_open( std::string const& filename );

    static DiskImageResource* construct_create( std::string const& filename,
                                                ImageFormat const& format );

    char **get_metadata() const;

    void* get_write_dataset_ptr() const { return m_write_dataset_ptr; }
    void* get_read_dataset_ptr() const { return (*(m_dataset_cache_handle)).get_dataset(); }

  private:
    static vw::Cache& gdal_cache();
    void initialize_write_resource();
    Vector2i default_block_size();

    std::string m_filename;
    void* m_write_dataset_ptr;
    Matrix<double,3,3> m_geo_transform;
    bool m_convert_jp2;
    std::vector<PixelRGBA<uint8> > m_palette;
    Vector2i m_blocksize;
    Options m_options;
    Cache::Handle<GdalDatasetGenerator> m_dataset_cache_handle;
  };

} // namespace vw

#endif // __VW_FILEIO_DISKIMAGERESOUCEGDAL_H__
