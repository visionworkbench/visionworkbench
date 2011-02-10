// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

// Forward declarations
class GDALDataset;
namespace vw {
  class Mutex;
}

namespace vw {

  class DiskImageResourceGDAL : public DiskImageResource {
    bool nodata_read_ok(double& value) const;
  public:

    typedef std::map<std::string,std::string> Options;

    DiskImageResourceGDAL( std::string const& filename )
      : DiskImageResource( filename )
    {
      open( filename );
    }

    DiskImageResourceGDAL( std::string const& filename,
                           ImageFormat const& format,
                           Vector2i block_size = Vector2i(-1,-1) )
      : DiskImageResource( filename )
    {
      create( filename, format, block_size );
    }

    DiskImageResourceGDAL( std::string const& filename,
                           ImageFormat const& format,
                           Vector2i block_size,
                           Options const& options )
      : DiskImageResource( filename )
    {
      create( filename, format, block_size, options );
    }

    virtual ~DiskImageResourceGDAL();

    /// Returns the type of disk image resource.
    static std::string type_static() { return "GDAL"; }
    static void set_gdal_cache_size(int size);  // Set GDAL cache size in bytes

    /// Returns the type of disk image resource.
    virtual std::string type() { return type_static(); }

    virtual void read( ImageBuffer const& dest, BBox2i const& bbox ) const;
    virtual void write( ImageBuffer const& dest, BBox2i const& bbox );

    virtual bool has_block_read()   const {return true;}
    virtual bool has_block_write()  const {return true;}
    virtual bool has_nodata_read()  const;
    virtual bool has_nodata_write() const {return true;}

    virtual Vector2i block_write_size() const;
    virtual void set_block_write_size(const Vector2i&);
    virtual Vector2i block_read_size() const;

    virtual void set_nodata_write(double);
    virtual double nodata_read() const;

    virtual void flush();

    // Ask GDAL if it's compiled with support for this file
    static bool gdal_has_support(std::string const& filename);

    void open( std::string const& filename );
    void create( std::string const& filename,
                 ImageFormat const& format,
                 Vector2i block_size,
                 Options const& options );

    void create( std::string const& filename,
                 ImageFormat const& format,
                 Vector2i block_size = Vector2i(-1,-1) )
    {
      std::string extension = boost::to_lower_copy(boost::filesystem::extension(filename));
      if ( extension == ".tif" || extension == ".tiff" ) {
        // TIFF should use LZW by default
        Options tiff_options;
        tiff_options["COMPRESS"] = "LZW";
        create( filename, format, block_size, tiff_options );
      } else {
        create( filename, format, block_size, Options() );
      }
    }

    static DiskImageResource* construct_open( std::string const& filename );

    static DiskImageResource* construct_create( std::string const& filename,
                                                ImageFormat const& format );

    // These functions return pointers to internal data.  They exist
    // to allow users to access underlying special-purpose GDAL
    // features, but they should be used with caution.  Unlike the
    // rest of the public interface, they are not thread-safe and if
    // you use them you must be sure to acquire the global GDAL lock
    // (accessed via the global_lock() function) for the duration of
    // your use, up to and including the release of your shared
    // pointer to the dataset.
    boost::shared_ptr<GDALDataset> get_dataset_ptr() const;
    char **get_metadata() const;

    // Provides access to the global GDAL lock, for users who need to go
    // tinkering around in GDAL directly for some reason.
    static Mutex &global_lock();

  private:
    void initialize_write_resource_locked();
    Vector2i default_block_size();

    std::string m_filename;
    boost::shared_ptr<GDALDataset> m_write_dataset_ptr;
    std::vector<PixelRGBA<uint8> > m_palette;
    Vector2i m_blocksize;
    Options m_options;
    boost::shared_ptr<GDALDataset> m_read_dataset_ptr;
  };

  void UnloadGDAL();

} // namespace vw

#endif // __VW_FILEIO_DISKIMAGERESOUCEGDAL_H__
