// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file DiskImageResource.cc
///
/// An abstract base class referring to an image on disk.
///

#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <vw/config.h>

#include <iostream>
#include <set>
#include <map>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

// For RunOnce
#include <vw/Core/Thread.h>

#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageResourcePDS.h>
#include <vw/FileIO/DiskImageResourcePBM.h>

#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1
#include <vw/FileIO/DiskImageResourcePNG.h>
#endif

#if defined(VW_HAVE_PKG_JPEG) && VW_HAVE_PKG_JPEG==1
#include <vw/FileIO/DiskImageResourceJPEG.h>
#endif

#if defined(VW_HAVE_PKG_TIFF) && VW_HAVE_PKG_TIFF==1
#include <vw/FileIO/DiskImageResourceTIFF.h>
#endif

#if defined(VW_HAVE_PKG_OPENEXR) && VW_HAVE_PKG_OPENEXR==1
#if defined(VW_ENABLE_EXCEPTIONS) && VW_ENABLE_EXCEPTIONS==1
#include <vw/FileIO/DiskImageResourceOpenEXR.h>
#endif
#endif

#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
#include <vw/FileIO/DiskImageResourceGDAL.h>
#endif

#include <vw/FileIO/DiskImageResource_internal.h>

static void register_default_file_types_impl();
namespace {
  typedef std::map<std::string,vw::DiskImageResource::construct_open_func> OpenMapType;
  typedef std::map<std::string,vw::DiskImageResource::construct_create_func> CreateMapType;
  OpenMapType *open_map = 0;
  CreateMapType *create_map = 0;

  // This extra class helps to ensure that register_file_type() is only run once.
  vw::RunOnce rdft_once = VW_RUNONCE_INIT;
  void register_default_file_types_internal() {
    rdft_once.run( register_default_file_types_impl );
  }

  // this one avoids calling the registration function, so it can be called
  // from INSIDE the registration function.
  void register_file_type_internal( std::string const& extension,
      std::string const& /*disk_image_resource_type*/,
      vw::DiskImageResource::construct_open_func open_func,
      vw::DiskImageResource::construct_create_func create_func ) {

    // This will create the entries if they don't exist
    (*open_map)[extension]   = open_func;
    (*create_map)[extension] = create_func;
  }
}

bool vw::DiskImageResource::default_rescale = true;

void vw::DiskImageResource::set_rescale(bool rescale) {
  m_rescale = rescale;
}

void vw::DiskImageResource::set_default_rescale(bool rescale) {
  default_rescale = rescale;
}

namespace vw {
  namespace internal {

void foreach_ext(std::string const& prefix, ExtTestFunction const& callback, std::set<std::string> const& exclude)
{
  OpenMapType::const_iterator oi;
  register_default_file_types_internal();

  for (oi = open_map->begin(); oi != open_map->end(); ++oi)
  {
    if (exclude.find(oi->first.substr(1)) == exclude.end())
      callback(prefix + oi->first);
  }
}

}}

void vw::DiskImageResource::register_file_type( std::string const& extension,
                                                std::string const& disk_image_resource_type,
                                                vw::DiskImageResource::construct_open_func open_func,
                                                vw::DiskImageResource::construct_create_func create_func )
{
  register_default_file_types_internal();

  // Add the file to the list
  register_file_type_internal(boost::to_lower_copy(extension), disk_image_resource_type, open_func, create_func);
}

static void register_default_file_types_impl() {

  if( ! open_map ) open_map = new OpenMapType();
  if( ! create_map ) create_map = new CreateMapType();

// Let's cut the verbosity of this func just a bit.
#define REGISTER(ext, driver) register_file_type_internal( ext, vw::DiskImageResource ## driver::type_static(), &vw::DiskImageResource ## driver::construct_open, &vw::DiskImageResource ## driver::construct_create );

  // Give GDAL precedence in reading PDS images when this is supported.
#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
  if (vw::DiskImageResourceGDAL::gdal_has_support(".img") &&
      vw::DiskImageResourceGDAL::gdal_has_support(".pds") &&
      vw::DiskImageResourceGDAL::gdal_has_support(".lbl")) {

    REGISTER(".img", GDAL)
    REGISTER(".pds", GDAL)
    REGISTER(".lbl", GDAL)
  } else {
#endif
  REGISTER(".img", PDS)
  REGISTER(".pds", PDS)
  REGISTER(".lbl", PDS)
#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
  }
#endif

#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1
  REGISTER(".png", PNG)
#elif defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
  if (vw::DiskImageResourceGDAL::gdal_has_support(".png"))
    REGISTER(".png", GDAL)
  else
    vw::vw_throw(vw::IOErr() << "GDAL does not have PNG support.");
#endif

#if defined(VW_HAVE_PKG_JPEG) && VW_HAVE_PKG_JPEG==1
  REGISTER(".jpg", JPEG)
  REGISTER(".jpeg", JPEG)
#elif defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
  if (vw::DiskImageResourceGDAL::gdal_has_support(".jpg"))
    REGISTER(".jpg", GDAL)
  if (vw::DiskImageResourceGDAL::gdal_has_support(".jpeg"))
    REGISTER(".jpeg", GDAL)
#endif

#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
  if (vw::DiskImageResourceGDAL::gdal_has_support(".jp2"))
    REGISTER(".jp2", GDAL)
  if (vw::DiskImageResourceGDAL::gdal_has_support(".j2k"))
    REGISTER(".j2k", GDAL)
#endif

// This is a little hackish but it makes it so libtiff acts as a proper fallback
#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
  if (vw::DiskImageResourceGDAL::gdal_has_support(".tif") && vw::DiskImageResourceGDAL::gdal_has_support(".tiff")) {
    REGISTER(".tif", GDAL)
    REGISTER(".tiff", GDAL)
  } else {
#endif
#if defined(VW_HAVE_PKG_TIFF) && VW_HAVE_PKG_TIFF==1
    REGISTER(".tif", TIFF)
    REGISTER(".tiff", TIFF)
#endif
#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
  }
#endif

#if defined(VW_HAVE_PKG_OPENEXR) && VW_HAVE_PKG_OPENEXR==1
  REGISTER(".exr", OpenEXR)
#endif

  // Filetypes that are always supported
  REGISTER(".pbm", PBM)
  REGISTER(".pgm", PBM)
  REGISTER(".ppm", PBM)
#undef REGISTER
}

// Kill this function eventually.. it's marked as deprecated now.
void vw::DiskImageResource::register_default_file_types() {
  register_default_file_types_internal();
}

vw::DiskImageResource* vw::DiskImageResource::open( std::string const& filename ) {
  register_default_file_types_internal();
  std::string extension = boost::to_lower_copy(fs::extension(filename));

  if( open_map ) {
    OpenMapType::iterator i = open_map->find( extension );
    if( i != open_map->end() ) {
      DiskImageResource* rsrc = i->second( filename );
      vw_out(DebugMessage,"fileio") << "Produce DiskImageResource of type: " << rsrc->type() << "\n";
      return rsrc;
    }
  }

  // GDAL has support for many useful file formats, and we fall back
  // on it here in case none of the registered file handlers know how
  // to do the job.
#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
  if (vw::DiskImageResourceGDAL::gdal_has_support( extension ))
    return vw::DiskImageResourceGDAL::construct_open(filename);
#endif

  // If all attempts to find a suitable file driver fails, we throw an
  // exception.
  vw_throw( NoImplErr() << "Unsupported file format: " << filename );
  return 0; // never reached
}

/// Returns a disk image resource with the given filename.  The file
/// type is determined by the value in 'type'.
vw::DiskImageResource* vw::DiskImageResource::create( std::string const& filename, ImageFormat const& format, std::string const& type ) {
  register_default_file_types_internal();
  if( create_map ) {
    CreateMapType::iterator i = create_map->find( boost::to_lower_copy(type) );
    if( i != create_map->end() )
      return i->second( filename, format );
  }
  vw_throw( NoImplErr() << "Unsupported file type \"" << type << "\" for filename: " << filename );
  return 0; // never reached
}

/// Returns a disk image resource with the given filename.  The file
/// type is determined by the extension of the filename.
vw::DiskImageResource* vw::DiskImageResource::create( std::string const& filename, ImageFormat const& format ) {
  register_default_file_types_internal();
  if( create_map ) {
    CreateMapType::iterator i = create_map->find( boost::to_lower_copy(fs::extension( filename )) );
    if( i != create_map->end() )
      return i->second( filename, format );
  }
  vw_throw( NoImplErr() << "Unsupported file format: " << filename );
  return 0; // never reached
}
