// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
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

// For RunOnce
#include <vw/Core/Thread.h>

#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageResourcePDS.h>

#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1
#include <vw/FileIO/DiskImageResourcePNG.h>
#endif

#if defined(VW_HAVE_PKG_JPEG) && VW_HAVE_PKG_JPEG==1
#include <vw/FileIO/DiskImageResourceJPEG.h>
#endif

#if defined(VW_HAVE_PKG_JPEG2K) && VW_HAVE_PKG_JPEG2K==1
#include <vw/FileIO/DiskImageResourceJP2.h>
#endif

#if defined(VW_HAVE_PKG_TIFF) && VW_HAVE_PKG_TIFF==1
#include <vw/FileIO/DiskImageResourceTIFF.h>
#endif

#if defined(VW_HAVE_PKG_OPENEXR) && VW_HAVE_PKG_OPENEXR==1
#if ! ( defined(VW_NO_EXCEPTIONS) && VW_NO_EXCEPTIONS==1 )
#include <vw/FileIO/DiskImageResourceOpenEXR.h>
#endif
#endif

#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
#include <vw/FileIO/DiskImageResourceGDAL.h>
#endif

#include <vw/FileIO/DiskImageResource_internal.h>

namespace {
  typedef std::map<std::string,vw::DiskImageResource::construct_open_func> OpenMapType;
  typedef std::map<std::string,vw::DiskImageResource::construct_create_func> CreateMapType;
  OpenMapType *open_map = 0;
  CreateMapType *create_map = 0;
}

namespace vw {
  namespace internal {

void foreach_ext(std::string fn, ExtTestFunction func, std::set<std::string> exclude)
{
  OpenMapType::const_iterator oi;
  vw::DiskImageResource::register_default_file_types();

  for (oi = open_map->begin(); oi != open_map->end(); ++oi)
  {
    if (exclude.find(oi->first.substr(1)) == exclude.end())
      func(fn + oi->first);
  }
}

}}

void vw::DiskImageResource::register_file_type( std::string const& extension,
                                                std::string const& disk_image_resource_type,
                                                vw::DiskImageResource::construct_open_func open_func,
                                                vw::DiskImageResource::construct_create_func create_func )
{
  if( ! open_map ) open_map = new OpenMapType();
  if( ! create_map ) create_map = new CreateMapType();
  //  std::cout << "REGISTERING DiskImageResource " << disk_image_resource_type << " for extension " << extension << std::endl;

  OpenMapType::iterator om_iter = open_map->find( extension );
  if ( om_iter != open_map->end() )
    om_iter->second = open_func;
  else
    open_map->insert( std::make_pair( extension, open_func ) );

  CreateMapType::iterator cm_iter = create_map->find( extension );
  if ( cm_iter != create_map->end() )
    cm_iter->second = create_func;
  else
    create_map->insert( std::make_pair( extension, create_func ) );
    
}

static std::string file_extension( std::string const& filename ) {
  std::string::size_type dot = filename.find_last_of('.');
  if (dot == std::string::npos)
    vw_throw( vw::IOErr() << "DiskImageResource: Cannot infer file format from filename with no file extension." );
  std::string extension = filename.substr( dot );
  boost::to_lower( extension );
  return extension;
}

static vw::RunOnce rdft_once = VW_RUNONCE_INIT;

static void register_default_file_types_impl() {
  vw::DiskImageResource::register_file_type( ".img", vw::DiskImageResourcePDS::type_static(), &vw::DiskImageResourcePDS::construct_open, &vw::DiskImageResourcePDS::construct_create );
  vw::DiskImageResource::register_file_type( ".pds", vw::DiskImageResourcePDS::type_static(), &vw::DiskImageResourcePDS::construct_open, &vw::DiskImageResourcePDS::construct_create );
  vw::DiskImageResource::register_file_type( ".lbl", vw::DiskImageResourcePDS::type_static(), &vw::DiskImageResourcePDS::construct_open, &vw::DiskImageResourcePDS::construct_create );

#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1
  vw::DiskImageResource::register_file_type( ".png", vw::DiskImageResourcePNG::type_static(), &vw::DiskImageResourcePNG::construct_open, &vw::DiskImageResourcePNG::construct_create );
#elif defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
  if (vw::DiskImageResourceGDAL::gdal_has_support(".png"))
    vw::DiskImageResource::register_file_type( ".png", vw::DiskImageResourceGDAL::type_static(), &vw::DiskImageResourceGDAL::construct_open, &vw::DiskImageResourceGDAL::construct_create );
#endif

#if defined(VW_HAVE_PKG_JPEG) && VW_HAVE_PKG_JPEG==1
  vw::DiskImageResource::register_file_type( ".jpg", vw::DiskImageResourceJPEG::type_static(), &vw::DiskImageResourceJPEG::construct_open, &vw::DiskImageResourceJPEG::construct_create );
  vw::DiskImageResource::register_file_type( ".jpeg", vw::DiskImageResourceJPEG::type_static(), &vw::DiskImageResourceJPEG::construct_open, &vw::DiskImageResourceJPEG::construct_create );
#elif
  if (vw::DiskImageResourceGDAL::gdal_has_support(".jpg"))
    vw::DiskImageResource::register_file_type( ".jpg", vw::DiskImageResourceGDAL::type_static(), &vw::DiskImageResourceGDAL::construct_open, &vw::DiskImageResourceGDAL::construct_create );
  if (vw::DiskImageResourceGDAL::gdal_has_support(".jpeg"))
    vw::DiskImageResource::register_file_type( ".jpeg", vw::DiskImageResourceGDAL::type_static(), &vw::DiskImageResourceGDAL::construct_open, &vw::DiskImageResourceGDAL::construct_create );
#endif

#if defined(VW_HAVE_PKG_JPEG2K) && VW_HAVE_PKG_JPEG2K==1 && 0
  // A file with a .jp2 extension is a full fledged JPEG2000 image
  // with acquisition metadata. A file with a .j2k extension has only
  // the "raw" encoded image, with image encoding and size specified
  // in a small header. A file with a .jpf extension is a full fledged
  // JPEG2000 image with acquisition and (possibly) GML metadata.
  vw::DiskImageResource::register_file_type(".jp2", vw::DiskImageResourceJP2::type_static(), &vw::DiskImageResourceJP2::construct_open, &vw::DiskImageResourceJP2::construct_create );
  vw::DiskImageResource::register_file_type(".j2k", vw::DiskImageResourceJP2::type_static(), &vw::DiskImageResourceJP2::construct_open, &vw::DiskImageResourceJP2::construct_create );
  vw::DiskImageResource::register_file_type(".jpf", vw::DiskImageResourceJP2::type_static(), &vw::DiskImageResourceJP2::construct_open, &vw::DiskImageResourceJP2::construct_create );
#elif defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
  if (vw::DiskImageResourceGDAL::gdal_has_support(".jp2"))
    vw::DiskImageResource::register_file_type(".jp2", vw::DiskImageResourceGDAL::type_static(), &vw::DiskImageResourceGDAL::construct_open, &vw::DiskImageResourceGDAL::construct_create );
  if (vw::DiskImageResourceGDAL::gdal_has_support(".j2k"))
    vw::DiskImageResource::register_file_type(".j2k", vw::DiskImageResourceGDAL::type_static(), &vw::DiskImageResourceGDAL::construct_open, &vw::DiskImageResourceGDAL::construct_create );
#endif

// This is a little hackish but it makes it so libtiff acts as a proper fallback
#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
  if (vw::DiskImageResourceGDAL::gdal_has_support(".tif") && vw::DiskImageResourceGDAL::gdal_has_support(".tiff")) {
    vw::DiskImageResource::register_file_type( ".tif", vw::DiskImageResourceGDAL::type_static(), &vw::DiskImageResourceGDAL::construct_open, &vw::DiskImageResourceGDAL::construct_create );
    vw::DiskImageResource::register_file_type( ".tiff", vw::DiskImageResourceGDAL::type_static(), &vw::DiskImageResourceGDAL::construct_open, &vw::DiskImageResourceGDAL::construct_create );
  } else {
#endif
#if defined(VW_HAVE_PKG_TIFF) && VW_HAVE_PKG_TIFF==1
    vw::DiskImageResource::register_file_type( ".tif", vw::DiskImageResourceTIFF::type_static(), &vw::DiskImageResourceTIFF::construct_open, &vw::DiskImageResourceTIFF::construct_create );
    vw::DiskImageResource::register_file_type( ".tiff", vw::DiskImageResourceTIFF::type_static(), &vw::DiskImageResourceTIFF::construct_open, &vw::DiskImageResourceTIFF::construct_create );
#endif
#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
  }
#endif

#if defined(VW_HAVE_PKG_OPENEXR) && VW_HAVE_PKG_OPENEXR==1
  vw::DiskImageResource::register_file_type( ".exr", vw::DiskImageResourceOpenEXR::type_static(), &vw::DiskImageResourceOpenEXR::construct_open, &vw::DiskImageResourceOpenEXR::construct_create );
#endif

#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
  if (vw::DiskImageResourceGDAL::gdal_has_support(".grd"))
    vw::DiskImageResource::register_file_type( ".grd", vw::DiskImageResourceGDAL::type_static(), &vw::DiskImageResourceGDAL::construct_open, &vw::DiskImageResourceGDAL::construct_create );
  if (vw::DiskImageResourceGDAL::gdal_has_support(".dem"))
    vw::DiskImageResource::register_file_type( ".dem", vw::DiskImageResourceGDAL::type_static(), &vw::DiskImageResourceGDAL::construct_open, &vw::DiskImageResourceGDAL::construct_create );
  if (vw::DiskImageResourceGDAL::gdal_has_support(".bil"))
    vw::DiskImageResource::register_file_type( ".bil", vw::DiskImageResourceGDAL::type_static(), &vw::DiskImageResourceGDAL::construct_open, &vw::DiskImageResourceGDAL::construct_create );
  if (vw::DiskImageResourceGDAL::gdal_has_support(".cub"))
    vw::DiskImageResource::register_file_type( ".cub", vw::DiskImageResourceGDAL::type_static(), &vw::DiskImageResourceGDAL::construct_open, &vw::DiskImageResourceGDAL::construct_create );
  if (vw::DiskImageResourceGDAL::gdal_has_support(".gif"))
    vw::DiskImageResource::register_file_type( ".gif", vw::DiskImageResourceGDAL::type_static(), &vw::DiskImageResourceGDAL::construct_open, &vw::DiskImageResourceGDAL::construct_create );
#endif

}

// This extra class helps to ensure that register_file_type() is only run once.
void vw::DiskImageResource::register_default_file_types() {
  rdft_once.run( register_default_file_types_impl );
}

vw::DiskImageResource* vw::DiskImageResource::open( std::string const& filename ) {
  register_default_file_types();
  if( open_map ) {
    OpenMapType::iterator i = open_map->find( file_extension( filename ) );
    if( i != open_map->end() ) 
      return i->second( filename );
  }
  vw_throw( NoImplErr() << "Unsuppported file format: " << filename );
  return 0; // never reached
}

vw::DiskImageResource* vw::DiskImageResource::create( std::string const& filename, ImageFormat const& format ) {
  register_default_file_types();
  if( create_map ) {
    CreateMapType::iterator i = create_map->find( file_extension( filename ) );
    if( i != create_map->end() )
      return i->second( filename, format );
  }
  vw_throw( NoImplErr() << "Unsuppported file format: " << filename );
  return 0; // never reached
}
