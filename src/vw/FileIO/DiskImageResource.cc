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
#include <map>
#include <boost/algorithm/string.hpp>

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


namespace {
  typedef std::map<std::string,vw::DiskImageResource::construct_open_func> OpenMapType;
  typedef std::map<std::string,vw::DiskImageResource::construct_create_func> CreateMapType;
  OpenMapType *open_map = 0;
  CreateMapType *create_map = 0;
}

void vw::DiskImageResource::register_file_type( std::string const& extension,
                                                vw::DiskImageResource::construct_open_func open_func,
                                                vw::DiskImageResource::construct_create_func create_func )
{
  if( ! open_map ) open_map = new OpenMapType();
  if( ! create_map ) create_map = new CreateMapType();
  open_map->insert( std::make_pair( extension, open_func ) );
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

static void register_default_file_types() {
  static bool already = false;
  if( already ) return;
  already = true;
  vw::DiskImageResource::register_file_type( ".img", &vw::DiskImageResourcePDS::construct_open, &vw::DiskImageResourcePDS::construct_create );

#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1
  vw::DiskImageResource::register_file_type( ".png", &vw::DiskImageResourcePNG::construct_open, &vw::DiskImageResourcePNG::construct_create );
#endif
#if defined(VW_HAVE_PKG_JPEG) && VW_HAVE_PKG_JPEG==1
  vw::DiskImageResource::register_file_type( ".jpg", &vw::DiskImageResourceJPEG::construct_open, &vw::DiskImageResourceJPEG::construct_create );
  vw::DiskImageResource::register_file_type( ".jpeg", &vw::DiskImageResourceJPEG::construct_open, &vw::DiskImageResourceJPEG::construct_create );
#endif

#if defined(VW_HAVE_PKG_JPEG2K) && VW_HAVE_PKG_JPEG2K==1
  // A file with a .jp2 extension is a full fledged JPEG2000 image
  // with acquisition metadata. A file with a .j2k extension has only
  // the "raw" encoded image, with image encoding and size specified
  // in a small header. A file with a .jpf extension is a full fledged
  // JPEG2000 image with acquisition and (possibly) GML metadata.
  std::cout << "------------- REGISTERING JP2 -----------------" << std::endl;
  vw::DiskImageResource::register_file_type(".jp2", &vw::DiskImageResourceJP2::construct_open, &vw::DiskImageResourceJP2::construct_create);

  vw::DiskImageResource::register_file_type(".j2k", &vw::DiskImageResourceJP2::construct_open, &vw::DiskImageResourceJP2::construct_create);
  
  vw::DiskImageResource::register_file_type(".jpf", &vw::DiskImageResourceJP2::construct_open, &vw::DiskImageResourceJP2::construct_create);

#endif

#if defined(VW_HAVE_PKG_TIFF) && VW_HAVE_PKG_TIFF==1
  vw::DiskImageResource::register_file_type( ".tif", &vw::DiskImageResourceTIFF::construct_open, &vw::DiskImageResourceTIFF::construct_create );
  vw::DiskImageResource::register_file_type( ".tiff", &vw::DiskImageResourceTIFF::construct_open, &vw::DiskImageResourceTIFF::construct_create );
#endif
#if defined(VW_HAVE_PKG_OPENEXR) && VW_HAVE_PKG_OPENEXR==1
  vw::DiskImageResource::register_file_type( ".exr", &vw::DiskImageResourceOpenEXR::construct_open, &vw::DiskImageResourceOpenEXR::construct_create );
#endif
}

vw::DiskImageResource* vw::DiskImageResource::open( std::string const& filename ) {
  register_default_file_types();
  if( open_map ) {
    OpenMapType::const_iterator i = open_map->find( file_extension( filename ) );
    if( i != open_map->end() ) return i->second( filename );
  }
  vw_throw( NoImplErr() << "Unsuppported file format: " << filename );
  return 0; // never reached
}

vw::DiskImageResource* vw::DiskImageResource::create( std::string const& filename, ImageFormat const& format ) {
  register_default_file_types();
  if( create_map ) {
    CreateMapType::const_iterator i = create_map->find( file_extension( filename ) );
    if( i != create_map->end() ) return i->second( filename, format );
  }
  vw_throw( NoImplErr() << "Unsuppported file format: " << filename );
  return 0; // never reached
}
