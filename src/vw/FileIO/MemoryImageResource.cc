// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/config.h>
#include <vw/FileIO/MemoryImageResource.h>

#if defined(VW_HAVE_PKG_JPEG) && VW_HAVE_PKG_JPEG==1
#  include <vw/FileIO/MemoryImageResourceJPEG.h>
#endif

#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1
#  include <vw/FileIO/MemoryImageResourcePNG.h>
#endif

#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
#  include <vw/FileIO/MemoryImageResourceGDAL.h>
#endif

#include <boost/assign/list_of.hpp>
#include <boost/lambda/construct.hpp>
#include <boost/function.hpp>
#include <boost/algorithm/string.hpp>
#include <map>


namespace {
  typedef boost::function<vw::SrcMemoryImageResource*(const boost::shared_array<const vw::uint8>, size_t)> open_func;
  typedef boost::function<vw::DstMemoryImageResource*(const vw::ImageFormat&)> create_func;

  typedef std::map<std::string, open_func> open_map_t;
  typedef std::map<std::string, create_func> create_map_t;

#define OPEN(Name, Type) (Name, boost::lambda::new_ptr<vw::SrcMemoryImageResource ## Type>())
#define CREAT(Name, Type) (Name, boost::lambda::new_ptr<vw::DstMemoryImageResource ## Type>())

  open_map_t open_map = boost::assign::list_of<std::pair<std::string, open_func> >
#if defined(VW_HAVE_PKG_JPEG) && VW_HAVE_PKG_JPEG==1
    OPEN("jpg",        JPEG)
    OPEN("jpeg",       JPEG)
    OPEN("image/jpeg", JPEG)
#endif
#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1
    OPEN("png",        PNG)
    OPEN("image/png",  PNG)
#endif
#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
    OPEN("tif",        GDAL)
    OPEN("tiff",       GDAL)
    OPEN("image/tiff", GDAL)
#endif
    ;

  create_map_t create_map = boost::assign::list_of<std::pair<std::string, create_func> >
#if defined(VW_HAVE_PKG_JPEG) && VW_HAVE_PKG_JPEG==1
    CREAT("jpg",        JPEG)
    CREAT("jpeg",       JPEG)
    CREAT("image/jpeg", JPEG)
#endif
#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1
    CREAT("png",        PNG)
    CREAT("image/png",  PNG)
#endif
#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
    CREAT("tif",        GDAL)
    CREAT("tiff",       GDAL)
    CREAT("image/tiff", GDAL)
#endif
    ;

   std::string clean_type(const std::string& type) {
     return boost::to_lower_copy(boost::trim_left_copy_if(type, boost::is_any_of(".")));
   }
}

namespace vw {

  SrcMemoryImageResource* SrcMemoryImageResource::open( const std::string& type, const uint8* data, size_t len ) {
    boost::shared_array<const uint8> p(data, NOP());
    return SrcMemoryImageResource::open(type, p, len);
  }

  SrcMemoryImageResource* SrcMemoryImageResource::open( const std::string& type, boost::shared_array<const uint8> data, size_t len ) {
    open_map_t::const_iterator i = open_map.find(clean_type(type));
    if (i == open_map.end())
      vw_throw( NoImplErr() << "Unsupported file format: " << type );
    return i->second(data, len);
  }

  DstMemoryImageResource* DstMemoryImageResource::create( const std::string& type, const ImageFormat& format ) {
    create_map_t::const_iterator i = create_map.find(clean_type(type));
    if (i == create_map.end())
      vw_throw( NoImplErr() << "Unsupported file format: " << type );
    return i->second(format);
  }

} // namespace vw
