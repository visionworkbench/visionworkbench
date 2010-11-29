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

#include <boost/assign/list_of.hpp>
#include <boost/lambda/construct.hpp>
#include <boost/function.hpp>
#include <boost/algorithm/string.hpp>
#include <map>


namespace {
  typedef boost::function<vw::SrcMemoryImageResource*(const vw::uint8*, size_t)> open_func;
  typedef boost::function<vw::DstMemoryImageResource*(std::vector<vw::uint8>*, const vw::ImageFormat&)> create_func;

  typedef std::map<std::string, open_func> open_map_t;
  typedef std::map<std::string, create_func> create_map_t;

#define ADD1(Name, Type) (Name, boost::lambda::new_ptr<vw::SrcMemoryImageResource ## Type>())
#define ADD2(Name, Type) (Name, boost::lambda::new_ptr<vw::DstMemoryImageResource ## Type>())

  open_map_t open_map = boost::assign::list_of<std::pair<std::string, open_func> >
#if defined(VW_HAVE_PKG_JPEG) && VW_HAVE_PKG_JPEG==1
    ADD1("jpg",        JPEG)
    ADD1("jpeg",       JPEG)
    ADD1("image/jpeg", JPEG)
#endif
#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1
    ADD1("png",        PNG)
    ADD1("image/png",  PNG)
#endif
    ;

  create_map_t create_map = boost::assign::list_of<std::pair<std::string, create_func> >
#if defined(VW_HAVE_PKG_JPEG) && VW_HAVE_PKG_JPEG==1
    ADD2("jpg",        JPEG)
    ADD2("jpeg",       JPEG)
    ADD2("image/jpeg", JPEG)
#endif
#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1
    ADD2("png",        PNG)
    ADD2("image/png",  PNG)
#endif
    ;

   std::string clean_type(const std::string& type) {
     return boost::to_lower_copy(boost::trim_left_copy_if(type, boost::is_any_of(".")));
   }
}

namespace vw {

  SrcMemoryImageResource* SrcMemoryImageResource::open( const std::string& type, const uint8* data, size_t len ) {
    open_map_t::const_iterator i = open_map.find(clean_type(type));
    if (i == open_map.end())
      vw_throw( NoImplErr() << "Unsupported file format: " << type );
    return i->second(data, len);
  }

  DstMemoryImageResource* DstMemoryImageResource::create( const std::string& type, std::vector<uint8>* data, const ImageFormat& format ) {
    create_map_t::const_iterator i = create_map.find(clean_type(type));
    if (i == create_map.end())
      vw_throw( NoImplErr() << "Unsupported file format: " << type );
    return i->second(data, format);
  }


} // namespace vw
