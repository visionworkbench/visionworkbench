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


#include <vw/config.h>
#include <vw/FileIO/MemoryImageResource.h>

#if defined(VW_HAVE_PKG_JPEG)
#  include <vw/FileIO/MemoryImageResourceJPEG.h>
#endif

#if defined(VW_HAVE_PKG_PNG)
#  include <vw/FileIO/MemoryImageResourcePNG.h>
#endif

#if defined(VW_HAVE_PKG_GDAL)
#  include <vw/FileIO/MemoryImageResourceGDAL.h>
#endif

#if defined(VW_HAVE_PKG_OPENEXR)
#  include <vw/FileIO/MemoryImageResourceOpenEXR.h>
#endif

#include <map>

#include <boost/assign/list_of.hpp>
#include <boost/lambda/construct.hpp>
#include <boost/function.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/trim.hpp>

namespace {
  typedef boost::function<vw::SrcMemoryImageResource*(const boost::shared_array<const vw::uint8>, size_t)> open_func;
  typedef boost::function<vw::DstMemoryImageResource*(const vw::ImageFormat&)> create_func;

  typedef std::map<std::string, open_func> open_map_t;
  typedef std::map<std::string, create_func> create_map_t;

#define OPEN(Name, Type) (Name, boost::lambda::new_ptr<vw::SrcMemoryImageResource ## Type>())
#define CREAT(Name, Type) (Name, boost::lambda::new_ptr<vw::DstMemoryImageResource ## Type>())

  open_map_t open_map = boost::assign::list_of<std::pair<std::string, open_func> >
#if defined(VW_HAVE_PKG_JPEG)
    OPEN("jpg",        JPEG)
    OPEN("jpeg",       JPEG)
    OPEN("image/jpeg", JPEG)
#endif
#if defined(VW_HAVE_PKG_PNG)
    OPEN("png",        PNG)
    OPEN("image/png",  PNG)
#endif
#if defined(VW_HAVE_PKG_GDAL)
    OPEN("tif",        GDAL)
    OPEN("tiff",       GDAL)
    OPEN("image/tiff", GDAL)
#endif
#if defined(VW_HAVE_PKG_OPENEXR)
    OPEN("exr",        OpenEXR)
    OPEN("image/exr",  OpenEXR)
#endif
    ;

  create_map_t create_map = boost::assign::list_of<std::pair<std::string, create_func> >
#if defined(VW_HAVE_PKG_JPEG)
    CREAT("jpg",        JPEG)
    CREAT("jpeg",       JPEG)
    CREAT("image/jpeg", JPEG)
#endif
#if defined(VW_HAVE_PKG_PNG)
    CREAT("png",        PNG)
    CREAT("image/png",  PNG)
#endif
#if defined(VW_HAVE_PKG_GDAL)
    CREAT("tif",        GDAL)
    CREAT("tiff",       GDAL)
    CREAT("image/tiff", GDAL)
#endif
#if defined(VW_HAVE_PKG_OPENEXR)
    CREAT("exr",        OpenEXR)
    CREAT("image/exr",  OpenEXR)
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
