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


#include <vw/vw_config.h>
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

// Turn off warnings about things we can't control
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#include <boost/assign/list_of.hpp>
#include <boost/lambda/construct.hpp>
#include <boost/function.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/trim.hpp>
#pragma GCC diagnostic pop

#include <map>

namespace {
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

    std::string ctype = clean_type(type);
    
#if defined(VW_HAVE_PKG_JPEG)
    if (ctype == "jpg" || ctype == "jpeg" || ctype == "image/jpeg")
      return new SrcMemoryImageResourceJPEG(data, len);
#endif
    
#if defined(VW_HAVE_PKG_PNG)
    if (ctype == "png" || ctype == "image/png")
      return new SrcMemoryImageResourcePNG(data, len);
#endif

#if defined(VW_HAVE_PKG_GDAL)
    if (ctype == "tif" || ctype == "tiff" || ctype == "image/tiff")
      return new SrcMemoryImageResourceGDAL(data, len);
#endif
    
#if defined(VW_HAVE_PKG_OPENEXR)
    if (ctype == "exr" || ctype == "image/exr")
      return new SrcMemoryImageResourceOpenEXR(data, len);
#endif
    
    vw_throw( NoImplErr() << "Unsupported file format: " << type );
    return NULL;
  }

  DstMemoryImageResource* DstMemoryImageResource::create( const std::string& type, const ImageFormat& format ) {

    std::string ctype = clean_type(type);
    
#if defined(VW_HAVE_PKG_JPEG)
    if (ctype == "jpg" || ctype == "jpeg" || ctype == "image/jpeg")
      return new DstMemoryImageResourceJPEG(format);
#endif
    
#if defined(VW_HAVE_PKG_PNG)
    if (ctype == "png" || ctype == "image/png")
      return new DstMemoryImageResourcePNG(format);
#endif

#if defined(VW_HAVE_PKG_GDAL)
    if (ctype == "tif" || ctype == "tiff" || ctype == "image/tiff")
      return new DstMemoryImageResourceGDAL(format);
#endif
    
#if defined(VW_HAVE_PKG_OPENEXR)
    if (ctype == "exr" || ctype == "image/exr")
      return new DstMemoryImageResourceOpenEXR(format);
#endif
    
    vw_throw( NoImplErr() << "Unsupported file format: " << type );
    return NULL;
  }

} // namespace vw
