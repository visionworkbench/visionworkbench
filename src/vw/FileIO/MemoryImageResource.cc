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


#include <vw/FileIO/MemoryImageResource.h>

#  include <vw/FileIO/MemoryImageResourceJPEG.h>

#  include <vw/FileIO/MemoryImageResourcePNG.h>

#  include <vw/FileIO/MemoryImageResourceGDAL.h>


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

  SrcMemoryImageResource* SrcMemoryImageResource::open(const std::string& type, const uint8* data, size_t len) {
    boost::shared_array<const uint8> p(data, NOP());
    return SrcMemoryImageResource::open(type, p, len);
  }

  SrcMemoryImageResource* SrcMemoryImageResource::open(const std::string& type, boost::shared_array<const uint8> data, size_t len) {

    std::string ctype = clean_type(type);
    
    if (ctype == "jpg" || ctype == "jpeg" || ctype == "image/jpeg")
      return new SrcMemoryImageResourceJPEG(data, len);
    
    if (ctype == "png" || ctype == "image/png")
      return new SrcMemoryImageResourcePNG(data, len);

    if (ctype == "tif" || ctype == "tiff" || ctype == "image/tiff")
      return new SrcMemoryImageResourceGDAL(data, len);
    
    
    vw_throw(NoImplErr() << "Unsupported file format: " << type);
    return NULL;
  }

  DstMemoryImageResource* DstMemoryImageResource::create(const std::string& type, const ImageFormat& format) {

    std::string ctype = clean_type(type);
    
    if (ctype == "jpg" || ctype == "jpeg" || ctype == "image/jpeg")
      return new DstMemoryImageResourceJPEG(format);
    
    if (ctype == "png" || ctype == "image/png")
      return new DstMemoryImageResourcePNG(format);

    if (ctype == "tif" || ctype == "tiff" || ctype == "image/tiff")
      return new DstMemoryImageResourceGDAL(format);
    
    
    vw_throw(NoImplErr() << "Unsupported file format: " << type);
    return NULL;
  }

} // namespace vw
