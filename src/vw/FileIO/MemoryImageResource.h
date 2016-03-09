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


/// \file MemoryImageResource.h An abstract base class referring to an encoded image in memory.

#ifndef __VW_FILEIO_MEMORYIMAGERESOURCE_H__
#define __VW_FILEIO_MEMORYIMAGERESOURCE_H__

#include <vw/Image/ImageResource.h>

namespace vw {

  class SrcMemoryImageResource : public SrcImageResource {
    public:
      // constructs the appropriate subclass for the type
      // Does not take ownership of data
      static SrcMemoryImageResource* open( const std::string& type, const uint8* data, size_t len);
      // Takes ownership of data
      static SrcMemoryImageResource* open( const std::string& type, boost::shared_array<const uint8> data, size_t len);
  };

  class DstMemoryImageResource : public DstImageResource {
    public:
      // constructs the appropriate subclass for the type
      static DstMemoryImageResource* create( const std::string& type, const ImageFormat& format );
      virtual const uint8* data() const = 0;
      virtual size_t size() const = 0;
  };


} // namespace vw

#endif
