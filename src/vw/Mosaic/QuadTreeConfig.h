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


#ifndef __VW_MOSAIC_QUADTREECONFIG_H__
#define __VW_MOSAIC_QUADTREECONFIG_H__

#include <vw/Cartography/GeoReference.h>
#include <boost/shared_ptr.hpp>

namespace vw {
namespace mosaic {

  class QuadTreeGenerator;

  /// ???
  class QuadTreeConfig {
  public:
    virtual ~QuadTreeConfig() {}
    /// 
    virtual void configure( QuadTreeGenerator& qtree ) const = 0;
    /// 
    virtual cartography::GeoReference output_georef(uint32 xresolution, uint32 yresolution = 0) = 0;
    
    /// Creates a new QuadTreeConfig object of the specified type
    static boost::shared_ptr<QuadTreeConfig> make(const std::string& type);
  };

}} // namespace vw::mosaic

#endif
