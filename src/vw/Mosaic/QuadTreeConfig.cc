// __BEGIN_LICENSE__
//  Copyright (c) 2006-2026, United States Government as represented by the
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

#include <vw/Mosaic/QuadTreeConfig.h>
#include <vw/Mosaic/KMLQuadTreeConfig.h>

#include <boost/algorithm/string.hpp>

using namespace vw;
using namespace vw::mosaic;

boost::shared_ptr<QuadTreeConfig> QuadTreeConfig::make(const std::string& type) {
  typedef boost::shared_ptr<QuadTreeConfig> ptr_t;
  std::string utype = boost::to_upper_copy(type);

  if (utype == "KML")
    return ptr_t(new KMLQuadTreeConfig());
  else
    vw_throw(NoImplErr() << "Unknown quad tree type: " << utype);
}
