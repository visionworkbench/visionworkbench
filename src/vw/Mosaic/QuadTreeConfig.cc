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


#include <vw/Mosaic/QuadTreeConfig.h>
#include <vw/Mosaic/CelestiaQuadTreeConfig.h>
#include <boost/algorithm/string.hpp>
#include <vw/Mosaic/GigapanQuadTreeConfig.h>
#include <vw/Mosaic/GMapQuadTreeConfig.h>
#include <vw/Mosaic/KMLQuadTreeConfig.h>
#include <vw/Mosaic/TMSQuadTreeConfig.h>
#include <vw/Mosaic/UniviewQuadTreeConfig.h>

using namespace vw;
using namespace vw::mosaic;

boost::shared_ptr<QuadTreeConfig> QuadTreeConfig::make(const std::string& type) {
  typedef boost::shared_ptr<QuadTreeConfig> ptr_t;
  std::string utype = boost::to_upper_copy(type);

  if (utype == "CELESTIA")
    return ptr_t(new CelestiaQuadTreeConfig());
  else if (utype == "GIGAPAN" || utype == "GIGAPAN_NOPROJ")
    return ptr_t(new GigapanQuadTreeConfig());
  else if (utype == "GMAP")
    return ptr_t(new GMapQuadTreeConfig());
  else if (utype == "KML")
    return ptr_t(new KMLQuadTreeConfig());
  else if (utype == "TMS")
    return ptr_t(new TMSQuadTreeConfig());
  else if (utype == "UNIVIEW")
    return ptr_t(new UniviewQuadTreeConfig());
  else
    vw_throw(NoImplErr() << "Unknown quad tree type: " << utype);
}
