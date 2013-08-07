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


#ifndef __VW_CARTOGRAPHY_GEOREFERENCEHELPERPDS_H__
#define __VW_CARTOGRAPHY_GEOREFERENCEHELPERPDS_H__

#include <vw/FileIO/DiskImageResourcePDS.h>
#include <vw/Cartography/GeoReference.h>

// Boost
#include <boost/algorithm/string.hpp>

/// \file GeoReferenceResourcePDS.h Functions to read/write a \ref vw::cartography::GeoReference from PDS.

namespace vw {
namespace cartography {

  bool read_pds_georeference( GeoReference& georef, DiskImageResourcePDS const& resource );
  // We do not support writing PDS images at this time.
  // void write_pds_georeference( DiskImageResourcePDS& resource, GeoReference const& georef );

}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_GEOREFERENCEHELPERPDS_H__
