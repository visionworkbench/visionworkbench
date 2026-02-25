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


#include <vw/Math/Matrix.h>
#include <vw/Core/Exception.h>
#include <vw/Core/Log.h>

#include <vw/FileIO/DiskImageResourcePDS.h>
#include <vw/Cartography/GeoReferenceResourcePDS.h>
#include <vw/Cartography/Datum.h>
#include <vw/Cartography/GeoReference.h>

// For boost::split
#include <boost/algorithm/string.hpp>

namespace vw {
namespace cartography {

  bool read_pds_georeference(GeoReference& georef, DiskImageResourcePDS const& resource) {

    std::string key, value;
    bool failed = false;

    Datum datum;
    datum.name() = "PDS Datum";
    datum.spheroid_name() = "PDS Spheroid";
    // Collect datum information
    key = "A_AXIS_RADIUS";
    if (resource.query(key,value)) { datum.set_semi_major_axis(atof(value.c_str()) * 1000); } else {failed = true; }
    key = "B_AXIS_RADIUS";
    if (resource.query(key,value)) { datum.set_semi_minor_axis(atof(value.c_str()) * 1000); } else {failed = true; }
    georef.set_datum(datum);

    // Collect projection information
    key = "MAP_PROJECTION_TYPE";
    if (resource.query(key,value)) {
      if (value != "\"SIMPLE CYLINDRICAL\"") {
        vw_out(ErrorMessage, "console") << "Unsupported map projection type in PDS header.";
        vw_out(ErrorMessage, "cartography") << "Unsupported map projection type in PDS header.";
        failed = true;
      }
      datum.set_semi_minor_axis(atof(value.c_str()) * 1000);
    } else {failed = true; }


    // Set affine transform
    Matrix<double,3,3> transform;
    transform.set_identity();

    key = "WESTERNMOST_LONGITUDE";
    if (resource.query(key,value)) { transform(0,2) = atof(value.c_str()); } else {failed = true; }
    key = "MAXIMUM_LATITUDE";
    if (resource.query(key,value)) { transform(1,2) = atof(value.c_str()); } else {failed = true; }
    key = "MAP_RESOLUTION";
    if (resource.query(key,value)) {
      transform(0,0) = 1/atof(value.c_str());
      transform(1,1) = -1/atof(value.c_str());
    } else {failed = true; }

    georef.set_transform(transform);

    return !failed;
  }

}} // namespace vw::cartography
