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


#ifndef __VW_CARTOGRAPHY_PROJECTION_H__
#define __VW_CARTOGRAPHY_PROJECTION_H__

/// \file Projection.h Multiple classes for geographic projections.

namespace vw {
namespace cartography {

  /// Abstract base class for geographic projections
  struct GeoProjection {
    std::string name;
    std::string proj4_string;
    GeoProjection() : name("Geographic Projection"),
                      proj4_string("+proj=latlon") {}
    virtual ~GeoProjection() {};
  };

  inline std::ostream& operator<<(std::ostream& os, const GeoProjection& projection) {
    os << "Projection       --> Name: " << projection.name;
    return os;
  }

  /// Geographic Projection
  struct GeographicProjection : public GeoProjection {

    GeographicProjection() { name = "Geographic Projection"; }
  };

  /// Sinusoidal Projection
  struct SinusoidalProjection : public GeoProjection {

    SinusoidalProjection() : center_longitude(0), false_easting(0), false_northing(0) {
      name = "Sinusoidal Projection";
    }

    double center_longitude;
    double false_easting;
    double false_northing;
  };

  /// Cylindrical Equal Area Projection
  struct CEAProjection : public GeoProjection {

    CEAProjection() : std_p1(0), central_meridian(0),
                      false_easting(0), false_northing(0) {
      name = "Cylindrical Equal-Area Projection";
    }

    double std_p1;
    double central_meridian;
    double false_easting;
    double false_northing;
  };


}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_PROJECTION_H__
