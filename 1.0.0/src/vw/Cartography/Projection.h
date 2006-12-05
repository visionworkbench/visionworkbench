// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__
#ifndef __VW_CARTOGRAPHY_PROJECTION_H__
#define __VW_CARTOGRAPHY_PROJECTION_H__

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

  static std::ostream& operator<<(std::ostream& os, const GeoProjection& projection) {
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
