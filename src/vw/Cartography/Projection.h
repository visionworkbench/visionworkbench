// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
