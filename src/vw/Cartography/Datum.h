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


#ifndef __VW_CARTOGRAPHY_DATUM_H__
#define __VW_CARTOGRAPHY_DATUM_H__

#include <string>
#include <ostream>
#include <cmath>

#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>

#include <ogr_spatialref.h>

/// \file Datum.h Planetary ellipsoidal coordinate system.
/// It is a wrapper around the GDAL OGRSpatialReference class.

namespace vw {
namespace cartography {

  /// A geodetic datum, i.e. a reference ellipsoid coordinate system
  /// for a planetary body.
  /// To express a spherical datum, set the semi-major axis equal to
  /// the semi-minor axis.  All angles are measured in degrees, and
  /// all distances are measured in meters.  This class incorporates a
  /// prime meridian offset, which is not usually strictly considered
  /// part of the datum but which has no better place to be right now.
  class Datum {
    std::string m_name;
    std::string m_spheroid_name;
    std::string m_meridian_name;
    double      m_semi_major_axis;
    double      m_semi_minor_axis;
    double      m_meridian_offset;
    OGRSpatialReference m_ogr_datum;
    std::string m_proj_str;

  public:
    /// The default constructor creates a WGS84 datum.
    Datum() {
      set_well_known_datum("WGS84");
    }

    /// Constructs a well-known datum by name.
    /// - Supported names: WGS84, WGS72, NAD83, NAD27, D_MOON, D_MARS, MOLA.
    /// - Note that this class does not fully support Earth datums other than WGS84!
    Datum(std::string const& name) {
      set_well_known_datum(name);
    }

    /// This constructor allows the user to create a custom datum.
    Datum(std::string const& name,
          std::string const& spheroid_name,
          std::string const& meridian_name,
          double             semi_major_axis,
          double             semi_minor_axis,
          double             meridian_offset);

    // A wrapper around the GDAL tool for extracting the datum information
    void set_datum_from_spatial_ref(OGRSpatialReference const& gdal_spatial_ref);
  
    // Set the datum using a WKT string
    void set_wkt(std::string const& wkt);

    // Set the datum using a Proj.4 string. It is preferred to set the datum
    // using a WKT.
    void set_datum_from_proj_str(std::string const& proj_str);
    
    /// See Datum(name)
    void set_well_known_datum(std::string const& name);

    // Basic accessors
    std::string      & name()       { return m_name; }
    std::string const& name() const { return m_name; }

    std::string      & spheroid_name()       { return m_spheroid_name; }
    std::string const& spheroid_name() const { return m_spheroid_name; }

    std::string      & meridian_name()       { return m_meridian_name; }
    std::string const& meridian_name() const { return m_meridian_name; }

    void   set_semi_major_axis(double val);
    double semi_major_axis() const { return m_semi_major_axis; }

    void   set_semi_minor_axis(double val);
    double semi_minor_axis() const { return m_semi_minor_axis; }

    double &meridian_offset()       { return m_meridian_offset; }
    double  meridian_offset() const { return m_meridian_offset; }

    std::string      & proj4_str()       { return m_proj_str; }
    std::string const& proj4_str() const { return m_proj_str; }
    std::string get_wkt() const;

    // Use set_datum_from_spatial_ref() to set the OGRSpatialReference.
    // This will only fetch the underlying object but not set it.
    OGRSpatialReference const& ogr_datum() const { return m_ogr_datum; }
    
    /// Returns the radius (distance from center of the body) at the given lat/lon
    double radius(double lon, double lat) const;

    /// Return geocentric latitude corresponding to geodetic lat:
    /// - Geocentric latitude is measured with body's center rather than tangent to surface.
    /// - For spherical datums these values are identical.
    double geocentric_latitude(double lat) const;

    /// Return radius of curvature in the prime vertical.
    double radius_of_curvature(double lon, double lat) const;

    /// Return distance from the center of the Earth (lat,lon are geodetic,
    //// alt is ellipsoidal or geodetic height).
    double geocentric_radius(double lon, double lat, double alt = 0.0) const;

    double inverse_flattening() const;

    // NED to ECEF transform. This takes into account that the planet is an
    // ellipsoid and the curent elevation above it. 
    Vector3 geodetic_to_cartesian( Vector3 const& llh ) const;

    /// Return the rotation matrix for converting from NED to ECEF.
    /// This can handle the planet being an ellipsoid.
    Matrix3x3 lonlat_to_ned_matrix(Vector3 const& llh) const;

    Vector3 cartesian_to_geodetic( Vector3 const& xyz ) const;
  };

  std::ostream& operator<<(std::ostream& os, const Datum& datum);

  // A function which reads a datum from a string. It reverses
  // what operator<< does.
  bool read_datum_from_str(std::string const& str, Datum & datum);
  
  // Free associated functions
  vw::Vector3 datum_intersection(double semi_major_axis, double semi_minor_axis,
                             vw::Vector3 camera_ctr, vw::Vector3 camera_vec);
  vw::Vector3 datum_intersection( Datum const& datum,
                              vw::Vector3 camera_ctr, vw::Vector3 camera_vec );
  
}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_DATUM_H__
