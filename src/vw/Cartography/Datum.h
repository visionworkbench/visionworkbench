// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_CARTOGRAPHY_DATUM_H__
#define __VW_CARTOGRAPHY_DATUM_H__

#include <string>
#include <ostream>
#include <cmath>

#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>

// For an excellent discussion of the various concepts, terms, and
// issues involved in global coordinate systems, see:
// http://www.posc.org/Epicentre.2_2/DataModel/ExamplesofUsage/eu_cs.html

namespace vw {
namespace cartography {

  /// A geodetic datum, i.e. a reference ellipsoid coordinate system
  /// for a planetary body.  This implementation assumes a relatively
  /// modern notion of a datum, ie. a geocentric bi-axial ellipsoidal
  /// model.
  ///
  /// To express a spherical datum, set the semi-major axis equal to
  /// the semi-minor axis.  All angles are measured in degrees, and
  /// all distances are measured in meters.  This class incorporates a
  /// prime meridian offset, which is not usually strictly considered
  /// part of the datum but which has no better place to be right now.
  class Datum {
    std::string m_name;
    std::string m_spheroid_name;
    std::string m_meridian_name;
    double m_semi_major_axis;
    double m_semi_minor_axis;
    double m_meridian_offset;
    bool m_geocentric;
    std::string m_proj_str;

  public:
    /// The default constructor creates a WGS84 datum.
    Datum() {
      set_well_known_datum("WGS84");
    }

    /// Constructs a well-known datum by name.
    Datum( std::string const& name ) {
      set_well_known_datum( name );
    }

    /// This constructor allows the user to create a custom datum.
    Datum(std::string const& name,
          std::string const& spheroid_name,
          std::string const& meridian_name,
          double semi_major_axis,
          double semi_minor_axis,
          double meridian_offset);

    /// Options include: WGS84, WGS72, NAD27, or NAD83.
    void set_well_known_datum(std::string const& name);

    std::string &name() { return m_name; }
    std::string const& name() const { return m_name; }

    std::string &spheroid_name() { return m_spheroid_name; }
    std::string const& spheroid_name() const { return m_spheroid_name; }

    std::string &meridian_name() { return m_meridian_name; }
    std::string const& meridian_name() const { return m_meridian_name; }

    void set_semi_major_axis(double val);
    double semi_major_axis() const { return m_semi_major_axis; }

    void set_semi_minor_axis(double val);
    double semi_minor_axis() const { return m_semi_minor_axis; }

    double &meridian_offset() { return m_meridian_offset; }
    double meridian_offset() const { return m_meridian_offset; }

    void set_geocentric(bool val);
    bool geocentric() const { return m_geocentric; }

    std::string &proj4_str() { return m_proj_str; }
    std::string const& proj4_str() const { return m_proj_str; }

    double radius(double lon, double lat) const;

    /// return geocentric latitude corresponding to geodetic lat:
    double geocentric_latitude(double lat) const;

    /// return radius of curvature in the prime vertical.
    double radius_of_curvature(double lon, double lat) const;

    /// return distance from the center of the Earth (lat,lon are geodetic,
    //// alt is ellipsoidal or geodetic height).
    double geocentric_radius(double lon, double lat, double alt = 0.0) const;

    double inverse_flattening() const;

    /// Return cartesian (ECEF) coordinates of geodetic coordinates p
    Vector3 geodetic_to_cartesian( Vector3 const& p ) const;

    /// Return rotation matrix for converting NED vectors
    /// to ECEF vectors.
    Matrix3x3 ecef_to_ned_matrix( Vector3 const& p) const;

    Vector3 cartesian_to_geodetic( Vector3 const& p ) const;
  };

  std::ostream& operator<<(std::ostream& os, const Datum& datum);


}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_DATUM_H__

