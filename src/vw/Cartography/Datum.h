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
#ifndef __VW_CARTOGRAPHY_DATUM_H__
#define __VW_CARTOGRAPHY_DATUM_H__

namespace vw {
namespace cartography {

  /// Abstract base class for a geodetic datum.  Subclasses of this class 
  /// must implement the radius() method for computing the radius of the 
  /// geoid at a given lat/lon.
  ///
  /// To express a spherical datum, set the semi-major axis equal to
  /// the semi-minor axis.
  class GeoDatum {
    std::string m_name;
    std::string m_spheroid_name;
    std::string m_meridian_name;
    double m_semi_major_axis;
    double m_semi_minor_axis;
    double m_meridian_offset;       /// given in angular units

  public:
    GeoDatum() : m_name("Unknown Datum"),
                 m_spheroid_name("Unknown Spheroid"),
                 m_meridian_name("Unknown Meridian"),
                 m_semi_major_axis(0),
                 m_semi_minor_axis(0),
                 m_meridian_offset(0) {}

    GeoDatum(std::string const& name,
             std::string const& spheroid_name,
             std::string const& meridian_name,
             double semi_major_axis,
             double semi_minor_axis,
             double meridian_offset) : m_name(name),
                                       m_spheroid_name(spheroid_name),
                                       m_meridian_name(spheroid_name),
                                       m_semi_major_axis(semi_major_axis),
                                       m_semi_minor_axis(semi_minor_axis),
                                       m_meridian_offset(meridian_offset) {}

    std::string &name() { return m_name; }
    std::string const& name() const { return m_name; }

    std::string &spheroid_name() { return m_spheroid_name; }
    std::string const& spheroid_name() const { return m_spheroid_name; }

    std::string &meridian_name() { return m_meridian_name; }
    std::string const& meridian_name() const { return m_meridian_name; }

    double &semi_major_axis() { return m_semi_major_axis; }
    double const& semi_major_axis() const { return m_semi_major_axis; }

    double &semi_minor_axis() { return m_semi_minor_axis; }
    double const &semi_minor_axis() const { return m_semi_minor_axis; }

    double &meridian_offset() { return m_meridian_offset; }
    double const& meridian_offset() const { return m_meridian_offset; }

    double radius(double lat, double lon) const {
      // Optimize in the case of spherical datum
      if (m_semi_major_axis == m_semi_minor_axis) {
        return m_semi_major_axis;
      } 
      
      // Bi-axial Ellpisoid datum
      double a = m_semi_major_axis;
      double b = m_semi_minor_axis;
      double t = atan((a/b) * tan(lat * M_PI / 180.0));
      double x = a * cos(t);
      double y = b * sin(t);
      return sqrt(x*x + y*y);
    }

    inline double inverse_flattening() const {
      return 1.0 / (1.0 - m_semi_minor_axis / m_semi_major_axis);
    }

    
  };

  static std::ostream& operator<<(std::ostream& os, const GeoDatum& datum) {
    os << "Geodeditic Datum --> Name: " << datum.name() << "  Spheroid: " << datum.spheroid_name() 
       << "  Semi-major: " << datum.semi_major_axis()
       << "  Semi-minor: " << datum.semi_minor_axis()
       << "  Meridian: " << datum.meridian_name()
       << "  at " << datum.meridian_offset();
    return os;
  }

}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_DATUM_H__

