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

#include <string>
#include <ostream>
#include <math.h>

namespace vw {
namespace cartography {

  /// Class for a bi-axial ellipsoid geodetic datum.
  ///
  /// To express a spherical datum, set the semi-major axis equal to
  /// the semi-minor axis.
  class Datum {
    std::string m_name;
    std::string m_spheroid_name;
    std::string m_meridian_name;
    double m_semi_major_axis;
    double m_semi_minor_axis;
    double m_meridian_offset;       /// given in angular units
    std::string m_proj_str;

  public:
    /// The default constructor creates a WGS84 datum.
    Datum() {
      set_well_known_datum("WGS84");
    }

    /// This constructor allows the user to create a custom datum.
    Datum(std::string const& name,
          std::string const& spheroid_name,
          std::string const& meridian_name,
          double semi_major_axis,
          double semi_minor_axis,
          double meridian_offset) : m_name(name),
                                    m_spheroid_name(spheroid_name),
                                    m_meridian_name(spheroid_name),
                                    m_semi_major_axis(semi_major_axis),
                                    m_semi_minor_axis(semi_minor_axis),
                                    m_meridian_offset(meridian_offset) {
      std::ostringstream strm;
      strm << "+a=" << semi_major_axis << " +b=" << semi_minor_axis;
      m_proj_str = strm.str();
    }

    /// Options include: WGS84, WGS72, NAD27, or NAD83. 
    void set_well_known_datum(std::string const& name) {
      m_meridian_name = "Grenwich";
      m_meridian_offset = 0;
      if (name == "WGS84") {        
        m_name = "WGS_1984";
        m_spheroid_name="WGS 84";
        m_semi_major_axis = 6378137.0;
        m_semi_minor_axis = 6356752.3;
        m_proj_str = "+ellps=WGS84 +datum=WGS84";

      } else if (name == "WGS72") {
        m_name="WGS_1972";
        m_spheroid_name="WGS 72";
        m_semi_major_axis = 6378135.0;
        m_semi_minor_axis = 6356750.5;
        m_proj_str = "+ellps=WGS72 +towgs84=0,0,4.5,0,0,0.554,0.2263";

      } else if (name == "NAD83") {
        m_name="North_American_Datum_1983";
        m_spheroid_name="GRS 1980";
        m_semi_major_axis = 6378137;
        m_semi_minor_axis = 6356752.3;
        m_proj_str = "+ellps=GRS80 +datum=NAD83";

      } else if (name == "NAD27") {
        m_name="North_American_Datum_1927";
        m_spheroid_name="Clarke 1866";
        m_semi_major_axis = 6378206.4;
        m_semi_minor_axis = 6356583.8;
        m_proj_str = "+ellps=clrk66 +datum=NAD27";
      }
    }
    
    std::string name() { return m_name; }
    std::string const& name() const { return m_name; }

    std::string &spheroid_name() { return m_spheroid_name; }
    std::string const& spheroid_name() const { return m_spheroid_name; }

    std::string &meridian_name() { return m_meridian_name; }
    std::string const& meridian_name() const { return m_meridian_name; }

    void set_semi_major_axis(double val) { 
      m_semi_major_axis = val;  
      std::ostringstream strm;
      strm << "+a=" << m_semi_major_axis << " +b=" << m_semi_minor_axis;
      m_proj_str = strm.str();
    }
    double const& semi_major_axis() const { return m_semi_major_axis; }

    void set_semi_minor_axis(double val) { 
      m_semi_minor_axis = val;  
      std::ostringstream strm;
      strm << "+a=" << m_semi_major_axis << " +b=" << m_semi_minor_axis;
      m_proj_str = strm.str();
    }
    double const &semi_minor_axis() const { return m_semi_minor_axis; }

    double &meridian_offset() { return m_meridian_offset; }
    double const& meridian_offset() const { return m_meridian_offset; }

    std::string proj4_str() const { return m_proj_str; }

    double radius(double lon, double lat) const {

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

  static std::ostream& operator<<(std::ostream& os, const Datum& datum) {
    os << "Geodeditic Datum --> Name: " << datum.name() << "  Spheroid: " << datum.spheroid_name() 
       << "  Semi-major: " << datum.semi_major_axis()
       << "  Semi-minor: " << datum.semi_minor_axis()
       << "  Meridian: " << datum.meridian_name()
       << "  at " << datum.meridian_offset();
    return os;
  }


  template <class ElemT>
  class SubtractDatumFunctor : public UnaryReturnSameType {
    Datum m_datum;

  public:
    SubtractDatumFunctor(Datum const& datum) : m_datum(datum) {}    
    Vector<ElemT,3> operator()(Vector<ElemT,3> const& p) const {
      if (p != Vector<ElemT,3>()) 
        return Vector<ElemT,3>(p[0], p[1], p[2]-m_datum.radius(p[0], p[1]));
      else 
        return p;
    }
  };
  

  /// Takes an ImageView of Vector<ElemT,3> in lon/lat/radius and
  /// returns an ImageView of vectors that contains the lol, lat, and
  /// altitude of that point referenced to the given datum.  For
  /// consistency with cartographic convention, angular values must be
  /// given in degrees rather than radians.
  ///
  /// Note: The following assumes latitude is measured from the
  /// equatorial plane with north positive. This is different than
  /// normal spherical coordinate conversion where the equivalent
  /// angle is measured from the positive z axis.
  //
  /// Note: notice that the order of the returned triple is longitude,
  /// latitude, radius.  This ordering of lon/lat is consistent with
  /// the notion of horizontal (x) and vertical (y) coordinates in an
  /// image.
  template <class ImageT>
  UnaryPerPixelView<ImageT, SubtractDatumFunctor<typename ImageT::pixel_type::value_type> >
  inline subtract_datum( ImageViewBase<ImageT> const& image, Datum const& datum) {
    typedef typename ImageT::pixel_type::value_type vector_value_type;
    return UnaryPerPixelView<ImageT,SubtractDatumFunctor<vector_value_type> >( image.impl(), SubtractDatumFunctor<vector_value_type>(datum) );
  }
  

}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_DATUM_H__

