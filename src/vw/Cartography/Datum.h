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

#include <vw/Math/Vector.h>
#include <vw/Image/PerPixelViews.h>

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
    double const& semi_major_axis() const { return m_semi_major_axis; }

    void set_semi_minor_axis(double val);
    double const &semi_minor_axis() const { return m_semi_minor_axis; }

    double &meridian_offset() { return m_meridian_offset; }
    double const& meridian_offset() const { return m_meridian_offset; }


    std::string proj4_str() const { return m_proj_str; }

    double radius(double lat, double lon) const;

    double inverse_flattening() const;

    Vector3 geodetic_to_cartesian( Vector3 const& p ) const;
    Vector3 cartesian_to_geodetic( Vector3 const& p ) const;
  };

  std::ostream& operator<<(std::ostream& os, const Datum& datum);


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

  template <class ElemT>
  class AddDatumFunctor : public UnaryReturnSameType {
    Datum m_datum;

  public:
    AddDatumFunctor(Datum const& datum) : m_datum(datum) {}    
    Vector<ElemT,3> operator()(Vector<ElemT,3> const& p) const {
      if (p != Vector<ElemT,3>()) 
        return Vector<ElemT,3>(p[0], p[1], p[2]+m_datum.radius(p[0], p[1]));
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
  
  /// The inverse of subtract_datum(), above.
  template <class ImageT>
  UnaryPerPixelView<ImageT, AddDatumFunctor<typename ImageT::pixel_type::value_type> >
  inline add_datum( ImageViewBase<ImageT> const& image, Datum const& datum) {
    typedef typename ImageT::pixel_type::value_type vector_value_type;
    return UnaryPerPixelView<ImageT,AddDatumFunctor<vector_value_type> >( image.impl(), AddDatumFunctor<vector_value_type>(datum) );
  }

}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_DATUM_H__

