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
#ifndef __VW_CARTOGRAPHY_POINTIMAGE_MANIPLULATION_H__
#define __VW_CARTOGRAPHY_POINTIMAGE_MANIPLULATION_H__

#include <string>
#include <ostream>
#include <math.h>

#include <vw/Math/Vector.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Cartography/Datum.h>
#include <vw/Cartography/GeoReferenceBase.h>

/// \file PointImageManipulation.h
/// 
/// Contains routines for manipulating ImageViews with a pixel type of
/// Vector3 that contain xyz or lon, lat, altitude triples.

namespace vw {
namespace cartography {

  // ---------------- XYZ to LON LAT ALT CONVERSION ------------------

  
  template <class ElemT>
  class XYZtoLonLatFunctor : public UnaryReturnSameType {
    bool m_east_positive;
  public:
    XYZtoLonLatFunctor(bool east_positive = true) : m_east_positive(east_positive) {}
    
    Vector<ElemT,3> operator()(Vector<ElemT,3> const& p) const {
      return this->apply(p, m_east_positive);
    } 

    static inline Vector<ElemT,3> apply(Vector<ElemT,3> const& p, bool east_positive = true)  {

      // Deal with "missing pixels"
      if (p == Vector<ElemT,3>()) { return p; }
      
      double radius = norm_2(p);
      double sin_lat = p.z() / radius;
      double lat = asin(sin_lat);

      double lon;
      if (east_positive) 
        lon = atan2(p.y(), p.x());
      else // West positive longitude
        lon = atan2(-p.y(), p.x()); 

      // For consistency-sake, we always return a longitude in the range +/-180.
      if (lon > M_PI)
        lon -= 2*M_PI;
      if (lon < -M_PI)
        lon += 2*M_PI;

      return Vector<ElemT,3> (lon * 180.0 / M_PI, lat * 180.0 / M_PI, radius);
    }
  };

  template <class ElemT>
  class LonLatToXYZFunctor : public UnaryReturnSameType {
    bool m_east_positive;
  public:
    LonLatToXYZFunctor(bool east_positive = true) : m_east_positive(east_positive) {}

    // Convert from lon, lat, radius to x,y,z:
    //
    // x = r * cos(longitude)
    // y = r * sin(longitude)
    // z = r * sin(latitude)
    //
    Vector<ElemT,3> operator()(Vector<ElemT,3> const& p) const {
      return this->apply(p, m_east_positive);
    }

    static inline Vector<ElemT,3> apply(Vector<ElemT,3> const& p, bool east_positive = true)  {
      ElemT z = p(2) * sin(p(1)*M_PI/180.0);
      ElemT sqrt_x_sqr_plus_y_sqr = p(2) * cos(p(1)*M_PI/180.0);

      if (east_positive) {
        return Vector3( sqrt_x_sqr_plus_y_sqr * cos(p(0)*M_PI/180.0),
                        sqrt_x_sqr_plus_y_sqr * sin(p(0)*M_PI/180.0),
                        z);
      } else {
        return Vector3( sqrt_x_sqr_plus_y_sqr * cos(-p(0)*M_PI/180.0),
                        sqrt_x_sqr_plus_y_sqr * sin(-p(0)*M_PI/180.0),
                        z );
      }
    }
  };
  

  /// Takes an ImageView of Vector<ElemT,3> in cartesian 3 space and
  /// returns a ImageView of vectors that contains the lat, lon, and
  /// radius of that point.  For consistency with cartographic
  /// convention, angular values are return in degrees rather than
  /// radians.
  ///
  /// Note: The following assumes latitude is measured from the
  /// equatorial plane with north positive. This is different than
  /// normal spherical coordinate conversion where the equivalent
  /// angle is measured from the positive z axis.
  ///
  /// Note: notice that the order of the returned triple is longitude,
  /// latitude, radius.  This ordering of lon/lat is consistent with
  /// the notion of horizontal (x) and vertical (y) coordinates in an
  /// image.
  template <class ImageT>
  UnaryPerPixelView<ImageT, XYZtoLonLatFunctor<typename ImageT::pixel_type::value_type> >
  inline xyz_to_lon_lat_radius( ImageViewBase<ImageT> const& image, bool east_positive = true ) {
    typedef typename ImageT::pixel_type::value_type vector_value_type;
    return UnaryPerPixelView<ImageT,XYZtoLonLatFunctor<vector_value_type> >( image.impl(), XYZtoLonLatFunctor<vector_value_type>(east_positive) );
  }

  template <class ElemT>
  inline Vector<ElemT,3> xyz_to_lon_lat_radius( Vector<ElemT,3> const& xyz, bool east_positive = true ) {
    return XYZtoLonLatFunctor<double>::apply(xyz, east_positive);
  }
    

  /// Takes an ImageView of Vector<ElemT,3> that contains longitude,
  /// latitude, and radius, and an ImageView of vectors that are in
  /// cartesian 3-space.  For consistency with cartographic
  /// convention, angular values are expected to be in degrees rather
  /// than radians.
  ///
  /// Note: The following assumes latitude is measured from the
  /// equatorial plane with north positive. This is different than
  /// normal spherical coordinate conversion where the equivalent
  /// angle is measured from the positive z axis.
  ///
  /// Note: notice that the order of the returned triple is longitude,
  /// latitude, radius.  This ordering of lon/lat is consistent with
  /// the notion of horizontal (x) and vertical (y) coordinates in an
  /// image.
  template <class ImageT>
  UnaryPerPixelView<ImageT, LonLatToXYZFunctor<typename ImageT::pixel_type::value_type> >
  inline lon_lat_radius_to_xyz( ImageViewBase<ImageT> const& image, bool east_positive = true ) {
    typedef typename ImageT::pixel_type::value_type vector_value_type;
    return UnaryPerPixelView<ImageT,LonLatToXYZFunctor<vector_value_type> >( image.impl(), LonLatToXYZFunctor<vector_value_type>(east_positive) );
  }

  template <class ElemT>
  inline Vector<ElemT,3> lon_lat_radius_to_xyz( Vector<ElemT,3> const& xyz, bool east_positive = true ) {
    return LonLatToXYZFunctor<ElemT>::apply(xyz, east_positive);
  }

  // --------------------- CHANGE OF PROJECTION ----------------------

  /// Move from a source projection to a destination projection.  The
  /// altitude value (that is, the third value in the vector) is left
  /// untouched.
  template <class ElemT>
  class ReprojectPointFunctor : public UnaryReturnSameType {
    GeoReference m_src, m_dst;

  public:
    ReprojectPointFunctor(GeoReference const& src, GeoReference const& dst) : m_src(src), m_dst(dst) {}    
    Vector<ElemT,3> operator()(Vector<ElemT,3> const& p) const {
      Vector<ElemT,2> lon_lat = m_dst.lonlat_to_point(m_src.point_to_lonlat(subvector(p,0,2)));
      return Vector<ElemT,3>(lon_lat(0), lon_lat(1), 
                             p(2)+m_src.datum().radius(lon_lat(0), lon_lat(1))-m_dst.datum().radius(lon_lat(0), lon_lat(1)));
    }
  };

  // This version of the functor assumes that the point inputs are
  // unprojected (lon, lat, radius)
  template <class ElemT>
  class ProjectPointFunctor : public UnaryReturnSameType {
    GeoReference m_dst;

  public:
    ProjectPointFunctor(GeoReference const& dst) : m_dst(dst) {}    
    Vector<ElemT,3> operator()(Vector<ElemT,3> const& p) const {
      if (p == Vector3()) 
        return p;
      Vector<ElemT,2> lon_lat = m_dst.lonlat_to_point(subvector(p,0,2));
      return Vector<ElemT,3>(lon_lat(0), lon_lat(1), 
                             p(2)-m_dst.datum().radius(lon_lat(0), lon_lat(1)));
    }
  };


  /// Takes an ImageView of Vector<ElemT,3> in some source projected space
  /// with (lon,lat,alt) or (x,y,alt) and returns an ImageView of
  /// vectors that are in the destination projection.
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
  UnaryPerPixelView<ImageT, ReprojectPointFunctor<typename ImageT::pixel_type::value_type> >
  inline reproject_point_image( ImageViewBase<ImageT> const& image, GeoReference const& src_georef, GeoReference const& dst_georef) {
    typedef typename ImageT::pixel_type::value_type vector_value_type;
    return UnaryPerPixelView<ImageT,ReprojectPointFunctor<vector_value_type> >( image.impl(), ReprojectPointFunctor<vector_value_type>(src_georef, dst_georef) );
  }

  // This variant, which only accepts a destination projection,
  // assumes that the source points are [lon, lat, radius] values.
  template <class ImageT>
  UnaryPerPixelView<ImageT, ProjectPointFunctor<typename ImageT::pixel_type::value_type> >
  inline project_point_image( ImageViewBase<ImageT> const& image, GeoReference const& dst_georef) {
    typedef typename ImageT::pixel_type::value_type vector_value_type;
    return UnaryPerPixelView<ImageT,ProjectPointFunctor<vector_value_type> >( image.impl(), ProjectPointFunctor<vector_value_type>(dst_georef) );
  }

}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_POINTIMAGE_MANIPLULATION_H__

