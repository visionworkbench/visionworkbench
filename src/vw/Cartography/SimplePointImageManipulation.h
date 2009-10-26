// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_CARTOGRAPHY_SIMPLEPOINTIMAGEMANIPLULATION_H__
#define __VW_CARTOGRAPHY_SIMPLEPOINTIMAGEMANIPLULATION_H__

#include <cmath>
#include <vw/Math/Vector.h>
#include <vw/Math/Functors.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/ImageViewBase.h>

/// \file SimplePointImageManipulation.h
///
/// Contains routines for manipulating ImageViews with a pixel type of
/// Vector3 that contain xyz or lon, lat, altitude triples. These functions
/// should stand alone, and not use GeoReferences.

namespace vw {
namespace cartography {

  // ---------------- XYZ to LON LAT ALT CONVERSION ------------------


  class XYZtoLonLatRadFunctor : public UnaryReturnSameType {
    bool m_east_positive;
  public:
    XYZtoLonLatRadFunctor(bool east_positive = true) : m_east_positive(east_positive) {}

    template <class T>
    T operator()(T const& p) const {
      return this->apply(p, m_east_positive);
    }

    template <class T>
    static inline T apply(T const& p, bool east_positive = true)  {

      // Deal with "missing pixels"
      if (p == T()) { return p; }

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

      return T (lon * 180.0 / M_PI, lat * 180.0 / M_PI, radius);
    }
  };

  class LonLatRadToXYZFunctor : public UnaryReturnSameType {
    bool m_east_positive;
  public:
    LonLatRadToXYZFunctor(bool east_positive = true) : m_east_positive(east_positive) {}

    // Convert from lon, lat, radius to x,y,z:
    //
    // x = r * cos(longitude)
    // y = r * sin(longitude)
    // z = r * sin(latitude)
    //
    template <class T>
    T operator()(T const& p) const {
      return this->apply(p, m_east_positive);
    }

    template <class T>
    static inline T apply(T const& p, bool east_positive = true)  {
      typename T::value_type z = p(2) * sin(p(1)*M_PI/180.0);
      typename T::value_type sqrt_x_sqr_plus_y_sqr = p(2) * cos(p(1)*M_PI/180.0);

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
  UnaryPerPixelView<ImageT, XYZtoLonLatRadFunctor>
  inline xyz_to_lon_lat_radius( ImageViewBase<ImageT> const& image, bool east_positive = true ) {
    return UnaryPerPixelView<ImageT,XYZtoLonLatRadFunctor>( image.impl(), XYZtoLonLatRadFunctor(east_positive) );
  }

  template <class ElemT>
  inline Vector<ElemT,3> xyz_to_lon_lat_radius( Vector<ElemT,3> const& xyz, bool east_positive = true ) {
    return XYZtoLonLatRadFunctor::apply(xyz, east_positive);
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
  UnaryPerPixelView<ImageT, LonLatRadToXYZFunctor>
  inline lon_lat_radius_to_xyz( ImageViewBase<ImageT> const& image, bool east_positive = true ) {
    return UnaryPerPixelView<ImageT,LonLatRadToXYZFunctor>( image.impl(), LonLatRadToXYZFunctor(east_positive) );
  }

  template <class ElemT>
  inline Vector<ElemT,3> lon_lat_radius_to_xyz( Vector<ElemT,3> const& lon_lat_alt, bool east_positive = true ) {
    return LonLatRadToXYZFunctor::apply(lon_lat_alt, east_positive);
  }

}} // namespace vw::cartography

#endif
