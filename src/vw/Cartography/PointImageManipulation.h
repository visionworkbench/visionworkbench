// __BEGIN_LICENSE__
//  Copyright (c) 2006-2024, United States Government as represented by the
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


#ifndef __VW_CARTOGRAPHY_POINTIMAGEMANIPLULATION_H__
#define __VW_CARTOGRAPHY_POINTIMAGEMANIPLULATION_H__

#include <vw/Math/Vector.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/UtilityViews.h>
#include <vw/Cartography/GeoReference.h>


#include <cmath>
#include <vw/Math/Functors.h>
#include <vw/Image/ImageViewBase.h>

/// \file PointImageManipulation.h
///
/// Contains routines for manipulating ImageViews with a pixel type of
/// Vector3 that contain xyz or lon, lat, altitude triples.

namespace vw {
namespace cartography {

  // Functors. View operations are lower in this file.

  /// Functor to convert split col/row/alt values to lon/lat/alt using a georef.
  template <class PixelT>
  class DemToGeodetic: public ReturnFixedType<Vector3> {
    // Use an alias, which is apparently thread-safe, as a copy is very expensive
    // for dynamic CRS.
    GeoReference const& m_georef;
  public:
    DemToGeodetic(GeoReference const& georef): m_georef(georef) {}

    Vector3 operator()(Vector2 const& loc, PixelT alt) const {
      if (is_transparent(alt))
        return Vector3(0,0,std::numeric_limits<double>::quiet_NaN());

      Vector3 result;
      subvector(result, 0, 2) = m_georef.pixel_to_lonlat(loc);
      result.z() = alt;

      return result;
    }
  };

  /// Functor to convert lon/lat/alt to GCC x/y/z using a datum.
  class GeodeticToCartesian: public ReturnFixedType<Vector3> {
    Datum m_datum;
  public:
    GeodeticToCartesian(Datum const& d): m_datum(d) {}

    Vector3 operator()(Vector3 const& v) const;
  };

  /// Functor to convert GCC x/y/z to lon/lat/alt using a datum.
  class CartesianToGeodetic: public ReturnFixedType<Vector3> {
    Datum m_datum;
  public:
    CartesianToGeodetic(Datum const& d): m_datum(d) {}

    Vector3 operator()(Vector3 const& v) const;
  };

  /// Functor to convert lon/lat/alt to projected pixel col/row/alt using a georef.
  class GeodeticToProjection: public ReturnFixedType<Vector3> {
    // Use an alias, which is apparently thread-safe, as a copy is very expensive
    // for dynamic CRS.
    GeoReference const& m_reference;
  public:
    GeodeticToProjection(GeoReference const& r): m_reference(r) {}

    Vector3 operator()(Vector3 const& v) const;
  };

  /// Functor to convert projected pixel col/row/alt to lon/lat/alt using a georef.
  class ProjectionToGeodetic: public ReturnFixedType<Vector3> {
    // Use an alias, which is apparently thread-safe, as a copy is very expensive
    // for dynamic CRS.
    GeoReference const& m_reference;
  public:
    ProjectionToGeodetic(GeoReference const& r): m_reference(r) {}

    Vector3 operator()(Vector3 const& v) const;
  };

  /// Functor to convert lon/lat/alt to projected x/y/alt using a georef.
  class GeodeticToPoint: public ReturnFixedType<Vector3> {
    // Use an alias, which is apparently thread-safe, as a copy is very expensive
    // for dynamic CRS.
    GeoReference const& m_reference;
  public:
    GeodeticToPoint(GeoReference const& r): m_reference(r) {}

    Vector3 operator()(Vector3 const& v) const;
  };

  /// Functor to convert projected x/y/alt to lon/lat/alt using a georef.
  class PointToGeodetic: public ReturnFixedType<Vector3> {
    // Use an alias, which is apparently thread-safe, as a copy is very expensive
    // for dynamic CRS.
    GeoReference const& m_reference;
  public:
    PointToGeodetic(GeoReference const& r): m_reference(r) {}

    Vector3 operator()(Vector3 const& v) const;
  };

  /// Functor to convert col/row/alt values to projected point and alt
  template <class PixelT>
  class DemToProj: public ReturnFixedType<Vector3> {
    // Use an alias, which is apparently thread-safe, as a copy is very expensive
    // for dynamic CRS.
    GeoReference const& m_georef;
  public:
    DemToProj(GeoReference const& georef): m_georef(georef) {}

    Vector3 operator()(Vector2 const& loc, PixelT alt) const {
      if (is_transparent(alt))
        return Vector3(0,0,std::numeric_limits<double>::quiet_NaN());

      Vector3 result;
      subvector(result, 0, 2) = m_georef.pixel_to_point(loc);
      result.z() = alt;

      return result;
    }
  };

  // Image View operations


  template <class ImageT>
  BinaryPerPixelView<PerPixelIndexView<VectorIndexFunctor>, ImageT, DemToGeodetic<typename ImageT::pixel_type>>
  inline dem_to_geodetic(ImageViewBase<ImageT> const& dem, GeoReference const& georef) {
    typedef DemToGeodetic<typename ImageT::pixel_type> func_type;
    typedef BinaryPerPixelView<PerPixelIndexView<VectorIndexFunctor>, ImageT, func_type> result_type;
    func_type func(georef);
    return result_type(pixel_index_view(dem), dem.impl(), func);
  }

  template <class ImageT>
  UnaryPerPixelView<ImageT,GeodeticToCartesian>
  inline geodetic_to_cartesian(ImageViewBase<ImageT> const& lla_image, Datum const& d) {
    typedef UnaryPerPixelView<ImageT,GeodeticToCartesian> result_type;
    return result_type(lla_image.impl(), GeodeticToCartesian(d));
  }

  template <class ImageT>
  UnaryPerPixelView<ImageT,CartesianToGeodetic>
  inline cartesian_to_geodetic(ImageViewBase<ImageT> const& xyz_image, Datum const& d) {
    typedef UnaryPerPixelView<ImageT,CartesianToGeodetic> result_type;
    return result_type(xyz_image.impl(), CartesianToGeodetic(d));
  }

  template <class ImageT>
  UnaryPerPixelView<ImageT,GeodeticToCartesian>
  inline geodetic_to_cartesian(ImageViewBase<ImageT> const& lla_image, GeoReference const& r) {
    typedef UnaryPerPixelView<ImageT,GeodeticToCartesian> result_type;
    return result_type(lla_image.impl(), GeodeticToCartesian(r.datum()));
  }

  template <class ImageT>
  UnaryPerPixelView<ImageT,CartesianToGeodetic>
  inline cartesian_to_geodetic(ImageViewBase<ImageT> const& xyz_image, GeoReference const& r) {
    typedef UnaryPerPixelView<ImageT,CartesianToGeodetic> result_type;
    return result_type(xyz_image.impl(), CartesianToGeodetic(r.datum()));
  }

  template <class ImageT>
  UnaryPerPixelView<ImageT,GeodeticToProjection>
  inline geodetic_to_projection(ImageViewBase<ImageT> const& lla_image, GeoReference const& r) {
    typedef UnaryPerPixelView<ImageT,GeodeticToProjection> result_type;
    return result_type(lla_image.impl(), GeodeticToProjection(r));
  }

  template <class ImageT>
  UnaryPerPixelView<ImageT,ProjectionToGeodetic>
  inline projection_to_geodetic(ImageViewBase<ImageT> const& prj_image, GeoReference const& r) {
    typedef UnaryPerPixelView<ImageT,ProjectionToGeodetic> result_type;
    return result_type(prj_image.impl(), ProjectionToGeodetic(r));
  }

  /// Function to convert col/row/alt values to projected point and alt
  template <class ImageT>
  BinaryPerPixelView<PerPixelIndexView<VectorIndexFunctor>, ImageT, DemToProj<typename ImageT::pixel_type>>
  inline dem_to_proj(ImageViewBase<ImageT> const& dem, GeoReference const& georef) {
    typedef DemToProj<typename ImageT::pixel_type> func_type;
    typedef BinaryPerPixelView<PerPixelIndexView<VectorIndexFunctor>, ImageT, func_type> result_type;
    func_type func(georef);
    return result_type(pixel_index_view(dem), dem.impl(), func);
  }

  // Point is the intermediate projection step that is in units of the
  // projection not pixels. The only difference between projection and
  // point is an affine transform.
  template <class ImageT>
  UnaryPerPixelView<ImageT,GeodeticToPoint>
  inline geodetic_to_point(ImageViewBase<ImageT> const& lla_image, GeoReference const& r) {
    typedef UnaryPerPixelView<ImageT,GeodeticToPoint> result_type;
    return result_type(lla_image.impl(), GeodeticToPoint(r));
  }

  template <class ImageT>
  UnaryPerPixelView<ImageT,PointToGeodetic>
  inline point_to_geodetic(ImageViewBase<ImageT> const& point_image, GeoReference const& r) {
    typedef UnaryPerPixelView<ImageT,PointToGeodetic> result_type;
    return result_type(point_image.impl(), PointToGeodetic(r));
  }

  // XYZ to LON LAT ALT CONVERSION
  // WARNING: These functions are estimations, they do not produce accurate results.
  // They seem to assume the datum to be spherical.
  /// GCC to GDC conversion with elevation being distance from 0,0,0
  class XYZtoLonLatRadEstimateFunctor: public UnaryReturnSameType {
    bool m_east_positive;
    bool m_centered_on_zero; // Use the range [-180,180] otherwise [0,360]
  public:
    XYZtoLonLatRadEstimateFunctor(bool east_positive = true, bool centered_on_zero = true): m_east_positive(east_positive), m_centered_on_zero(centered_on_zero) {}

    template <class T>
    T operator()(T const& p) const {
      return this->apply(p, m_east_positive, m_centered_on_zero);
    }

    template <class T>
    static inline T apply(T const& p, bool east_positive = true,
                          bool centered_on_zero = true) {

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
      if (centered_on_zero) {
        if (lon > M_PI)
          lon -= 2*M_PI;
        if (lon < -M_PI)
          lon += 2*M_PI;
      } else {
        if (lon < 0)
          lon += 2*M_PI;
        if (lon > 2*M_PI)
          lon -= 2*M_PI;
      }

      return T (lon * 180.0 / M_PI, lat * 180.0 / M_PI, radius);
    }
  };

  /// GDC to GCC conversion with elevation being distance from 0,0,0
  class LonLatRadToXYZEstimateFunctor: public UnaryReturnSameType {
    bool m_east_positive;
  public:
    LonLatRadToXYZEstimateFunctor(bool east_positive = true): m_east_positive(east_positive) {}

    // Convert from lon, lat, radius to x,y,z:
    //
    // East positive:
    // x = r * cos(latitude) * cos(longitude)
    // y = r * cos(latitude) * sin(longitude)
    // z = r * sin(latitude)
    //
    // West positive:
    // x = r * cos(latitude) * cos(-longitude)
    // y = r * cos(latitude) * sin(-longitude)
    // z = r * sin(latitude)

    template <class T>
    T operator()(T const& p) const {
      return this->apply(p, m_east_positive);
    }

    template <class T>
    static inline T apply(T const& p, bool east_positive = true) {
      typename T::value_type z = p(2) * sin(p(1)*M_PI/180.0);
      typename T::value_type sqrt_x_sqr_plus_y_sqr = p(2) * cos(p(1)*M_PI/180.0);

      if (east_positive) {
        return Vector3(sqrt_x_sqr_plus_y_sqr * cos(p(0)*M_PI/180.0),
                        sqrt_x_sqr_plus_y_sqr * sin(p(0)*M_PI/180.0),
                        z);
      } else {
        return Vector3(sqrt_x_sqr_plus_y_sqr * cos(-p(0)*M_PI/180.0),
                        sqrt_x_sqr_plus_y_sqr * sin(-p(0)*M_PI/180.0),
                        z);
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
  UnaryPerPixelView<ImageT, XYZtoLonLatRadEstimateFunctor>
  inline xyz_to_lon_lat_radius_estimate(ImageViewBase<ImageT> const& image,
                                         bool east_positive    = true,
                                         bool centered_on_zero = true) {
    return UnaryPerPixelView<ImageT,XYZtoLonLatRadEstimateFunctor>(image.impl(), XYZtoLonLatRadEstimateFunctor(east_positive, centered_on_zero));
  }

  template <class ElemT>
  inline Vector<ElemT,3> xyz_to_lon_lat_radius_estimate(Vector<ElemT,3> const& xyz,
                                                          bool east_positive = true,
                                                          bool centered_on_zero = true) {
    return XYZtoLonLatRadEstimateFunctor::apply(xyz, east_positive, centered_on_zero);
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
  UnaryPerPixelView<ImageT, LonLatRadToXYZEstimateFunctor>
  inline lon_lat_radius_to_xyz_estimate(ImageViewBase<ImageT> const& image, bool east_positive = true) {
    return UnaryPerPixelView<ImageT,LonLatRadToXYZEstimateFunctor>(image.impl(), LonLatRadToXYZEstimateFunctor(east_positive));
  }

  template <class ElemT>
  inline Vector<ElemT,3> lon_lat_radius_to_xyz_estimate(Vector<ElemT,3> const& lon_lat_alt, bool east_positive = true) {
    return LonLatRadToXYZEstimateFunctor::apply(lon_lat_alt, east_positive);
  }

}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_POINTIMAGEMANIPLULATION_H__
