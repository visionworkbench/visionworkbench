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


#ifndef __VW_CARTOGRAPHY_POINTIMAGEMANIPLULATION_H__
#define __VW_CARTOGRAPHY_POINTIMAGEMANIPLULATION_H__

#include <vw/Math/Vector.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/UtilityViews.h>
#include <vw/Cartography/GeoReference.h>

// This include is here to keep compat (the contents of that header used to be
// here, and was split up to break the Camera<=>Cartography circular dep).
#include <vw/Cartography/SimplePointImageManipulation.h>

/// \file PointImageManipulation.h
///
/// Contains routines for manipulating ImageViews with a pixel type of
/// Vector3 that contain xyz or lon, lat, altitude triples.

namespace vw {
namespace cartography {

  // Functors. View operations are lower in this file.
  template <class PixelT>
  class DemToGeodetic : public ReturnFixedType<Vector3> {
    GeoReference m_georef;
  public:
    DemToGeodetic(GeoReference const& georef) : m_georef(georef) {}

    Vector3 operator()(Vector2 const& loc, PixelT alt) const {
      if (is_transparent(alt))
        return Vector3(0,0,std::numeric_limits<double>::quiet_NaN());

      Vector3 result;
      subvector(result, 0, 2) = m_georef.pixel_to_lonlat(loc);
      result.z() = alt;

      return result;
    }
  };

  class GeodeticToCartesian : public ReturnFixedType<Vector3> {
    Datum m_datum;
  public:
    GeodeticToCartesian(Datum const& d) : m_datum(d) {}

    Vector3 operator()( Vector3 const& v ) const;
  };

  class CartesianToGeodetic : public ReturnFixedType<Vector3> {
    Datum m_datum;
  public:
    CartesianToGeodetic(Datum const& d) : m_datum(d) {}

    Vector3 operator()( Vector3 const& v ) const;
  };

  class GeodeticToProjection : public ReturnFixedType<Vector3> {
    GeoReference m_reference;
  public:
    GeodeticToProjection( GeoReference const& r ) : m_reference(r) {}

    Vector3 operator()( Vector3 const& v ) const;
  };

  class ProjectionToGeodetic : public ReturnFixedType<Vector3> {
    GeoReference m_reference;
  public:
    ProjectionToGeodetic( GeoReference const& r ) : m_reference(r) {}

    Vector3 operator()( Vector3 const& v ) const;
  };

  class GeodeticToPoint : public ReturnFixedType<Vector3> {
    GeoReference m_reference;
  public:
    GeodeticToPoint( GeoReference const& r ) : m_reference(r) {}

    Vector3 operator()( Vector3 const& v ) const;
  };

  class PointToGeodetic : public ReturnFixedType<Vector3> {
    GeoReference m_reference;
  public:
    PointToGeodetic( GeoReference const& r ) : m_reference(r) {}

    Vector3 operator()( Vector3 const& v ) const;
  };

  // Image View operations
  template <class ImageT>
  BinaryPerPixelView<PerPixelIndexView<VectorIndexFunctor>, ImageT, DemToGeodetic<typename ImageT::pixel_type> >
  inline dem_to_geodetic(ImageViewBase<ImageT> const& dem, GeoReference const& georef) {
    typedef DemToGeodetic<typename ImageT::pixel_type> func_type;
    typedef BinaryPerPixelView<PerPixelIndexView<VectorIndexFunctor>, ImageT, func_type> result_type;
    func_type func(georef);
    return result_type(pixel_index_view(dem), dem.impl(), func);
  }

  template <class ImageT>
  UnaryPerPixelView<ImageT,GeodeticToCartesian>
  inline geodetic_to_cartesian( ImageViewBase<ImageT> const& lla_image, Datum const& d ) {
    typedef UnaryPerPixelView<ImageT,GeodeticToCartesian> result_type;
    return result_type(lla_image.impl(), GeodeticToCartesian(d) );
  }

  template <class ImageT>
  UnaryPerPixelView<ImageT,CartesianToGeodetic>
  inline cartesian_to_geodetic( ImageViewBase<ImageT> const& xyz_image, Datum const& d ) {
    typedef UnaryPerPixelView<ImageT,CartesianToGeodetic> result_type;
    return result_type(xyz_image.impl(), CartesianToGeodetic(d) );
  }

  template <class ImageT>
  UnaryPerPixelView<ImageT,GeodeticToCartesian>
  inline geodetic_to_cartesian( ImageViewBase<ImageT> const& lla_image, GeoReference const& r ) {
    typedef UnaryPerPixelView<ImageT,GeodeticToCartesian> result_type;
    return result_type(lla_image.impl(), GeodeticToCartesian(r.datum()) );
  }

  template <class ImageT>
  UnaryPerPixelView<ImageT,CartesianToGeodetic>
  inline cartesian_to_geodetic( ImageViewBase<ImageT> const& xyz_image, GeoReference const& r ) {
    typedef UnaryPerPixelView<ImageT,CartesianToGeodetic> result_type;
    return result_type(xyz_image.impl(), CartesianToGeodetic(r.datum()) );
  }

  template <class ImageT>
  UnaryPerPixelView<ImageT,GeodeticToProjection>
  inline geodetic_to_projection( ImageViewBase<ImageT> const& lla_image, GeoReference const& r ) {
    typedef UnaryPerPixelView<ImageT,GeodeticToProjection> result_type;
    return result_type(lla_image.impl(), GeodeticToProjection(r) );
  }

  template <class ImageT>
  UnaryPerPixelView<ImageT,ProjectionToGeodetic>
  inline projection_to_geodetic( ImageViewBase<ImageT> const& prj_image, GeoReference const& r ) {
    typedef UnaryPerPixelView<ImageT,ProjectionToGeodetic> result_type;
    return result_type(prj_image.impl(), ProjectionToGeodetic(r));
  }

  // Point is the intermediate projection step that is in units of the
  // projection not pixels. The only difference between projection and
  // point is an affine transform.
  template <class ImageT>
  UnaryPerPixelView<ImageT,GeodeticToPoint>
  inline geodetic_to_point( ImageViewBase<ImageT> const& lla_image, GeoReference const& r ) {
    typedef UnaryPerPixelView<ImageT,GeodeticToPoint> result_type;
    return result_type(lla_image.impl(), GeodeticToPoint(r) );
  }

  template <class ImageT>
  UnaryPerPixelView<ImageT,PointToGeodetic>
  inline point_to_geodetic( ImageViewBase<ImageT> const& point_image, GeoReference const& r ) {
    typedef UnaryPerPixelView<ImageT,PointToGeodetic> result_type;
    return result_type(point_image.impl(), PointToGeodetic(r));
  }
}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_POINTIMAGEMANIPLULATION_H__

