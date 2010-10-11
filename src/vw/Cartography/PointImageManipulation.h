// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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

  // --------------------- CHANGE OF PROJECTION ----------------------

  /// Move from a source projection to a destination projection.  The
  /// altitude value (that is, the third value in the vector) is left
  /// untouched.
  class ReprojectPointFunctor : public UnaryReturnSameType {
    GeoReference m_src, m_dst;

  public:
    ReprojectPointFunctor(GeoReference const& src, GeoReference const& dst) : m_src(src), m_dst(dst) {}

    template <class T>
    T operator()(T const& p) const {
      Vector<typename T::value_type,2> lon_lat = m_dst.lonlat_to_point(m_src.point_to_lonlat(subvector(p,0,2)));
      return T(lon_lat(0), lon_lat(1),
               p(2)+m_src.datum().radius(lon_lat(0), lon_lat(1))-m_dst.datum().radius(lon_lat(0), lon_lat(1)));
    }
  };

  // This version of the functor assumes that the point inputs are
  // unprojected (lon, lat, radius) and then projects them when forward = true.
  // When forward = false, the point inputs are assumed to be projected, and then
  // unprojects them.
  class ProjectPointFunctor : public UnaryReturnSameType {
    GeoReference m_dst;
    bool m_forward;

  public:
    ProjectPointFunctor(GeoReference const& dst, bool forward = true) : m_dst(dst), m_forward(forward) {}

    template <class T>
    T operator()(T const& p) const {
      if (p == T()) return p;
      Vector<typename T::value_type,2> lon_lat = m_dst.lonlat_to_point(subvector(p,0,2));
      typename T::value_type offset = m_dst.datum().radius(lon_lat(0), lon_lat(1));
      return T(lon_lat(0), lon_lat(1),
               m_forward ? p(2) - offset : p(2) + offset);
    }
  };

  template <class PixelT>
  class DemToPointImageFunctor : public ReturnFixedType<Vector3> {
      GeoReference m_georef;
    public:
      DemToPointImageFunctor(GeoReference georef) : m_georef(georef) {}

      Vector3 operator()(Vector2 loc, PixelT alt) const {
        if (is_transparent(alt))
          return Vector3();

        Vector3 result;
        subvector(result, 0, 2) = m_georef.pixel_to_lonlat(loc);
        result.z() = alt;

        return result;
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
  UnaryPerPixelView<ImageT, ReprojectPointFunctor>
  inline reproject_point_image( ImageViewBase<ImageT> const& image, GeoReference const& src_georef, GeoReference const& dst_georef) {
    return UnaryPerPixelView<ImageT,ReprojectPointFunctor>( image.impl(), ReprojectPointFunctor(src_georef, dst_georef) );
  }

  // This variant, which only accepts a destination projection,
  // assumes that the source points are [lon, lat, radius] values.
  template <class ImageT>
  UnaryPerPixelView<ImageT, ProjectPointFunctor>
  inline project_point_image( ImageViewBase<ImageT> const& image, GeoReference const& dst_georef, bool forward=true) {
    return UnaryPerPixelView<ImageT,ProjectPointFunctor>( image.impl(), ProjectPointFunctor(dst_georef, forward) );
  }

  // This utility function converts a DEM to a point image
  template <class ImageT>
  BinaryPerPixelView<PerPixelIndexView<VectorIndexFunctor>, ImageT, DemToPointImageFunctor<typename ImageT::pixel_type> >
  inline dem_to_point_image(ImageViewBase<ImageT> const& dem, GeoReference georef) {
    typedef DemToPointImageFunctor<typename ImageT::pixel_type> func_type;
    typedef BinaryPerPixelView<PerPixelIndexView<VectorIndexFunctor>, ImageT, func_type> result_type;
    func_type func(georef);
    return result_type(pixel_index_view(dem), dem.impl(), func);
  }
}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_POINTIMAGEMANIPLULATION_H__

