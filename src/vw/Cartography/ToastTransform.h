// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file ToastTransform.h
///
/// A transform class that transforms images from a standard
/// cartographic projection into a TOAST projection.
///
#ifndef __VW_CARTOGRAPHY_TOASTTRANSFORM_H__
#define __VW_CARTOGRAPHY_TOASTTRANSFORM_H__

// THE TOAST PROJECTION
//
// TOAST stands for Tessellated Octahedral Adaptive Subdivision
// Transform.  I wouldn't necessarily call the transform itself
// "adaptive", but it does make for a cute acronym.
//
// Each octant of the sphere is mapped onto one triangle of an
// octahedron, and this octahedron is then unfolded and mapped
// onto a square, like so:
//
//             90W
//    +---------+---------+
//    |       / | \       |
//    |     /   |   \     |
//    |   /     |     \   |
//    | /       |       \ |
// 0E +---------+---------+ 180E
//    | \       |       / |
//    |   \     |     /   |
//    |     \   |   /     |
//    |       \ | /       |
//    +---------+---------+
//             90E
//
// The north pole is at the center of the square, and the south pole
// appears at all four corners.  The topology is somewhat unusual:
// each edge of the square connects to itself with a 180-degree
// rotation.
//
// Within each octant, the mapping between the triangle and the sphere
// is based on iterative tesselation of the octahedron.  At each
// iteration each edge of each triangle is divided in half, and the
// midpoints are lifted up to the sphere.  This divides each triangle
// into four smaller triangles, and the new edges are likewise lifted
// to the sphere as great circles.  Lather, rinse, repeat.
//
// A more complete description is available here:
// http://research.microsoft.com/en-us/um/people/dinos/spheretoaster.pdf
//
// The TOAST transform is not cheap, and we work around that by using
// the lookup-table-based approximation capabilities of TransformView.
// We could get back some precision, and possibly make things even
// faster, by exploiting the fact that the first many iterations will
// be identical for nearby points.  We typically compute the reverse
// transform for a large number of nearby points at a time.  We could
// either cache the stack used in the most recent calculation, or
// directly expose a the ability to reverse-transform a block of
// points at a time.

#include <vw/Math/Vector.h>
#include <vw/Math/BBox.h>
#include <vw/Image/Transform.h>
#include <vw/Image/SparseImageCheck.h>
#include <vw/Cartography/GeoReference.h>

namespace vw {
namespace cartography {

  class ToastTransform : public TransformHelper<ToastTransform,ContinuousFunction,ContinuousFunction> {
    GeoReference m_georef;
    int32 m_resolution;

    // A helper function to convert a point on the unit sphere to
    // a lon/lat vector.
    Vector2 unitvec_to_lonlat(Vector3 const& vec) const;

    // A helper function to convert a lon/lat vector to a point on the
    // unit sphere.
    Vector3 lonlat_to_unitvec(Vector2 const& lonlat) const;

    // Maps the unit right triangle onto the first octant of the unit
    // sphere.
    Vector3 octant_point_to_unitvec(double x, double y) const;

    // Maps the first octant of the unit sphere onto the unit right
    // triangle.
    Vector2 octant_unitvec_to_point(Vector3 const& vec) const;

    // Convert a normalized point located in octant 0 to lon/lat in degrees
    inline Vector2 octant_point_to_lonlat(double x, double y) const {
      return unitvec_to_lonlat(octant_point_to_unitvec(x, y));
    }

    // Convert a lon/lat point located in octant 0 to a normalized point
    inline Vector2 octant_lonlat_to_point(double lon, double lat) const {
      return octant_unitvec_to_point(lonlat_to_unitvec(Vector2(lon,lat)));
    }

  public:
    ToastTransform(GeoReference const& georef, int32 resolution)
      : m_georef(georef), m_resolution(resolution)
    {
      // We enable approximation of the TOAST transform by
      // linear interpolation into lookup tables, because it is
      // quite slow.  The 0.2 pixel tolerance is fairly wide; it
      // would be great to speed up the algorithm enough to at
      // least use the 0.1 pixel tolerance that is used by default
      // for GeoTransform.
      set_tolerance(0.2);
    }

    virtual Vector2 forward( Vector2 const& point ) const;
    virtual Vector2 reverse( Vector2 const& point ) const;

    virtual BBox2i forward_bbox( BBox2i const& bbox ) const;

    // We override reverse_bbox so it understands to check if the image crosses
    // the poles or not.  Pass in 'approximate' to compute reverse bounding
    // boxes by transforming corner coordinates only. (This is much faster, but
    // not 100% accurate.)

    // We need to have one prototype that matches the Transform prototype or we
    // get a warning.
    virtual BBox2i reverse_bbox( BBox2i const& bbox) const {
      return this->reverse_bbox( bbox, false );
    }
    virtual BBox2i reverse_bbox( BBox2i const& bbox, bool approximate ) const;

    // Attempt to expand the given bounding box in the source pixel space
    // to include any poles containd in the given bounding box in the
    // destiation (TOAST) pixel space.  Public so it can easily be used
    // by SparseImageCheck.
    void reverse_bbox_poles( BBox2 const& dst_bbox, BBox2 & src_bbox ) const;

    // A heuristic that back-projects a line on to a conservative
    // bounding box.  Used by SparseImageCheck.
    BBox2 reverse_line( Vector2 const& a, Vector2 const& b, int num_divisions ) const;
  };

} // namespace vw::cartography

  template <class ChildT>
  class SparseImageCheck<TransformView<ChildT, cartography::ToastTransform> > {

    TransformView<ChildT, cartography::ToastTransform> m_view;

  public:
    SparseImageCheck(TransformView<ChildT, cartography::ToastTransform> const& source)
      : m_view(source) {}

    bool operator()( BBox2i const& bbox ) const {
      // We could just call ToastTransform::reverse_bbox() here, but we
      // anticipate getting called with very large bounding boxes, so
      // even walking around the edges could be unacceptably slow.
      // Since this function is just an optimization heuristic, we only
      // sample a small number of points around the border of the tile,
      // and then pad out the resulting bounding box aggressively based
      // on a measurement of the sampling error.  This will result in
      // some false positives, but it will be orders of magnitude faster
      // than reverse_bbox() for large tiles.

      cartography::ToastTransform const& txform = m_view.transform();

      BBox2 src_bbox;
      src_bbox.grow( txform.reverse_line( Vector2( bbox.min().x(), bbox.min().y() ), Vector2( bbox.max().x(), bbox.min().y() ), 8 ) );
      src_bbox.grow( txform.reverse_line( Vector2( bbox.min().x(), bbox.min().y() ), Vector2( bbox.min().x(), bbox.max().y() ), 8 ) );
      src_bbox.grow( txform.reverse_line( Vector2( bbox.max().x(), bbox.min().y() ), Vector2( bbox.max().x(), bbox.max().y() ), 8 ) );
      src_bbox.grow( txform.reverse_line( Vector2( bbox.min().x(), bbox.max().y() ), Vector2( bbox.max().x(), bbox.max().y() ), 8 ) );
      txform.reverse_bbox_poles( bbox, src_bbox );

      return SparseImageCheck<ChildT>(m_view.child())( src_bbox );
    }
  };

} // namespace vw

#endif // __VW_CARTOGRAPHY_TOASTTRANSFORM_H__
