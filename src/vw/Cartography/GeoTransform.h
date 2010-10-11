// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_CARTOGRAPHY_GEOTRANSFORM_H__
#define __VW_CARTOGRAPHY_GEOTRANSFORM_H__

#include <sstream>
#include <string>

#include <vw/Math/Vector.h>
#include <vw/Image/Transform.h>
#include <vw/Cartography/GeoReference.h>

namespace vw {
namespace cartography {

  class GeoTransform : public TransformHelper<GeoTransform,ContinuousFunction,ContinuousFunction> {

    GeoReference m_src_georef;
    GeoReference m_dst_georef;
    boost::shared_ptr<ProjContext> m_src_datum;
    boost::shared_ptr<ProjContext> m_dst_datum;
    bool m_skip_map_projection;
    bool m_skip_datum_conversion;

    /* Converts between datums. The parameter 'forward' specifies whether
     * we convert forward (true) or reverse (false).
    */
    Vector2 datum_convert(Vector2 const& v, bool forward) const;

  public:
    /// Normal constructor
    GeoTransform(GeoReference const& src_georef, GeoReference const& dst_georef);

    /// Given a pixel coordinate of an image in a destination
    /// georeference frame, this routine computes the corresponding
    /// pixel from an image in the source georeference frame.
    Vector2 reverse(Vector2 const& v) const {
      if (m_skip_map_projection)
        return m_src_georef.point_to_pixel(m_dst_georef.pixel_to_point(v));
      if(m_skip_datum_conversion)
        return m_src_georef.lonlat_to_pixel(m_dst_georef.pixel_to_lonlat(v));
      Vector2 dst_lonlat = m_dst_georef.pixel_to_lonlat(v);
      dst_lonlat = datum_convert(dst_lonlat, false);
      return m_src_georef.lonlat_to_pixel(dst_lonlat);
    }

    /// Given a pixel coordinate of an image in a source
    /// georeference frame, this routine computes the corresponding
    /// pixel the destination (transformed) image.
    Vector2 forward(Vector2 const& v) const {
      if (m_skip_map_projection)
        return m_dst_georef.point_to_pixel(m_src_georef.pixel_to_point(v));
      if(m_skip_datum_conversion)
        return m_dst_georef.lonlat_to_pixel(m_src_georef.pixel_to_lonlat(v));
      Vector2 src_lonlat = m_src_georef.pixel_to_lonlat(v);
      src_lonlat = datum_convert(src_lonlat, true);
      return m_dst_georef.lonlat_to_pixel(src_lonlat);
    }

    // We override forward_bbox so it understands to check if the image
    // crosses the poles or not.
    BBox2i forward_bbox( BBox2i const& bbox ) const;

    // We do the same for reverse_bbox
    BBox2i reverse_bbox( BBox2i const& bbox ) const;
  };


  // ---------------------------------------------------------------------------
  // Functional API
  // ---------------------------------------------------------------------------

  /// Returns a transformed image view.  The user can specify the type
  /// of interpolation and edge extension to be done by supplying the
  /// appropriate functors in the last two arguments.  For example:
  /// See the transform() function in Transform.h for more details.
  template <class ImageT, class EdgeT, class InterpT>
  typename boost::disable_if<IsScalar<InterpT>, TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, GeoTransform> >::type
  inline geo_transform( ImageViewBase<ImageT> const& v,
                    GeoReference const& src_georef,
                    GeoReference const& dst_georef,
                    EdgeT const& edge_func,
                    InterpT const& interp_func ) {
    return TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, GeoTransform>
      (interpolate(v, interp_func, edge_func), GeoTransform(src_georef,dst_georef));
  }

  /// Convenience function: transform with a default interpolation
  /// scheme of bilinear interpolation.
  template <class ImageT, class EdgeT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, GeoTransform>
  inline geo_transform( ImageViewBase<ImageT> const& v,
                    GeoReference const& src_georef,
                    GeoReference const& dst_georef,
                    EdgeT const& edge_func ) {
    return TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, GeoTransform>
      (interpolate(v, BilinearInterpolation(), edge_func), GeoTransform(src_georef,dst_georef));
  }

  /// Convenience function: transform with a default scheme of
  /// bilinear interpolation and zero edge extension.
  template <class ImageT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, ZeroEdgeExtension>, BilinearInterpolation>, GeoTransform>
  inline geo_transform( ImageViewBase<ImageT> const& v,
                    GeoReference const& src_georef,
                    GeoReference const& dst_georef ) {
    return TransformView<InterpolationView<EdgeExtensionView<ImageT, ZeroEdgeExtension>, BilinearInterpolation>, GeoTransform>
      (interpolate(v, BilinearInterpolation(), ZeroEdgeExtension()), GeoTransform(src_georef,dst_georef));
  }


  /// This variant of transform allows the user to specify the
  /// dimensions of the transformed image.  The upper left hand point
  /// (0,0) stays fixed.  For a more flexible method of cropping to an
  /// arbitrary bounding box, use one of the transform methods defined
  /// below.
  template <class ImageT, class EdgeT, class InterpT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, GeoTransform>
  inline geo_transform( ImageViewBase<ImageT> const& v,
                    GeoReference const& src_georef,
                    GeoReference const& dst_georef,
                    int32 width,
                    int32 height,
                    EdgeT const& edge_func,
                    InterpT const& interp_func ) {
    return TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, GeoTransform>
      (interpolate(v, interp_func, edge_func), GeoTransform(src_georef,dst_georef), width, height);
  }

  /// Convenience function: transform with a default interpolation
  /// scheme of bilinear interpolation. The user can specify the
  /// dimensions of the output image.
  template <class ImageT, class EdgeT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, GeoTransform>
  inline geo_transform( ImageViewBase<ImageT> const& v,
                    GeoReference const& src_georef,
                    GeoReference const& dst_georef,
                    int32 width,
                    int32 height,
                    EdgeT const& edge_func ) {
    return TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, GeoTransform>
      (interpolate(v, BilinearInterpolation(), edge_func), GeoTransform(src_georef,dst_georef), width, height);
  }

  /// Convenience function: transform with a default scheme of
  /// bilinear interpolation and zero edge extension.  The user can
  /// specify the dimensions of the output image.
  template <class ImageT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, ZeroEdgeExtension>, BilinearInterpolation>, GeoTransform>
  inline geo_transform( ImageViewBase<ImageT> const& v,
                    GeoReference const& src_georef,
                    GeoReference const& dst_georef,
                    int32 width,
                    int32 height ) {
    return TransformView<InterpolationView<EdgeExtensionView<ImageT, ZeroEdgeExtension>, BilinearInterpolation>, GeoTransform>
      (interpolate(v, BilinearInterpolation(), ZeroEdgeExtension()), GeoTransform(src_georef,dst_georef), width, height);
  }


  /// Reproject an image whose pixels contain 3D points (usually in
  /// some spherical coordinate system).  Important note: it is
  /// assumed here that the 3D points already have the affine
  /// transform applied to them (they correspond to real 3D
  /// coordinates and not pixel coordinates in an image), therefore
  /// the affine transform portion of the georeference is completely
  /// ignored by the function.  It does not matter what affine
  /// transform you are using in the src_georef or dst_georef.
  ///
  /// Important Note: The convention here is that the Vector3 contains
  /// the ordered triple: (longitude, latitude, altitude).
  void reproject_point_image(ImageView<Vector3> const& point_image,
                             GeoReference const& src_georef,
                             GeoReference const& dst_georef);

}} // namespace vw::cartography

#endif // __GEO_TRANSFORM_H__
