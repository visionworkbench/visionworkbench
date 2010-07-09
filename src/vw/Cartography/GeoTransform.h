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

    class BresenhamLine {
      int32 x0, y0, x1, y1;
      int32 x, y;
      bool steep;
      int32 deltax, deltay, error, ystep;
    public:
      BresenhamLine( Vector2i const& start, Vector2i const& stop ) :
      x0(start[0]), y0(start[1]), x1(stop[0]), y1(stop[1]), x(start[0]), y(start[1]) {
        steep = abs(y1-y0) > abs(x1-x0);
        if (steep) {
          std::swap(x0,y0);
          std::swap(x1,y1);
        }
        deltax = x1 - x0;
        deltay = abs(y1-y0);
        error = deltax / 2;
        ystep = y0 < y1 ? 1 : -1;
      }

      Vector2i operator*() {
        if (steep)
          return Vector2i(y,x);
        else
          return Vector2i(x,y);
      }

      void operator++() {
        x++;
        error -= deltay;
        if ( error < 0 ) {
          y += ystep;
          error += deltax;
        }
      }

      bool is_good() { return x != x1 && y != y1; }
    };

    // We override forward_bbox so it understands to check if the image 
    // crosses the poles or not.
    BBox2i forward_bbox( BBox2i const& bbox ) const {
      BBox2 r = TransformHelper<GeoTransform,ContinuousFunction,ContinuousFunction>::forward_bbox(bbox);
      BresenhamLine l1( bbox.min(), bbox.max() );
      while ( l1.is_good() ) {
        try {
          r.grow( this->forward( *l1 ) );
        } catch ( cartography::ProjectionErr const& e ) {}
        ++l1;
      }
      BresenhamLine l2( bbox.min() + Vector2i(bbox.width(),0),
                        bbox.max() + Vector2i(-bbox.width(),0) );
      while ( l2.is_good() ) {
        try {
          r.grow( this->forward( *l2 ) );
        } catch ( cartography::ProjectionErr const& e ) {}
        ++l2;
      }

      return grow_bbox_to_int(r);
    }

    // We do the same for reverse_bbox
    BBox2i reverse_bbox( BBox2i const& bbox ) const {
      BBox2 r = TransformHelper<GeoTransform,ContinuousFunction,ContinuousFunction>::reverse_bbox(bbox);
      BresenhamLine l1( bbox.min(), bbox.max() );
      while ( l1.is_good() ) {
        try {
          r.grow( this->reverse( *l1 ) );
        } catch ( cartography::ProjectionErr const& e ) {}
        ++l1;
      }
      BresenhamLine l2( bbox.min() + Vector2i(bbox.width(),0),
                        bbox.max() + Vector2i(-bbox.width(),0) );
      while ( l2.is_good() ) {
        try {
          r.grow( this->reverse( *l2 ) );
        } catch ( cartography::ProjectionErr const& e ) {}
        ++l2;
      }

      return grow_bbox_to_int(r);
    }

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
