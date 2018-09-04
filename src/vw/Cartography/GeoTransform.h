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


#ifndef __VW_CARTOGRAPHY_GEOTRANSFORM_H__
#define __VW_CARTOGRAPHY_GEOTRANSFORM_H__

#include <sstream>
#include <string>

#include <vw/Core/Thread.h>
#include <vw/Math/Vector.h>
#include <vw/Image/Transform.h>
#include <vw/Cartography/GeoReference.h>

/// \file GeoTransform.h Tools for converting points and images between geographic datums.

namespace vw {
namespace cartography {


  /// A structure to allow the fast and accurate conversion of coordinates between
  /// two GeoReference objects.
  /// - Use of this class is the ONLY safe method to convert coordinates between two
  ///   GeoReference objects that fully handles all variables.
  /// - Be very careful using this class to convert between datums!  Not all datums will correctly convert.
  ///   To be safe verify that your datums convert properly or use a dedicated tool such as 
  ///   ASP's datum_convert to do the conversions.
  class GeoTransform : public TransformHelper<GeoTransform,ContinuousFunction,ContinuousFunction> {

    GeoReference  m_src_georef;
    GeoReference  m_dst_georef;
    BBox2         m_src_bbox,
                  m_dst_bbox;
    ProjContext   m_src_datum_proj,
                  m_dst_datum_proj;
    bool          m_skip_map_projection;
    bool          m_skip_datum_conversion;
    mutable Mutex m_mutex; // Used to control access to the ProjContext objects

  public:
  
    /// Default constructor, does not generate a usable object.
    GeoTransform() {}

    GeoTransform(GeoTransform const& other);

    /// Normal constructor
    GeoTransform(GeoReference const& src_georef, GeoReference const& dst_georef,
                 BBox2 const& src_bbox = BBox2i(0, 0, 0, 0),
                 BBox2 const& dst_bbox = BBox2i(0, 0, 0, 0));

    GeoTransform& operator=(GeoTransform const& other);

    //---------------------------------------------------------------
    // These functions implement the Transform interface and allow
    //  this class to be passed in to functions expecting a Transform object.

    /// Given a pixel coordinate of an image in a source
    /// georeference frame, this routine computes the corresponding
    /// pixel in the destination (transformed) image.
    Vector2 forward(Vector2 const& v) const {return pixel_to_pixel(v);}

    /// Given a pixel coordinate of an image in a destination
    /// georeference frame, this routine computes the corresponding
    /// pixel from an image in the source georeference frame.
    Vector2 reverse(Vector2 const& v) const;

    /// Convert a pixel bounding box in the source image to
    ///  a pixel bounding box in the destination image.
    /// - This function handles the case where the image crosses the poles.
    BBox2i forward_bbox( BBox2i const& bbox ) const;

    /// As forward_bbox, but from destination to source.
    BBox2i reverse_bbox( BBox2i const& bbox ) const;


    //------------------------------------------------------------
    // These functions do not implement the Transform interface.

    /// Convert a pixel in the source to a pixel in the destination.
    Vector2 pixel_to_pixel( Vector2 const& v ) const;

    /// Convert a point in the source to a point (projected coords) in the destination.
    Vector2 point_to_point( Vector2 const& v ) const;

    /// Convert a pixel in the source to a point (projected coords) in the destination.
    Vector2 pixel_to_point( Vector2 const& v ) const;

    /// Convert a point in the source to a pixel in the destination.
    Vector2 point_to_pixel( Vector2 const& v ) const;

    /// Converts lonlat coords, taking the datums into account.
    /// - The parameter 'forward' specifies whether we convert forward (true) or reverse (false).
    Vector2 lonlat_to_lonlat(Vector2 const& lonlat, bool forward=true) const;

    /// Converts lonlatalt coords, taking the datums into account.
    /// - The parameter 'forward' specifies whether we convert forward (true) or reverse (false).
    Vector3 lonlatalt_to_lonlatalt(Vector3 const& lonlatalt, bool forward=true) const;

    /// Returns true if bounding box conversions wrap around the output
    ///  georeference, creating a very large bounding box.
    bool check_bbox_wraparound() const;


    /// Convert a lonlat bounding box in the source to a lonlat in the destination.
    BBox2 lonlat_to_lonlat_bbox( BBox2 const& pixel_bbox ) const;

    /// Convert a pixel bounding box in the source to a point (projected coords)
    ///  bounding box in the destination.
    BBox2 pixel_to_pixel_bbox( BBox2 const& pixel_bbox ) const {return(forward_bbox(pixel_bbox));}

    /// Convert a pixel bounding box in the source to a point (projected coords)
    ///  bounding box in the destination.
    BBox2 pixel_to_point_bbox( BBox2 const& pixel_bbox ) const;

    /// Convert a point bounding box in the source to a pixel bounding box in the destination.
    BBox2 point_to_pixel_bbox( BBox2 const& point_bbox ) const;

    
    friend std::ostream& operator<<(std::ostream& os, const GeoTransform& trans);
  }; // End class GeoTransform

  /// Format a GeoTransform to a text stream (for debugging)
  std::ostream& operator<<(std::ostream& os, const GeoTransform& trans);


  // ---------------------------------------------------------------------------
  // Image View Functions
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
  /// arbitrary bounding box, use one of the transform methods defined below.
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


  // ---------------------------------------------------------------------------
  // Miscellaneous Functions
  // ---------------------------------------------------------------------------


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
