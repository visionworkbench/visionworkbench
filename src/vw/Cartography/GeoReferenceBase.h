// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_CARTOGRAPHY_GEOREFERENCE_H__
#define __VW_CARTOGRAPHY_GEOREFERENCE_H__

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/Algorithms.h>
#include <vw/Cartography/Datum.h>

namespace vw {
namespace cartography {

  /// The georeference class contains the mapping from image
  /// coordinates (u,v) to geospatial coordinates (typically lat/lon,
  /// or possibly meters in a UTM grid cell, etc.).  It must also
  /// encode how to translate between this coordinate system and the
  /// "Geographic" coordinate system (lat,lon)
  class GeoReferenceBase {
  public:
    /// The affine transform converts from pixel space to geographic
    /// or projected space and vice versa.  Most often, this process
    /// entails interpolating based on floating point pixel
    /// coordinates in the image.  However, images are discrete
    /// samples of pixel space, so you must adopt a convention
    /// regarding how floating point pixel coordinates in your
    /// georeferenced image are to be interpreted.
    ///
    /// You have one of two choices: If you assume PixelAsArea, the
    /// upper left hand corner of the top left pixel is considered as
    /// the origin (0,0), and the center of the top left pixel is
    /// (0.5, 0.5).  This assumption is common when dealing with
    /// satellite imagery or maps.
    ///
    /// On the other hand, if you assume the PixelAsPoint, then the
    /// center of the upper left hand pixel is the origin (0,0), and
    /// the top left corner of that pixel is at (-0.5,-0.5) in pixel
    /// coordinates.  This mode is common when working with elevation
    /// data, etc.
    ///
    /// Note: The Vision Workbench *always* interprets floating point
    /// pixel location (0,0) as being at the _center_ of the upper
    /// left hand pixel.  If you choose the PixelAsArea option for
    /// this flag, the GeoTransform class will automatically adjust
    /// your affine transform my (0.5,0.5) to bring the coordinate
    /// system in line with the Vision Workbench internal
    /// representation.
    ///
    /// The default pixel interpretation for GeoReference is PixelAsArea
    enum PixelInterpretation { PixelAsArea = 0, PixelAsPoint = 1 };

  protected:
    PixelInterpretation m_pixel_interpretation;
    Datum m_datum;

  public:
    PixelInterpretation pixel_interpretation() const { return m_pixel_interpretation; }
    void set_pixel_interpretation(PixelInterpretation const& p) { m_pixel_interpretation = p; }

    /// Default Constructor
    GeoReferenceBase() : m_pixel_interpretation( GeoReferenceBase::PixelAsArea ) {}

    /// Takes a geodetic datum.
    GeoReferenceBase(Datum const& datum) : m_pixel_interpretation( GeoReferenceBase::PixelAsArea ), m_datum(datum) {}

    /// Destructor.
    virtual ~GeoReferenceBase() {}

    Datum const& datum() const { return m_datum; }
    void set_datum(Datum const& datum) { m_datum = datum; }

    /// For a given pixel coordinate, compute the position of that
    /// pixel in this georeferenced space.
    virtual Vector2 pixel_to_point(Vector2 pix) const = 0;

    /// For a given location 'loc' in projected space, compute the
    /// corresponding pixel coordinates in the image.
    virtual Vector2 point_to_pixel(Vector2 loc) const = 0;


    /// For a point in the projected space, compute the position of
    /// that point in unprojected (Geographic) coordinates (lat,lon).
    virtual Vector2 point_to_lonlat(Vector2 loc) const = 0;

    /// Given a position in geographic coordinates (lat,lon), compute
    /// the location in the projected coordinate system.
    virtual Vector2 lonlat_to_point(Vector2 lat_lon) const = 0;


    /// For a given pixel coordinate, compute the position of that
    /// pixel in Geographic coordinates (lat,lon).
    virtual Vector2 pixel_to_lonlat(Vector2 pix) const {
      return point_to_lonlat(pixel_to_point(pix));
    }

    /// Given a position in geographic coordinates (lat,lon), compute
    /// the location in pixel coordinates in this image that
    /// corresponds to the given geographic coordinates.
    virtual Vector2 lonlat_to_pixel(Vector2 lat_lon) const {
      return point_to_pixel(lonlat_to_point(lat_lon));
    }

    /// For a given pixel bbox, return the corresponding bbox in projected
    /// space
    virtual BBox2 pixel_to_point_bbox(BBox2i pixel_bbox) const {
      BBox2 point_bbox;
      point_bbox.grow(pixel_to_point(pixel_bbox.min()));
      point_bbox.grow(pixel_to_point(pixel_bbox.max()));
      point_bbox.grow(pixel_to_point(Vector2(pixel_bbox.min().x(), pixel_bbox.max().y())));
      point_bbox.grow(pixel_to_point(Vector2(pixel_bbox.max().x(), pixel_bbox.min().y())));
      return point_bbox;
    }

    /// For a bbox in projected space, return the corresponding bbox in
    /// pixels on the image
    virtual BBox2i point_to_pixel_bbox(BBox2 const& point_bbox) const;

    /// For a given pixel bbox, return the corresponding bbox in
    /// Geographic (lon, lat) coordinates
    virtual BBox2 pixel_to_lonlat_bbox(BBox2i const& pixel_bbox) const;

    /// For a given bbox in Geographic coordinates, return the corresponding
    /// bbox in image pixel coordinates
    virtual BBox2i lonlat_to_pixel_bbox(BBox2 const& lonlat_bbox) const;

    // TODO: Probably should implement point_to_lonlat_bbox and
    // lonlat_to_point_bbox for completeness sake, but right now they're
    // not needed for anything

    /// Return the box that bounds the area represented by the
    /// geotransform for the dimensions of the given image.
    template <class ViewT>
    BBox2 bounding_box(ImageViewBase<ViewT> const& view) const VW_DEPRECATED;

    /// Return the box that bounds the area represented by the
    /// geotransform for the dimensions of the given image.
    /// Note that this doesn't tell you whether the image takes the
    /// long path or the short path from the left longitude to the
    /// right longitude.
    ///
    /// Assumption: that the projection is continuous.
    template <class ViewT>
    BBox2 lonlat_bounding_box(ImageViewBase<ViewT> const& view) const VW_DEPRECATED;

  };

  template <class ViewT>
  BBox2 GeoReferenceBase::bounding_box(ImageViewBase<ViewT> const& view) const {
    return pixel_to_point_bbox(vw::bounding_box(view.impl()));
  }

  template <class ViewT>
  BBox2 GeoReferenceBase::lonlat_bounding_box(ImageViewBase<ViewT> const& view) const {
    return pixel_to_lonlat_bbox(vw::bounding_box(view.impl()));
  }

}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_GEOREFERENCE_H__
