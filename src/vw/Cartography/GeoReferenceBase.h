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
#ifndef __VW_CARTOGRAPHY_GEOREFERENCE_H__
#define __VW_CARTOGRAPHY_GEOREFERENCE_H__

#include <vw/Core/Exception.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/PerPixelViews.h>
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
    /// The default pixel interpretation for GeoReference is PixelAsPoint
    enum PixelInterpretation { PixelAsArea, PixelAsPoint };

  protected:
    PixelInterpretation m_pixel_interpretation;
    Datum m_datum;

  public:
    PixelInterpretation pixel_interpretation() const { return m_pixel_interpretation; }
    void set_pixel_interpretation(PixelInterpretation const& p) { m_pixel_interpretation = p; }

    /// Default Constructor
    GeoReferenceBase() : m_pixel_interpretation( GeoReferenceBase::PixelAsPoint ) {}

    /// Takes a geodetic datum.
    GeoReferenceBase(Datum const& datum) : m_pixel_interpretation( GeoReferenceBase::PixelAsPoint ), m_datum(datum) {}

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

    /// Return the box that bounds the area represented by the
    /// geotransform for the dimensions of the given image.
    template <class ViewT> 
    BBox2 bounding_box(ImageViewBase<ViewT> const& view) const {
      BBox2 bbox;
      bbox.grow(pixel_to_point(Vector2(0,0)));
      bbox.grow(pixel_to_point(Vector2(view.impl().cols(),0)));
      bbox.grow(pixel_to_point(Vector2(0,view.impl().rows())));
      bbox.grow(pixel_to_point(Vector2(view.impl().cols(), view.impl().rows())));
      return bbox;
    }

    /// Return the box that bounds the area represented by the
    /// geotransform for the dimensions of the given image.
    template <class ViewT> 
    BBox2 lonlat_bounding_box(ImageViewBase<ViewT> const& view) const {
      BBox2 bbox;
      bbox.grow(pixel_to_lonlat(Vector2(0,0)));
      bbox.grow(pixel_to_lonlat(Vector2(view.impl().cols(),0)));
      bbox.grow(pixel_to_lonlat(Vector2(0,view.impl().rows())));
      bbox.grow(pixel_to_lonlat(Vector2(view.impl().cols(), view.impl().rows())));
      return bbox;
    }
  };

}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_GEOREFERENCE_H__
