// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Cartography/GeoReferenceBase.h>
#include <vw/Image/Transform.h>

namespace vw {
namespace cartography {

  /// For a bbox in projected space, return the corresponding bbox in
  /// pixels on the image
  BBox2i GeoReferenceBase::point_to_pixel_bbox(BBox2 const& point_bbox) const {
    BBox2 pixel_bbox;
    pixel_bbox.grow(point_to_pixel(point_bbox.min()));
    pixel_bbox.grow(point_to_pixel(point_bbox.max()));
    pixel_bbox.grow(point_to_pixel(Vector2(point_bbox.min().x(), point_bbox.max().y())));
    pixel_bbox.grow(point_to_pixel(Vector2(point_bbox.max().x(), point_bbox.min().y())));
    return grow_bbox_to_int(pixel_bbox);
  }

  BBox2 GeoReferenceBase::pixel_to_lonlat_bbox(BBox2i const& pixel_bbox) const {
    // TODO: This should be tested with all the different projections

    BBox2 lonlat_bbox;

    // As all the projections are continuous, we can just walk the
    // edges to find the lonlat bounding box.

    // We take 10 samples per edge
    const int nsamples = 10;

    for(int i = 0; i < nsamples; i++) {
      // Walk the top & bottom (technically past the edge of pixel space) rows
      double x = pixel_bbox.min().x() + double(i) / nsamples * pixel_bbox.width();
      lonlat_bbox.grow(pixel_to_lonlat(Vector2(x,pixel_bbox.min().y())));
      lonlat_bbox.grow(pixel_to_lonlat(Vector2(x,pixel_bbox.max().y())));

      // Walk the left & right (technically past the edge of pixel space) columns
      double y = pixel_bbox.min().y() + double(i) / nsamples * pixel_bbox.height();
      lonlat_bbox.grow(pixel_to_lonlat(Vector2(pixel_bbox.min().x(),y)));
      lonlat_bbox.grow(pixel_to_lonlat(Vector2(pixel_bbox.max().x(),y)));
    }

    // Do we cross the north or south pole? Have to cover that case
    // specially. Fortunately it's easy, because (anything, 90) or
    // (anything, -90) will always be in the image.

    // TODO: Should this be done with Bresham lines through the
    // image instead?

    // North pole:
    if (pixel_bbox.contains(lonlat_to_pixel(Vector2(0, 90)))) {
        lonlat_bbox.grow(Vector2(0, 90));
    }

    // South pole:
    if (pixel_bbox.contains(lonlat_to_pixel(Vector2(0, -90)))) {
        lonlat_bbox.grow(Vector2(0, -90));
    }

    return lonlat_bbox;
  }

  BBox2i GeoReferenceBase::lonlat_to_pixel_bbox(BBox2 const& lonlat_bbox) const {
    // TODO: This should be tested with all the different projections

    BBox2 pixel_bbox;

    // As all the projections are continuous, we can just walk the
    // edges to find the lonlat bounding box.

    // We take 10 samples per edge
    const int nsamples = 10;

    for(int i = 0; i < nsamples; i++) {
      // Walk the top & bottom (technically past the edge of pixel space) rows
      double x = lonlat_bbox.min().x() + double(i) / nsamples * lonlat_bbox.width();
      pixel_bbox.grow(lonlat_to_pixel(Vector2(x,lonlat_bbox.min().y())));
      pixel_bbox.grow(lonlat_to_pixel(Vector2(x,lonlat_bbox.max().y())));


      // Walk the left & right (technically past the edge of pixel space) columns
      double y = lonlat_bbox.min().y() + double(i) / nsamples * lonlat_bbox.height();
      pixel_bbox.grow(lonlat_to_pixel(Vector2(lonlat_bbox.min().x(),y)));
      pixel_bbox.grow(lonlat_to_pixel(Vector2(lonlat_bbox.max().x(),y)));
    }

    // TODO: Do poles need to be taken into account?

    return grow_bbox_to_int(pixel_bbox);
  }

}} // namespace vw, cartography
