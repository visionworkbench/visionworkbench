// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Cartography/detail/BresenhamLine.h>
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
    BBox2 lonlat_bbox;

    // Testing the parameter of the pixel bbox
    for ( int32 x=pixel_bbox.min().x(); x<pixel_bbox.max().x(); ++x ) {
      try {
        lonlat_bbox.grow(pixel_to_lonlat( Vector2(x,pixel_bbox.min().y()) ));
        lonlat_bbox.grow(pixel_to_lonlat( Vector2(x,pixel_bbox.max().y()-1) ));
      } catch ( const cartography::ProjectionErr& e ) {}
    }
    for ( int32 y=pixel_bbox.min().y()+1; y<pixel_bbox.max().y()-1; ++y ) {
      try {
        lonlat_bbox.grow(pixel_to_lonlat( Vector2(pixel_bbox.min().x(),y) ));
        lonlat_bbox.grow(pixel_to_lonlat( Vector2(pixel_bbox.max().x()-1,y) ));
      } catch ( const cartography::ProjectionErr& e ) {}
    }

    // Drawing an X inside the bbox. This covers the poles. It will
    // produce a lonlat boundary that is within at least one pixel of
    // the pole. This will also help catch terminator boundaries from
    // orthographic projections.
    BresenhamLine l1( pixel_bbox.min(), pixel_bbox.max() );
    while ( l1.is_good() ) {
      try {
        lonlat_bbox.grow( pixel_to_lonlat( *l1 ) );
      } catch ( const cartography::ProjectionErr& e ) {}
      ++l1;
    }
    BresenhamLine l2( pixel_bbox.min() + Vector2i(pixel_bbox.width(),0),
                      pixel_bbox.max() - Vector2i(pixel_bbox.width(),0) );
    while ( l2.is_good() ) {
      try {
        lonlat_bbox.grow( pixel_to_lonlat( *l2 ) );
      } catch ( const cartography::ProjectionErr& e ) {}
      ++l2;
    }

    return lonlat_bbox;
  }

  BBox2i GeoReferenceBase::lonlat_to_pixel_bbox(BBox2 const& lonlat_bbox, size_t nsamples) const {
    // Alternatively this function could avoid the nsamples
    // option. The sample discrete step could just be this average
    // size of pixel in degrees.

    BBox2 pixel_bbox;

    Vector2 lower_fraction(lonlat_bbox.width()/double(nsamples),
                           lonlat_bbox.height()/double(nsamples));
    for(size_t i = 0; i < nsamples; i++) {
      try {
        // Walk the top & bottom (technically past the edge of pixel space) rows
        double x = lonlat_bbox.min().x() + double(i) * lower_fraction.x();
        pixel_bbox.grow(lonlat_to_pixel(Vector2(x,lonlat_bbox.min().y())));
        pixel_bbox.grow(lonlat_to_pixel(Vector2(x,lonlat_bbox.max().y())));


        // Walk the left & right (technically past the edge of pixel space) columns
        double y = lonlat_bbox.min().y() + double(i) * lower_fraction.y();
        pixel_bbox.grow(lonlat_to_pixel(Vector2(lonlat_bbox.min().x(),y)));
        pixel_bbox.grow(lonlat_to_pixel(Vector2(lonlat_bbox.max().x(),y)));
      } catch ( const cartography::ProjectionErr& e ) {}
    }

    // It is possible that this may not required. However in the
    // cartography it seems better to be rigorous than sorry.
    BresenhamLine l1( Vector2i(), Vector2i(nsamples,nsamples) );
    while ( l1.is_good() ) {
      try {
        pixel_bbox.grow( lonlat_to_pixel( elem_prod(Vector2(*l1),lower_fraction) + lonlat_bbox.min() ) );
      } catch ( const cartography::ProjectionErr& e ) {}
      ++l1;
    }
    BresenhamLine l2( Vector2i(nsamples,0), Vector2i(0,nsamples) );
    while ( l2.is_good() ) {
      try {
        pixel_bbox.grow( lonlat_to_pixel( elem_prod(Vector2(*l2),lower_fraction) + lonlat_bbox.min() ) );
      } catch ( const cartography::ProjectionErr& e ) {}
      ++l2;
    }

    return grow_bbox_to_int(pixel_bbox);
  }

}} // namespace vw, cartography
