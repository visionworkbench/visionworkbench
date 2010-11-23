// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Cartography/GeoTransform.h>

// Vision Workbench
#include <vw/Image/ImageView.h>

// Proj.4
#include <projects.h>

namespace {
  class BresenhamLine {
    vw::int32 x0, y0, x1, y1;
    vw::int32 x, y;
    bool steep;
    vw::int32 deltax, deltay, error, ystep;
  public:
    BresenhamLine( vw::Vector2i const& start, vw::Vector2i const& stop ) :
    x0(start[0]), y0(start[1]), x1(stop[0]), y1(stop[1]) {
      steep = abs(y1-y0) > abs(x1-x0);
      if (steep) {
        std::swap(x0,y0);
        std::swap(x1,y1);
      }
      if ( x0 > x1 ) {
        std::swap(x0,x1);
        std::swap(y0,y1);
      }
      deltax = x1 - x0;
      deltay = abs(y1-y0);
      error = deltax / 2;
      ystep = y0 < y1 ? 1 : -1;
      x = x0; y = y0;
    }

    vw::Vector2i operator*() const {
      if (steep)
        return vw::Vector2i(y,x);
      else
        return vw::Vector2i(x,y);
    }

    void operator++() {
      x++;
      error -= deltay;
      if ( error < 0 ) {
        y += ystep;
        error += deltax;
      }
    }

    bool is_good() const { return x < x1; }
  };
}


namespace vw {
namespace cartography {

  // Constructor
    GeoTransform::GeoTransform(GeoReference const& src_georef, GeoReference const& dst_georef) :
    m_src_georef(src_georef), m_dst_georef(dst_georef) {
    const std::string src_datum = m_src_georef.datum().proj4_str();
    const std::string dst_datum = m_dst_georef.datum().proj4_str();

    // This optimizes in the common case where the two images are
    // already in the same map projection, and we need only apply
    // the affine transform.  This will break, of course, as soon as
    // we have mare than one type of GeoReference object, but it
    // makes life faster for now. -mbroxton
    if (m_src_georef.proj4_str() == m_dst_georef.proj4_str())
      m_skip_map_projection = true;
    else
      m_skip_map_projection = false;

    // This optimizes the case where the two datums are the same,
    // and thus we don't need to call proj to convert between them
    // as we transform.
    if(src_datum == dst_datum) {
      m_skip_datum_conversion = true;
    } else {
      // Set up the various variables for proj.
      m_skip_datum_conversion = false;

      // The source proj4 context.
      std::stringstream ss_src;
      // We convert lat/long to lat/long regardless of what the
      // source or destination georef uses.
      ss_src << "+proj=latlong " << src_datum;
      m_src_datum = boost::shared_ptr<ProjContext>(new ProjContext(ss_src.str()));
      CHECK_PROJ_INIT_ERROR(ss_src.str().c_str());

      // The destination proj4 context.
      std::stringstream ss_dst;
      ss_dst << "+proj=latlong " << dst_datum;
      m_dst_datum = boost::shared_ptr<ProjContext>(new ProjContext(ss_dst.str()));
      CHECK_PROJ_INIT_ERROR(ss_dst.str().c_str());
    }
    // Because GeoTransform is typically very slow, we default to a tolerance
    // of 0.1 pixels to allow ourselves to be approximated.
    set_tolerance( 0.1 );
  }

  // Performs a forward or reverse datum conversion.
  Vector2 GeoTransform::datum_convert(Vector2 const& v, bool forward) const {
    double x = v[0];
    double y = v[1];
    double z = 0;

    if(forward)
      pj_transform(m_src_datum->proj_ptr(), m_dst_datum->proj_ptr(), 1, 0, &x, &y, &z);
    else
      pj_transform(m_dst_datum->proj_ptr(), m_src_datum->proj_ptr(), 1, 0, &x, &y, &z);
    CHECK_PROJ_ERROR;

    return Vector2(x, y);
  }

  BBox2i GeoTransform::forward_bbox( BBox2i const& bbox ) const {
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

  BBox2i GeoTransform::reverse_bbox( BBox2i const& bbox ) const {
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

  void reproject_point_image(ImageView<Vector3> const& point_image,
                             GeoReference const& src_georef,
                             GeoReference const& dst_georef) {

    GeoTransform gtx(src_georef, dst_georef);

    // Iterate over the image, transforming the first two coordinates
    // in the Vector one at a time.  The third coordinate is taken to
    // be the altitude value, and this value is not touched.
    for (int32 j=0; j < point_image.rows(); ++j) {
      for (int32 i=0; i < point_image.cols(); ++i) {
        if (point_image(i,j) != Vector3()) {
          Vector2 in(point_image(i,j)[0], point_image(i,j)[1]);
          Vector2 out = gtx.forward(in);
          point_image(i,j).x() = out[0];
          point_image(i,j).y() = out[1];
        }
      }
    }
  }


}} // namespace vw::cartography

