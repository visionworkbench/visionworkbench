// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Cartography/ToastTransform.h>


// A helper function to convert a point on the unit sphere to
// a lon/lat vector.
vw::Vector2 vw::cartography::ToastTransform::unitvec_to_lonlat(vw::Vector3 const& vec) const {
  double lat = (180.0/M_PI) * asin(vec.z());
  if( vec.x()==0 && vec.y()==0 ) {
    return Vector2(0, lat);
  }
  double lon = (180.0/M_PI) * atan2(vec.y(), vec.x());
  return Vector2(lon, lat);
}


// A helper function to convert a lon/lat vector to a point on the
// unit sphere.
vw::Vector3 vw::cartography::ToastTransform::lonlat_to_unitvec(vw::Vector2 const& lonlat) const {
  double lon = (M_PI/180)*lonlat.x(), lat = (M_PI/180)*lonlat.y();
  double sinlon = sin(lon), coslon = cos(lon);
  double sinlat = sin(lat), coslat = cos(lat);
  return Vector3( coslon*coslat, sinlon*coslat, sinlat );
}


// Maps the unit right triangle onto the first octant of the unit
// sphere.  In particular, maps (0,0) onto (0,0,1), (1,0) onto
// (1,0,0), and (0,1) onto (0,1,0).
//
// Operates by tesselating one octant of an icosohedron
// iteratively down to roughly the level of one pixel
// at the requested resolution; then linearly interpolates
// within the terminal triangle.
vw::Vector3 vw::cartography::ToastTransform::octant_point_to_unitvec(double x, double y) const {
  Vector3 c1(0,0,1), c2(1,0,0), c3(0,1,0);
  double epsilon = 1.0/m_resolution;
  while( epsilon < 1.0 ) {
    if( x < 0.5 ) {
      if( y < 0.5 ) {
        if( y < 0.5 - x ) {
          x = 2*x;
          y = 2*y;
          c2 = normalize(c1 + c2);
          c3 = normalize(c3 + c1);
        }
        else {
          x = 1-2*x;
          y = 1-2*y;
          Vector3 c12 = normalize(c1 + c2);
          Vector3 c31 = normalize(c3 + c1);
          c1 = normalize(c2 + c3);
          c2 = c31;
          c3 = c12;
        }
      }
      else {
        x = 2*x;
        y = 2*y-1;
        c2 = normalize(c2 + c3);
        c1 = normalize(c3 + c1);
      }
    }
    else {
      x = 2*x-1;
      y = 2*y;
      c1 = normalize(c1 + c2);
      c3 = normalize(c2 + c3);
    }
    epsilon *= 2;
  }
  return normalize(c1 + x*(c2-c1) + y*(c3-c1));
}


// Maps the first octant of the unit sphere onto the unit right
// triangle.  In particular, maps (0,0,1) onto (0,0), (1,0,0) onto
// (1,0), and (0,1,0) onto (0,1).
//
// Operates much like octant_point_to_unitvec().
vw::Vector2 vw::cartography::ToastTransform::octant_unitvec_to_point(vw::Vector3 const& vec) const {
  Vector3 c1(0,0,1), c2(1,0,0), c3(0,1,0);
  Vector2 p1(0,0), p2(1,0), p3(0,1);
  double epsilon = 1.0/m_resolution;
  while( epsilon < 1.0 ) {
    Vector3 c12 = normalize(c1+c2);
    Vector3 c13 = normalize(c1+c3);
    Vector3 c23 = normalize(c2+c3);
    if( dot_prod(cross_prod(c12,c23),vec) > 0 ) {
      if( dot_prod(cross_prod(c23,c13),vec) > 0 ) {
        if( dot_prod(cross_prod(c12,c13),vec) > 0 ) {
          c2 = c12;
          c3 = c13;
          p2 = (p1+p2)/2;
          p3 = (p1+p3)/2;
        }
        else {
          c1 = c23;
          c2 = c13;
          c3 = c12;
          Vector2 p12 = (p1+p2)/2;
          Vector2 p13 = (p1+p3)/2;
          p1 = (p2+p3)/2;
          p2 = p13;
          p3 = p12;
        }
      }
      else {
        c1 = c13;
        c2 = c23;
        p1 = (p1+p3)/2;
        p2 = (p2+p3)/2;
      }
    }
    else {
      c1 = c12;
      c3 = c23;
      p1 = (p1+p2)/2;
      p3 = (p2+p3)/2;
    }
    epsilon *= 2;
  }

  // We assume that vec now lies essentially in the plane defined by
  // the corner points c1,c2,c3.  We construct a local coordinate
  // system using those points and use it to interpolate the corner
  // point values.  In principle the result we'd get using an SVD
  // might be slightly better, since it would minimize the
  // square-error associated with the assumption of coplanarity, but
  // this is pretty good and avoids extra dependencies.
  Matrix3x3 M;
  select_col(M,0) = c2-c1;
  select_col(M,1) = c3-c1;
  select_col(M,2) = cross_prod(c2-c1,c3-c1);
  Vector3 coeffs = inverse(M) * (vec-c1);
  return p1 + (p2-p1)*coeffs.x() + (p3-p1)*coeffs.y();
}


// Forward-projects a pixel location in the projected source image
// space into a pixel location in the TOAST image space.
vw::Vector2 vw::cartography::ToastTransform::forward(vw::Vector2 const& point) const {
  // There is a fundamental eight-fold symmetry to the TOAST
  // projection which we exploit here.  We first determine which
  // top-level triangle (i.e. which octant) the requested point lies
  // within.

  // Compte sanitize lon/lat coordinates
  Vector2 lonlat = m_georef.pixel_to_lonlat(point);
  double lon = lonlat.x();
  if( lon > 180 ) lon -= 360;
  if( lon <= -180 ) lon += 360;
  double lat = lonlat.y();
  if( lat > 90 ) lat = 90;
  if( lat < -90 ) lat = -90;

  if( lon<90 && lon>-90 ) {
    if( lon >= 0 ) {
      // Lower left: 0 to 90E
      if( lat < 0 ) {
        Vector2 octpt = octant_lonlat_to_point(90-lonlat.x(), -lonlat.y());
        return Vector2( octpt.x()/2, 1.0-octpt.y()/2 ) * (m_resolution-1);
      }
      else {
        Vector2 octpt = octant_lonlat_to_point(lonlat.x(), lonlat.y());
        return Vector2( (1.0-octpt.x())/2, 1.0-(1.0-octpt.y())/2 ) * (m_resolution-1);
      }
    }
    else {
      // Upper left: 9 to 90W
      if( lat < 0 ) {
        Vector2 octpt = octant_lonlat_to_point(90+lonlat.x(), -lonlat.y());
        return Vector2( octpt.x()/2, 1.0-(2.0-octpt.y())/2 ) * (m_resolution-1);
      }
      else {
        Vector2 octpt = octant_lonlat_to_point(-lonlat.x(), lonlat.y());
        return Vector2( (1-octpt.x())/2, 1.0-(octpt.y()+1.0)/2 ) * (m_resolution-1);
      }
    }
  }
  else {
    if( lon >= 0 ) {
      // Lower right: 90E to 180
      if( lat < 0 ) {
        Vector2 octpt = octant_lonlat_to_point(lonlat.x()-90, -lonlat.y());
        return Vector2( (2.0-octpt.x())/2, 1.0-octpt.y()/2 ) * (m_resolution-1);
      }
      else {
        Vector2 octpt = octant_lonlat_to_point(180-lonlat.x(), lonlat.y());
        return Vector2( (octpt.x()+1)/2, 1.0-(1.0-octpt.y())/2 ) * (m_resolution-1);
      }
    }
    else {
      // Upper right: 90W to 180
      if( lat < 0 ) {
        Vector2 octpt = octant_lonlat_to_point(-90-lonlat.x(), -lonlat.y());
        return Vector2( (2.0-octpt.x())/2, 1.0-(2.0-octpt.y())/2 ) * (m_resolution-1);
      }
      else {
        Vector2 octpt = octant_lonlat_to_point(180+lonlat.x(), lonlat.y());
        return Vector2( (octpt.x()+1.0)/2, 1.0-(octpt.y()+1.0)/2 ) * (m_resolution-1);
      }
    }
  }
}


// Back-projects a pixel location in the TOAST image space into a
// pixel location in the projected source image space.
vw::Vector2 vw::cartography::ToastTransform::reverse(vw::Vector2 const& point) const {
  // There is a fundamental eight-fold symmetry to the TOAST
  // projection which we exploit here.  We first determine which
  // top-level triangle (i.e. which octant) the requested point lies
  // within.  We then map that octant onto a unit triangle, call
  // octant_point_to_lonlat to project onto the sphere, and then
  // rework the resulting lat/lon back into the proper octant.
  double x = point.x()/(m_resolution-1);
  double y = 1.0-point.y()/(m_resolution-1);
  if( x < 0.5 ) {
    if( y < 0.5 ) {
      // Lower left: 0 to 90E
      if( y < 0.5 - x ) {
        Vector2 lonlat = octant_point_to_lonlat(2*x, 2*y);
        return m_georef.lonlat_to_pixel(Vector2(90-lonlat.x(), -lonlat.y()));
      }
      else {
        Vector2 lonlat = octant_point_to_lonlat(1-2*x, 1-2*y);
        return m_georef.lonlat_to_pixel(Vector2(lonlat.x(), lonlat.y()));
      }
    }
    else {
      // Upper left: 0 to 90W
      if( y > 0.5 + x ) {
        Vector2 lonlat = octant_point_to_lonlat(2*x, 2-2*y);
        return m_georef.lonlat_to_pixel(Vector2(-90+lonlat.x(), -lonlat.y()));
      }
      else {
        Vector2 lonlat = octant_point_to_lonlat(1-2*x,2*y-1);
        return m_georef.lonlat_to_pixel(Vector2(-lonlat.x(), lonlat.y()));
      }
    }
  }
  else {
    // Lower right: 90E to 180
    if( y < 0.5 ) {
      if( y < x - 0.5 ) {
        Vector2 lonlat = octant_point_to_lonlat(2-2*x, 2*y);
        return m_georef.lonlat_to_pixel(Vector2(90+lonlat.x(), -lonlat.y()));
      }
      else {
        Vector2 lonlat = octant_point_to_lonlat(2*x-1, 1-2*y);
        return m_georef.lonlat_to_pixel(Vector2(180-lonlat.x(), lonlat.y()));
      }
    }
    else {
      // Upper right: 90W to 180
      if( y > 1.5 - x ) {
        Vector2 lonlat = octant_point_to_lonlat(2-2*x, 2-2*y);
        return m_georef.lonlat_to_pixel(Vector2(-90-lonlat.x(), -lonlat.y()));
      }
      else {
        Vector2 lonlat = octant_point_to_lonlat(2*x-1, 2*y-1);
        return m_georef.lonlat_to_pixel(Vector2(-180+lonlat.x(), lonlat.y()));
      }
    }
  }
}


// We override forward_bbox so it understands to check if the image
// crosses the poles or not.
vw::BBox2i vw::cartography::ToastTransform::forward_bbox( vw::BBox2i const& bbox ) const {
  // If the source bounding box contains the south pole, then the dest
  // bounding box is the entire TOAST projection space, since the
  // south pole is mapped to the four corners of TOAST.

  Vector2 south_pole_pixel;
  try {
    south_pole_pixel = m_georef.lonlat_to_pixel(Vector2i(0,-90));
  } catch (const ProjectionErr& e) {
    // We asked for a point not defined in the projection, most likely.  Assume
    // the point represents a notch discontinuity, and ask for a cross through it.
    // The center of that cross will hopefully be the correct point.
    Vector2 a,b,c,d;

    a = m_georef.lonlat_to_pixel(Vector2i(  0, -89.9));
    b = m_georef.lonlat_to_pixel(Vector2i( 90, -89.9));
    c = m_georef.lonlat_to_pixel(Vector2i(180, -89.9));
    d = m_georef.lonlat_to_pixel(Vector2i(270, -89.9));

    south_pole_pixel = (a + b + c + d)/4.;
  }

  if( bbox.contains(south_pole_pixel) )
    return BBox2i(0,0,m_resolution,m_resolution);

  BBox2 src_bbox = TransformHelper<ToastTransform,ContinuousFunction,ContinuousFunction>::forward_bbox(bbox);
  return grow_bbox_to_int(src_bbox);
}


// We override reverse_bbox so it understands to check if the image
// crosses the poles or not.
vw::BBox2i vw::cartography::ToastTransform::reverse_bbox( vw::BBox2i const& bbox, bool approximate ) const {

  BBox2 src_bbox;
  if (approximate) {
    src_bbox.grow( reverse( Vector2(bbox.min().x(),bbox.min().y()) ) ); // Top left
    src_bbox.grow( reverse( Vector2(bbox.max().x()-1,bbox.min().y()) ) ); // Top right
    src_bbox.grow( reverse( Vector2(bbox.min().x(),bbox.max().y()-1) ) ); // Bottom left
    src_bbox.grow( reverse( Vector2(bbox.max().x()-1,bbox.max().y()-1) ) ); // Bottom right
    grow_bbox_to_int( src_bbox );
  } else {
    src_bbox = TransformHelper<ToastTransform,ContinuousFunction,ContinuousFunction>::reverse_bbox(bbox);
  }
  reverse_bbox_poles( bbox, src_bbox );
  return grow_bbox_to_int(src_bbox);
}


// Attempt to expand the given bounding box in the source pixel space
// to include any poles containd in the given bounding box in the
// destiation (TOAST) pixel space.
void vw::cartography::ToastTransform::reverse_bbox_poles(vw::BBox2 const& dst_bbox, vw::BBox2 & src_bbox) const {
  if( dst_bbox.contains( Vector2(m_resolution/2, m_resolution/2) ) ) {
    src_bbox.grow( m_georef.lonlat_to_pixel(Vector2(-180,90)) );
    src_bbox.grow( m_georef.lonlat_to_pixel(Vector2(0,90)) );
    src_bbox.grow( m_georef.lonlat_to_pixel(Vector2(180,90)) );
  }
  if( dst_bbox.contains( Vector2(0,0) ) ||
      dst_bbox.contains( Vector2(m_resolution-1,0) ) ||
      dst_bbox.contains( Vector2(0,m_resolution-1) ) ||
      dst_bbox.contains( Vector2(m_resolution-1, m_resolution-1) ) ) {
    src_bbox.grow( m_georef.lonlat_to_pixel(Vector2(-180,-90)) );
    src_bbox.grow( m_georef.lonlat_to_pixel(Vector2(0,-90)) );
    src_bbox.grow( m_georef.lonlat_to_pixel(Vector2(180,-90)) );
  }
}


// A heuristic that back-projects a line on to a conservative
// bounding box.  Used by SparseImageCheck.
//
// Sample the line connecting the two points in TOAST image space
// using the given number of intervals.  For each interval, project
// the endpoints and midpoint back into the source image space,
// noting the deviation of the midpoint from a linear estimate.
// Returns the bounding box of all the reverse-projected points,
// padded by the greatest estimation error encountered.
vw::BBox2 vw::cartography::ToastTransform::reverse_line( Vector2 const& a, Vector2 const& b, int num_divisions ) const {
  Vector2 last_point = reverse( a );
  BBox2 bbox;
  Vector2 max_error;
  for( int i=0; i<num_divisions; ++i ) {
    Vector2 mid_point = reverse( a + (i+0.5)*(b-a)/num_divisions );
    bbox.grow( mid_point );
    Vector2 next_point = reverse( a + (i+1.0)*(b-a)/num_divisions );
    bbox.grow( next_point );
    Vector2 error = mid_point - (last_point+next_point)/2;
    if( fabs(error.x()) > max_error.x() ) max_error.x() = fabs(error.x());
    if( fabs(error.y()) > max_error.y() ) max_error.y() = fabs(error.y());
    last_point = next_point;
  }
  bbox.min() -= max_error;
  bbox.max() += max_error;
  return bbox;
}
