// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Cartography/CameraBBox.h>

using namespace vw;

// Return map projected point location (the intermediate between LLA
// and Pixel)
Vector2
cartography::geospatial_intersect( Vector2 pix,
                                   cartography::GeoReference const& georef,
                                   boost::shared_ptr<camera::CameraModel> camera_model,
                                   double z_scale, bool& did_intersect ) {
  Vector3 ccenter = camera_model->camera_center( pix );
  Vector3 cpoint = camera_model->pixel_to_vector( pix );
  ccenter.z() *= z_scale;
  cpoint.z() *= z_scale;
  cpoint = normalize( cpoint );

  double radius_2 = georef.datum().semi_major_axis()*
    georef.datum().semi_major_axis();
  double alpha = - dot_prod(ccenter,cpoint);
  Vector3 projection = ccenter + alpha * cpoint;
  if( norm_2_sqr(projection) > radius_2 ) {  // the ray does not intersect the sphere
    did_intersect = false;
    return Vector2();
  } else {
    did_intersect = true;
  }

  alpha -= sqrt( radius_2 -
                 norm_2_sqr(projection) );
  Vector3 intersection = ccenter + alpha * cpoint;
  intersection.z() /= z_scale;

  Vector3 llr = georef.datum().cartesian_to_geodetic( intersection );
  Vector2 geospatial_point = georef.lonlat_to_point( Vector2( llr.x(),
                                                              llr.y() ) );
  return geospatial_point;
}

// Compute the bounding box in points (georeference space) that is
// defined by georef. Scale is MPP as georeference space is in meters.
BBox2 cartography::camera_bbox( cartography::GeoReference const& georef,
                                boost::shared_ptr<camera::CameraModel> camera_model,
                                int32 cols, int32 rows, float &scale ) {

  // Testing to see if we should be centering on zero
  bool center_on_zero = true;
  Vector3 camera_llr =
    XYZtoLonLatRadFunctor::apply(camera_model->camera_center(Vector2()));
  if ( camera_llr[0] < -90 ||
       camera_llr[0] > 90 )
    center_on_zero = false;

  int32 step_amount = (2*cols+2*rows)/100;
  detail::CameraDatumBBoxHelper functor( georef, camera_model,
                                         center_on_zero );

  // Running the edges
  bresenham_apply( BresenhamLine(0,0,cols,0),
                   step_amount, functor );
  functor.last_valid = false;
  bresenham_apply( BresenhamLine(cols-1,0,cols-1,rows),
                   step_amount, functor );
  functor.last_valid = false;
  bresenham_apply( BresenhamLine(cols-1,rows-1,0,rows-1),
                   step_amount, functor );
  functor.last_valid = false;
  bresenham_apply( BresenhamLine(0,rows-1,0,0),
                   step_amount, functor );
  functor.last_valid = false;

  // Running once through the center
  bresenham_apply( BresenhamLine(0,0,cols,rows),
                   step_amount, functor );

  scale = functor.scale/double(step_amount);
  return functor.box;
}
