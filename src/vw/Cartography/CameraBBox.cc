// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Cartography/CameraBBox.h>
#include <vw/Cartography/PointImageManipulation.h>

using namespace vw;

Vector2 geospatial_intersect( Vector2 pix,
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

  double semi_major_axis = georef.datum().semi_major_axis();
  double semi_minor_axis = georef.datum().semi_minor_axis();
  double z_scale = semi_major_axis / semi_minor_axis;
  double radius = semi_major_axis;
  double radius_sqr = radius*radius;

  BBox2 georeference_space_bbox;
  bool first_sign=false, pole=true, last_valid=false;
  Vector3 last_intersection;
  Vector2 last_geospatial_point;
  scale = -1;

  // Top row
  for( int x=0; x<cols; ++x ) {
    Vector2 pix(x,0);
    bool test_intersect;
    Vector2 geospatial_point = geospatial_intersect( pix, georef,
                                                     camera_model,
                                                     z_scale,
                                                     test_intersect );
    if ( !test_intersect ) {
      last_valid = false;
      continue;
    }

    if( last_valid ) {
      double current_scale = norm_2( geospatial_point - last_geospatial_point );
      if ( current_scale < 0 ||
           current_scale < scale )
        scale = current_scale;
    }
    last_geospatial_point = geospatial_point;

    georeference_space_bbox.grow( geospatial_point );

    last_valid = true;
  }
  // Bottom row
  for( int x=cols-1; x>=0; --x ) {
    Vector2 pix(x,rows-1);

    bool test_intersect;
    Vector2 geospatial_point = geospatial_intersect( pix, georef,
                                                     camera_model,
                                                     z_scale,
                                                     test_intersect );
    if ( !test_intersect ) {
      last_valid = false;
      continue;
    }

    if( last_valid ) {
      double current_scale = norm_2( geospatial_point - last_geospatial_point );
      if ( current_scale < 0 ||
           current_scale < scale )
        scale = current_scale;
    }
    last_geospatial_point = geospatial_point;

    georeference_space_bbox.grow( geospatial_point );

    last_valid = true;
  }
  // Left side
  for( int y=rows-2; y>0; --y ) {
    Vector2 pix(0,y);

    bool test_intersect;
    Vector2 geospatial_point = geospatial_intersect( pix, georef,
                                                     camera_model,
                                                     z_scale,
                                                     test_intersect );
    if ( !test_intersect ) {
      last_valid = false;
      continue;
    }

    if( last_valid ) {
      double current_scale = norm_2( geospatial_point - last_geospatial_point );
      if ( current_scale < 0 ||
           current_scale < scale )
        scale = current_scale;
    }
    last_geospatial_point = geospatial_point;

    georeference_space_bbox.grow( geospatial_point );

    //bbox_180.grow( lonlat );
    //lonlat.x() += 360.0;
    //bbox_360.grow( lonlat );
    last_valid = true;
  }
  // Right side
  last_valid = false;
  for( int y=1; y<rows-1; ++y ) {
    Vector2 pix(cols-1,y);

    bool test_intersect;
    Vector2 geospatial_point = geospatial_intersect( pix, georef,
                                                     camera_model,
                                                     z_scale,
                                                     test_intersect );
    if ( !test_intersect ) {
      last_valid = false;
      continue;
    }

    if( last_valid ) {
      double current_scale = norm_2( geospatial_point - last_geospatial_point );
      if ( scale < 0 ||
           current_scale < scale )
        scale = current_scale;
    }
    last_geospatial_point = geospatial_point;

    georeference_space_bbox.grow( geospatial_point );

    last_valid = true;
  }

  BBox2 bbox = georeference_space_bbox;

  // Did we fail to find scale?
  if ( scale == -1 ) {
    Vector2 pix(cols,rows);
    Vector2 pix2 = pix + Vector2(1,1);

    bool test_intersect, test_intersect2;
    Vector2 geospatial_point = geospatial_intersect( pix, georef,
                                                     camera_model,
                                                     z_scale,
                                                     test_intersect );
    Vector2 geospatial_point2 = geospatial_intersect( pix2, georef,
                                                      camera_model,
                                                      z_scale,
                                                      test_intersect2 );
    if ( (!test_intersect) || (!test_intersect2) )
      return bbox;

    scale = norm_2( geospatial_point - geospatial_point2 );
  }

  // If the bbox_180 crosses the singularity at 180 degrees, then we
  // return the bbox_360 instead.
  // FIXME: We do not properly handle cases where the correct bounding
  // box includes both the prime meridian and the antimeridean but
  // does not include the pole.  This is a highly unlikely corner
  // case that can only occur under very particular circumstances
  // involving camera distortion and near alignment of one side of the
  // bounding box with the prime meridian over a pole.
  // if( (bbox_180.min().x() < 180.0 && bbox_180.max().x() > 180.0) ||
//       (bbox_180.min().x() < -180.0 && bbox_180.max().x() > -180.0) )
//     bbox = bbox_360;
//   if( pole ) {
//     bbox.min().x() = -180;
//     bbox.max().x() = 180;
//     if( bbox.min().y() > -bbox.max().y() ) bbox.max().y() = 90;
//     else bbox.min().y() = -90;
//   }

  return bbox;

}
