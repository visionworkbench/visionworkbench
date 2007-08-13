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

#include <vw/Cartography/CameraBBox.h>

vw::BBox2 vw::cartography::camera_bbox( vw::GeoReference const& georef, double min_alt, double max_alt,
                                        boost::shared_ptr<vw::camera::CameraModel> camera_model, int32 cols, int32 rows ) {
    
  double semi_major_axis = georef.datum().semi_major_axis();
  double semi_minor_axis = georef.datum().semi_minor_axis();
  double z_scale = semi_major_axis / semi_minor_axis;
  double radius = semi_major_axis;
  double radius_sqr = radius*radius;

  BBox2 bbox_180, bbox_360;
  bool first_sign=false, pole=true;
  Vector3 last_intersection;

  // Top row
  for( int x=0; x<cols; ++x ) {
    Vector2 pix(x,0);
    
    Vector3 center = camera_model->camera_center( pix );
    Vector3 vector = camera_model->pixel_to_vector( pix );
    center.z() *= z_scale;
    vector.z() *= z_scale;
    vector = normalize( vector );

    double alpha = - dot_prod(center,vector);
    Vector3 projection = center + alpha * vector;
    if( norm_2_sqr(projection) > radius_sqr ) continue; // the ray does not intersect the sphere
    alpha -= sqrt( radius_sqr - norm_2_sqr(projection) );
    Vector3 intersection = center + alpha * vector;
    intersection.z() /= z_scale;
      
    bool sign = ( (intersection.x()-last_intersection.x())*intersection.y() - (intersection.y()-last_intersection.y())*intersection.x() ) > 0;
    last_intersection = intersection;
    if( x==1 ) first_sign = sign;
    if( x>1 && sign!=first_sign ) pole = false;
    
    Vector3 llr = xyz_to_lon_lat_radius( intersection );
    bbox_180.grow( Vector2(llr.x(),llr.y()) );
    llr.x() += 180.0;
    bbox_360.grow( Vector2(llr.x(),llr.y()) );
  }
  // Bottom row
  for( int x=cols-1; x>=0; --x ) {
    Vector2 pix(x,rows-1);
    
    Vector3 center = camera_model->camera_center( pix );
    Vector3 vector = camera_model->pixel_to_vector( pix );
    center.z() *= z_scale;
    vector.z() *= z_scale;
    vector = normalize( vector );
    
    double alpha = - dot_prod(center,vector);
    Vector3 projection = center + alpha * vector;
    if( norm_2_sqr(projection) > radius_sqr ) continue; // the ray does not intersect the sphere
    alpha -= sqrt( radius_sqr - norm_2_sqr(projection) );
    Vector3 intersection = center + alpha * vector;
    intersection.z() /= z_scale;
    
    bool sign = ( (intersection.x()-last_intersection.x())*intersection.y() - (intersection.y()-last_intersection.y())*intersection.x() ) > 0;
    last_intersection = intersection;
    if( x>0 && sign!=first_sign ) pole = false;
    
    Vector3 llr = xyz_to_lon_lat_radius( intersection );
    bbox_180.grow( Vector2(llr.x(),llr.y()) );
    llr.x() += 180.0;
    bbox_360.grow( Vector2(llr.x(),llr.y()) );
  }
  // Left side
  for( int y=rows-2; y>0; --y ) {
    Vector2 pix(0,y);
    
    Vector3 center = camera_model->camera_center( pix );
    Vector3 vector = camera_model->pixel_to_vector( pix );
    center.z() *= z_scale;
    vector.z() *= z_scale;
    vector = normalize( vector );
    
    double alpha = - dot_prod(center,vector);
    Vector3 projection = center + alpha * vector;
    if( norm_2_sqr(projection) > radius_sqr ) continue; // the ray does not intersect the sphere
    alpha -= sqrt( radius_sqr - norm_2_sqr(projection) );
    Vector3 intersection = center + alpha * vector;
    intersection.z() /= z_scale;
    
    bool sign = ( (intersection.x()-last_intersection.x())*intersection.y() - (intersection.y()-last_intersection.y())*intersection.x() ) > 0;
    last_intersection = intersection;
    if( y<rows-2 && sign!=first_sign ) pole = false;
    
    Vector3 llr = xyz_to_lon_lat_radius( intersection );
    bbox_180.grow( Vector2(llr.x(),llr.y()) );
    llr.x() += 180.0;
    bbox_360.grow( Vector2(llr.x(),llr.y()) );
  }
  // Right side
  for( int y=1; y<rows-1; ++y ) {
    Vector2 pix(cols-1,y);
    
    Vector3 center = camera_model->camera_center( pix );
    Vector3 vector = camera_model->pixel_to_vector( pix );
    center.z() *= z_scale;
    vector.z() *= z_scale;
    vector = normalize( vector );
    
    double alpha = - dot_prod(center,vector);
    Vector3 projection = center + alpha * vector;
    if( norm_2_sqr(projection) > radius_sqr ) continue; // the ray does not intersect the sphere
    alpha -= sqrt( radius_sqr - norm_2_sqr(projection) );
    Vector3 intersection = center + alpha * vector;
    intersection.z() /= z_scale;
    
    bool sign = ( (intersection.x()-last_intersection.x())*intersection.y() - (intersection.y()-last_intersection.y())*intersection.x() ) > 0;
    last_intersection = intersection;
    if( y>1 && sign!=first_sign ) pole = false;
    
    Vector3 llr = xyz_to_lon_lat_radius( intersection );
    bbox_180.grow( Vector2(llr.x(),llr.y()) );
    llr.x() += 180.0;
    bbox_360.grow( Vector2(llr.x(),llr.y()) );
  }
  
  BBox2 bbox = bbox_180;
  // FIXME: This may break for bboxes greater than 180 degrees wide.
  if( bbox_180.min().x() < 0 && bbox_180.max().x() > 0 && bbox_360.max().x() < 360 ) bbox = bbox_360;
  if( pole ) {
    bbox.min().x() = -180;
    bbox.max().x() = 180;
    if( bbox.min().y() > -bbox.max().y() ) bbox.max().y() = 90;
    else bbox.min().y() = -90;
  }
  
  return bbox;
  
}
