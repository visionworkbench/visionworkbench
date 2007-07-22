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

#ifndef __VW_CARTOGRAPHY_CAMERABBOX_H__
#define __VW_CARTOGRAPHY_CAMERABBOX_H__

#if defined(VW_HAVE_PKG_CAMERA) && (VW_HAVE_PKG_CAMERA==1)

#include <vw/Cartography/GeoReference.h>
#include <vw/Camera/CameraModel.h>

#include <boost/shared_ptr.hpp>

namespace vw {
namespace cartography {

  BBox2 camera_bbox( GeoReference const& georef, double min_alt, double max_alt,
                     boost::shared_ptr<vw::camera::CameraModel> camera_model, int32 cols, int32 rows ) {
    
    double semi_major_axis = georef.datum().semi_major_axis();
    double semi_minor_axis = georef.datum().semi_minor_axis();
    double z_scale = semi_major_axis / semi_minor_axis;
    double radius = semi_major_axis;
    double radius_sqr = radius*radius;

    BBox2 bbox;

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
      
      Vector3 llr = xyz_to_lon_lat_radius( intersection );
      bbox.grow( georef.lonlat_to_pixel( Vector2(llr.x(),llr.y()) ) );
    }
    // Bottom row
    for( int x=0; x<cols; ++x ) {
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
      
      Vector3 llr = xyz_to_lon_lat_radius( intersection );
      bbox.grow( georef.lonlat_to_pixel( Vector2(llr.x(),llr.y()) ) );
    }
    // Left side
    for( int y=1; y<rows-1; ++y ) {
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
      
      Vector3 llr = xyz_to_lon_lat_radius( intersection );
      bbox.grow( georef.lonlat_to_pixel( Vector2(llr.x(),llr.y()) ) );
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
      
      Vector3 llr = xyz_to_lon_lat_radius( intersection );
      bbox.grow( georef.lonlat_to_pixel( Vector2(llr.x(),llr.y()) ) );
    }

    return bbox;

  }

  inline BBox2 camera_bbox( GeoReference const& georef, boost::shared_ptr<vw::camera::CameraModel> camera_model, int32 cols, int32 rows ) {
    return camera_bbox( georef, 0.0, 0.0, camera_model, cols, rows );
  }


} // namespace cartography
} // namespace vw

#endif // VW_HAVE_PKG_CAMERA

#endif // __VW_CARTOGRAPHY_CAMERABBOX_H__
