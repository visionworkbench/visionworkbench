// __BEGIN_LICENSE__
//  Copyright (c) 2006-2012, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#include <vw/Cartography/CameraBBox.h>

using namespace vw;


// Intersect the ray back-projected from the camera with the datum.
Vector3
cartography::datum_intersection( Datum const& datum,
                                 camera::CameraModel const* model,
                                 Vector2 const& pix ) {
  
  // The datum is a spheroid. To simplify the calculations, scale
  // everything in such a way that the spheroid becomes a
  // sphere. Scale back at the end of computation.
  
  double z_scale = datum.semi_major_axis() / datum.semi_minor_axis();
  
  Vector3 ccenter = model->camera_center( pix );
  Vector3 cvec  = model->pixel_to_vector( pix );
  ccenter.z() *= z_scale;
  cvec.z() *= z_scale;
  cvec = normalize(cvec);
  double radius_2 = datum.semi_major_axis() *
    datum.semi_major_axis();
  double alpha = -dot_prod(ccenter, cvec );
  Vector3 projection = ccenter + alpha*cvec;
  if ( norm_2_sqr(projection) > radius_2 ) {
    // did not intersect
    return Vector3();
  }
  
  alpha -= sqrt( radius_2 -
                 norm_2_sqr(projection) );
  Vector3 intersection = ccenter + alpha * cvec;
  intersection.z() /= z_scale;
  return intersection;
}

// Return the intersection between the ray emanating from the
// current camera pixel with the datum ellipsoid. The return value
// is a map projected point location (the intermediate between
// lon-lat-altitude and pixel).
Vector2
cartography::geospatial_intersect( Vector2 pix,
                                   cartography::GeoReference const& georef,
                                   boost::shared_ptr<camera::CameraModel> camera_model,
                                   bool& has_intersection ) {

  Vector3 intersection = cartography::datum_intersection(georef.datum(), camera_model.get(), pix);
  if (intersection == Vector3()){
    has_intersection = false;
    return Vector2();
  } else {
    has_intersection = true;
  }
  
  Vector3 llh = georef.datum().cartesian_to_geodetic( intersection );
  Vector2 geospatial_point = georef.lonlat_to_point( Vector2( llh.x(),
                                                              llh.y() ) );

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

void cartography::detail::CameraDatumBBoxHelper::operator()( Vector2 const& pixel ) {
  bool has_intersection;
  Vector2 geospatial_point =
    geospatial_intersect( pixel, m_georef, m_camera,
                          has_intersection );
  if ( !has_intersection ) {
    last_valid = false;
    return;
  }

  if (!m_georef.is_projected()){
    // If we don't use a projected coordinate system, then the
    // coordinates of this point are simply lon and lat.
    if ( center_on_zero && geospatial_point[0] > 180 )
      geospatial_point[0] -= 360.0;
  }
  
  if ( last_valid ) {
    double current_scale =
      norm_2( geospatial_point - m_last_intersect );
    if ( current_scale < scale )
      scale = current_scale;
  }
  m_last_intersect = geospatial_point;
  box.grow( geospatial_point );
  last_valid = true;
}
