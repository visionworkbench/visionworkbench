// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
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
using vw::math::BresenhamLine;

Vector3
cartography::datum_intersection( Datum const& datum,
                                 camera::CameraModel const* model,
                                 Vector2 const& pix ) {
  return datum_intersection(datum, model->camera_center(pix), model->pixel_to_vector(pix));
}

// Return the intersection between the ray emanating from the
// current camera pixel with the datum ellipsoid. The return value
// is a map projected point location (the intermediate between
// lon-lat-altitude and pixel).
Vector2
cartography::geospatial_intersect(cartography::GeoReference const& georef,
                                  Vector3 const& camera_ctr, Vector3 const& camera_vec,
                                  bool& has_intersection) {

  Vector3 intersection = datum_intersection(georef.datum(), camera_ctr, camera_vec);
  if (intersection == Vector3()){
    has_intersection = false;
    return Vector2();
  }
  
  has_intersection = true;
  
  Vector3 llh = georef.datum().cartesian_to_geodetic( intersection );
  Vector2 projected_point = georef.lonlat_to_point( Vector2( llh.x(),
                                                             llh.y() ) );

  return projected_point;
}

namespace vw { namespace cartography { namespace detail {

  /// A class to help identify the extent of an image when
  /// projected onto a datum.
  class CameraDatumBBoxHelper {
    GeoReference m_georef;
    boost::shared_ptr<camera::CameraModel> m_camera;
    Vector2      m_last_intersect;
    std::vector<Vector2> *m_coords;

  public:
    bool   last_valid, center_on_zero;
    BBox2  box;
    double scale;

    CameraDatumBBoxHelper( GeoReference const& georef,
			   boost::shared_ptr<camera::CameraModel> camera,
			   bool center=false,
			   std::vector<Vector2> *coords=0):
      m_georef(georef), m_camera(camera), m_coords(coords), last_valid(false),
      center_on_zero(center), scale( std::numeric_limits<double>::max() ) {
      if (m_coords)
        m_coords->clear();
    }

    void operator()( Vector2 const& pixel ) {
      bool has_intersection;
      Vector2 point = geospatial_intersect(m_georef,
			                                     m_camera->camera_center(pixel),
			                                     m_camera->pixel_to_vector(pixel),
			                                     has_intersection);
      if ( !has_intersection ) {
        last_valid = false;
        return;
      }
      
      if (!m_georef.is_projected()){
        // If we don't use a projected coordinate system, then the
        // coordinates of this point are simply lon and lat.
        if ( center_on_zero && point[0] > 180 )
          point[0] -= 360.0;
      }
      if (m_coords)
        m_coords->push_back(point);
      
      if ( last_valid ) {
        double current_scale = norm_2( point - m_last_intersect );
        if ( current_scale < scale )
          scale = current_scale;
      }
      m_last_intersect = point;
      box.grow( point );
      last_valid = true;
    }
    
  }; // End class CameraDatumBBoxHelper

  // TODO: Why is this not done by default?
  void recenter_point(bool center_on_zero, GeoReference const& georef, Vector2 & point){
    if (!georef.is_projected()){
      // If we don't use a projected coordinate system, then the
      // coordinates of this point are simply lon and lat.
      // - Normalize the longitude coordinate.
      if ( center_on_zero && point[0] > 180 )
        point[0] -= 360.0;
      else if ( center_on_zero && point[0] < -180 )
        point[0] += 360.0;
      else if ( !center_on_zero && point[0] < 0 )
        point[0] += 360.0;
      else if ( !center_on_zero && point[0] > 360 )
        point[0] -= 360.0;
    }
  }
  
}}}
  
BBox2 cartography::camera_bbox( cartography::GeoReference const& georef,
                                boost::shared_ptr<camera::CameraModel> camera_model,
                                int32 cols, int32 rows, float &scale,
                                std::vector<Vector2> *coords ) {

  // Testing to see if we should be centering on zero
  bool center_on_zero = true;
  Vector3 camera_llr =
    georef.datum().cartesian_to_geodetic(camera_model->camera_center(Vector2()));
  if ( camera_llr[0] < -90 || camera_llr[0] > 90 )
    center_on_zero = false;

  int32 step_amount = (2*cols+2*rows)/100;
  step_amount = std::min(step_amount, cols/4); // ensure at least 4 pts/col
  step_amount = std::min(step_amount, rows/4); // ensure at least 4 pts/row
  step_amount = std::max(step_amount, 1);      // step amount must be > 0
  detail::CameraDatumBBoxHelper functor(georef, camera_model, center_on_zero, coords);

  // Running the edges. Note: The last valid point on a
  // BresenhamLine is the last point before the endpoint.
  bresenham_apply( BresenhamLine(0,      0,      cols,   0     ), step_amount, functor );
  functor.last_valid = false;
  bresenham_apply( BresenhamLine(cols-1, 0,      cols-1, rows  ), step_amount, functor );
  functor.last_valid = false;
  bresenham_apply( BresenhamLine(cols-1, rows-1, 0,      rows-1), step_amount, functor );
  functor.last_valid = false;
  bresenham_apply( BresenhamLine(0,      rows-1, 0,      0     ), step_amount, functor );
  functor.last_valid = false;

  // Do the x pattern
  bresenham_apply( BresenhamLine(0, 0,      cols-1, rows-1), step_amount, functor );
  bresenham_apply( BresenhamLine(0, rows-1, cols-1, 0     ), step_amount, functor );

  scale = functor.scale/double(step_amount);
  return functor.box;
}

