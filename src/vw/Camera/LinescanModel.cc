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

#include <vw/Camera/LinescanModel.h>
#include <vw/Camera/CameraSolve.h>

namespace vw {
namespace camera {

Vector2 LinescanModel::point_to_pixel(Vector3 const& point) const {
  return point_to_pixel(point, -1); // Redirect to the function with no guess
}


Vector2 LinescanModel::point_to_pixel(Vector3 const& point, double starty) const {

  // Use the generic solver to find the pixel 
  // - This method will be slower but works for more complicated geometries
  CameraGenericLMA model( this, point );
  int status;
  Vector2 start = m_image_size / 2.0; // Use the center as the initial guess
  if (starty >= 0) // If the user provided a line number guess..
    start[1] = starty;

  // Solver constants
  const double ABS_TOL = 1e-16;
  const double REL_TOL = 1e-16;
  const int    MAX_ITERATIONS = 1e+5;

  Vector3 objective(0, 0, 0);
  Vector2 solution = math::levenberg_marquardt(model, start, objective, status,
                                               ABS_TOL, REL_TOL, MAX_ITERATIONS);
  VW_ASSERT( status > 0,
	     camera::PointToPixelErr() << "Unable to project point into Linescan model" );

  return solution;
}


Vector3 LinescanModel::pixel_to_vector(Vector2 const& pixel) const {
  try {
    // Compute local vector from the pixel out of the sensor
    // - m_detector_origin and m_focal_length have been converted into units of pixels
    Vector3 local_vec = get_local_pixel_vector(pixel);
    // Put the local vector in world coordinates using the pose information.
    Vector3 output_vector = camera_pose(pixel).rotate(local_vec);


    Vector3 cam_ctr = camera_center(pixel);
    if (!m_correct_atmospheric_refraction) 
      output_vector = apply_atmospheric_refraction_correction(cam_ctr, m_mean_earth_radius,
                                                              m_mean_surface_elevation, output_vector);

    if (!m_correct_velocity_aberration) 
      return output_vector;
    else
      return apply_velocity_aberration_correction(cam_ctr, camera_velocity(pixel),
                                                  m_mean_earth_radius, output_vector);

  } catch(const vw::Exception &e) {
    // Repackage any of our exceptions thrown below this point as a 
    //  pixel to ray exception that other code will be able to handle.
    vw_throw(vw::camera::PixelToRayErr() << e.what());
  }
}

/*
std::ostream& operator<<( std::ostream& os, LinescanModel const& camera_model) {
  os << "\n-------------------- Linescan Camera Model -------------------\n\n";
  os << " Number of Lines        :   " << camera_model.number_of_lines()        << "\n";
  os << " Samples per Line       :   " << camera_model.samples_per_line()       << "\n";
  os << "\n------------------------------------------------------------------------\n\n";
  return os;
}
*/

}} // namespace vw::camera

