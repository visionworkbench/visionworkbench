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

// WARNING: This currently only works for Earth!
Vector3 LinescanModel::get_rotation_corrected_velocity(Vector2 const& pixel,
                                                       Vector3 const& uncorrected_vector) const {
  // TODO: This could be a function that lives somewhere else!
  
  // 1. Find the distance from the camera to the first
  // intersection of the current ray with the Earth surface.
  Vector3 cam_ctr          = camera_center(pixel);
  double  earth_ctr_to_cam = norm_2(cam_ctr);
  double  cam_angle_cos    = dot_prod(uncorrected_vector, -normalize(cam_ctr));
  double  len_cos          = earth_ctr_to_cam*cam_angle_cos;
  double  earth_rad        = 6371000.0; // TODO: Vary by location?
  double  cam_to_surface   = len_cos - sqrt(earth_rad*earth_rad
                                            + len_cos*len_cos
                                            - earth_ctr_to_cam*earth_ctr_to_cam);
  // 2. Account for Earth's rotation  
  double seconds_in_day = 86164.0905;
  Vector3 earth_rotation_vec(0.0, 0.0, 2*M_PI/seconds_in_day);
  Vector3 cam_vel      = camera_velocity(pixel);
  Vector3 cam_vel_corr = cam_vel
    - cam_to_surface * cross_prod(earth_rotation_vec, uncorrected_vector);
  return cam_vel_corr;
}


Vector3 LinescanModel::
apply_velocity_aberration_correction(Vector2 const& pixel,
				     Vector3 const& uncorrected_vector) const {
  
  // 1. Correct the camera velocity due to the fact that the Earth
  // rotates around its axis.
  Vector3 cam_vel_corr1 = get_rotation_corrected_velocity(pixel, uncorrected_vector);

  // 2. Find the component of the camera velocity orthogonal to the
  // direction the camera is pointing to.
  Vector3 cam_vel_corr2 = cam_vel_corr1
    - dot_prod(cam_vel_corr1, uncorrected_vector) * uncorrected_vector;

  // 3. Correct direction for velocity aberration due to the speed of light.
  double light_speed = 299792458.0;
  Vector3 corrected_vector = uncorrected_vector - cam_vel_corr2/light_speed;
  return normalize(corrected_vector);

}

Vector3 LinescanModel::pixel_to_vector(Vector2 const& pixel) const {
  try {
    // Compute local vector from the pixel out of the sensor
    // - m_detector_origin and m_focal_length have been converted into units of pixels
    Vector3 local_vec = get_local_pixel_vector(pixel);
    // Put the local vector in world coordinates using the pose information.
    Vector3 uncorrected_vector = camera_pose(pixel).rotate(local_vec);

    if (!m_correct_velocity_aberration) 
      return uncorrected_vector;
    else
      return apply_velocity_aberration_correction(pixel, uncorrected_vector);
      
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

