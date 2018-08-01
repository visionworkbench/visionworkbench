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
  
  // 1. Find the distance from the camera to the first
  // intersection of the current ray with the Earth surface.
  Vector3 cam_ctr          = camera_center(pixel);
  double  earth_ctr_to_cam = norm_2(cam_ctr);
  double  cam_angle_cos    = dot_prod(uncorrected_vector, -normalize(cam_ctr));
  double  len_cos          = earth_ctr_to_cam*cam_angle_cos;
  double  earth_rad        = m_mean_earth_radius;
  double  cam_to_surface   = len_cos - sqrt(earth_rad*earth_rad
                                            + len_cos*len_cos
                                            - earth_ctr_to_cam*earth_ctr_to_cam);
  // 2. Account for Earth's rotation  
  const double SECONDS_PER_DAY = 86164.0905;
  Vector3 earth_rotation_vec(0.0, 0.0, 2*M_PI/SECONDS_PER_DAY);
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
  const double LIGHT_SPEED = 299792458.0;
  Vector3 corrected_vector = uncorrected_vector - cam_vel_corr2/LIGHT_SPEED;
  return normalize(corrected_vector);
}


/// Compute the ray correction for atmospheric refraction using the 
///  Saastamoinen equation.
double LinescanModel::saastamoinen_atmosphere_correction(
          double camera_alt, double ground_alt, double alpha) {

  const double METERS_PER_KM = 1000.0;
  double H      = camera_alt / METERS_PER_KM; // In kilometers
  double h      = ground_alt / METERS_PER_KM; // Height of target point over surface in km
  double h_diff = H - h;
  double p1     = (2335.0 / h_diff)*pow(1.0 - 0.02257*h, 5.256);
  double p2     = pow(0.8540, H-11.0) * (82.2 - 521.0/h_diff);
  double K      = (p1 - p2) * pow(10.0,-6.0);

  double delta_alpha = K * tan(alpha);
  return delta_alpha;
}

// WARNING: This currently only works for Earth!
Vector3 LinescanModel::
apply_atmospheric_refraction_correction(Vector2 const& pixel,
                                        Vector3 const& uncorrected_vector) const {
  // Correct for atmospheric refraction
  // - From Saastamoinen, J. (1972), Atmospheric correction for the 
  //   troposphere and stratosphere in radio ranging of satellites.

  // Get some information
  vw::Vector3 cam_ctr          = camera_center(pixel); // ECF camera coords
  vw::Vector3 cam_ctr_norm     = normalize(cam_ctr);
  vw::Vector3 cam_to_earth_center_unit = -1.0 * cam_ctr_norm;
  double  earth_ctr_to_cam     = norm_2(cam_ctr);    // Distance in meters from cam to earth center
  double  cam_to_earth_surface = earth_ctr_to_cam - m_mean_earth_radius;

  // Compute angle alpha and correction angle
  vw::Vector3 u     = normalize(uncorrected_vector);
  double      alpha = acos(dot_prod(cam_to_earth_center_unit,u)); // Both get angle of normalized vectors.

  // There are more sophisticated correction methods but this one is simple and does help.
  double delta_alpha = saastamoinen_atmosphere_correction(cam_to_earth_surface,
                                                          m_mean_surface_elevation, alpha);

  // Rotate the vector by delta_alpha
  vw::Vector3 rotation_axis = normalize(cross_prod(u, cam_to_earth_center_unit));
  vw::Quaternion<double> refraction_rotation(rotation_axis, delta_alpha);
  vw::Vector3 u_prime = refraction_rotation.rotate(u);

  return u_prime;
}

Vector3 LinescanModel::pixel_to_vector(Vector2 const& pixel) const {
  try {
    // Compute local vector from the pixel out of the sensor
    // - m_detector_origin and m_focal_length have been converted into units of pixels
    Vector3 local_vec = get_local_pixel_vector(pixel);
    // Put the local vector in world coordinates using the pose information.
    Vector3 output_vector = camera_pose(pixel).rotate(local_vec);

    if (!m_correct_atmospheric_refraction) 
      output_vector = apply_atmospheric_refraction_correction(pixel, output_vector);

    if (!m_correct_velocity_aberration) 
      return output_vector;
    else
      return apply_velocity_aberration_correction(pixel, output_vector);
      
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

