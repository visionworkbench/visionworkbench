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

#include <vw/Camera/OrbitalCorrections.h>

/// This file contains algorithms for velocity aberration and atmospheric
/// refraction correction. 

/// These currently only work for Earth.

namespace vw {
namespace camera {

Vector3 get_rotation_corrected_velocity(Vector3 const& camera_center,
                                        Vector3 const& camera_velocity,
                                        double mean_earth_radius,
                                        Vector3 const& uncorrected_vector) {
  
  // 1. Find the distance from the camera to the first
  // intersection of the current ray with the Earth surface.
  double  earth_ctr_to_cam = norm_2(camera_center);
  double  cam_angle_cos    = dot_prod(uncorrected_vector, -normalize(camera_center));
  double  len_cos          = earth_ctr_to_cam*cam_angle_cos;
  double  earth_rad        = mean_earth_radius;
  double  cam_to_surface   = len_cos - sqrt(earth_rad*earth_rad
                                            + len_cos*len_cos
                                            - earth_ctr_to_cam*earth_ctr_to_cam);
  // 2. Account for Earth's rotation  
  const double SECONDS_PER_DAY = 86164.0905;
  Vector3 earth_rotation_vec(0.0, 0.0, 2*M_PI/SECONDS_PER_DAY);
  Vector3 cam_vel_corr = camera_velocity
    - cam_to_surface * cross_prod(earth_rotation_vec, uncorrected_vector);
  return cam_vel_corr;
}

Vector3 apply_velocity_aberration_correction(Vector3 const& camera_center,
                                             Vector3 const& camera_velocity,
                                             double         mean_earth_radius,
                                             Vector3 const& uncorrected_vector,
                                             vw::Quaternion<double> & corr_rot) {

  // 1. Correct the camera velocity due to the fact that the Earth
  // rotates around its axis.
  Vector3 cam_vel_corr1 = get_rotation_corrected_velocity(camera_center, camera_velocity,
                                                          mean_earth_radius, uncorrected_vector);

  // 2. Find the component of the camera velocity orthogonal to the
  // direction the camera is pointing to.
  Vector3 cam_vel_corr2 = cam_vel_corr1
    - dot_prod(cam_vel_corr1, uncorrected_vector) * uncorrected_vector;

  // 3. Compute the correct direction for velocity aberration due to the speed of light.
  const double LIGHT_SPEED = 299792458.0;
  vw::Vector3 corr = -cam_vel_corr2/LIGHT_SPEED;
  
  // 4. Find the rotation that will rotate the uncorrected vector to the corrected one.
  // Axis of rotation
  vw::Vector3 rotation_axis = normalize(cross_prod(uncorrected_vector, corr));
  // Angle of rotation
  double angle = atan(norm_2(corr) / norm_2(uncorrected_vector));
  // Rotation quaternion
  corr_rot = vw::Quaternion<double>(rotation_axis, angle);

  // 5. Apply the correction  
  Vector3 corrected_vector = uncorrected_vector + corr;
  
  return normalize(corrected_vector);
}

/// Compute the ray correction for atmospheric refraction using the 
///  Saastamoinen equation.
double saastamoinen_atmosphere_correction(double camera_alt, double ground_alt, double alpha) {

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

Vector3 apply_atmospheric_refraction_correction(Vector3 const& camera_center,
                                                double         mean_earth_radius,
                                                double         mean_surface_elevation,
                                                Vector3 const& uncorrected_vector,
                                                vw::Quaternion<double> & corr_rot) {
  // Correct for atmospheric refraction
  // - From Saastamoinen, J. (1972), Atmospheric correction for the 
  //   troposphere and stratosphere in radio ranging of satellites.

  // Get some information
  vw::Vector3 cam_ctr_norm     = normalize(camera_center);
  vw::Vector3 cam_to_earth_center_unit = -1.0 * cam_ctr_norm;
  double  earth_ctr_to_cam     = norm_2(camera_center);   // Distance in meters from cam to earth center
  double  cam_to_earth_surface = earth_ctr_to_cam - mean_earth_radius;

  // Compute angle alpha and correction angle. Normalized vectors must be used for this.
  vw::Vector3 uncorr = normalize(uncorrected_vector);
  double alpha = acos(dot_prod(cam_to_earth_center_unit, uncorr));

  // There are more sophisticated correction methods but this one is simple and does help.
  double delta_alpha = saastamoinen_atmosphere_correction(cam_to_earth_surface,
                                                          mean_surface_elevation, alpha);

  // Rotate the vector by delta_alpha
  vw::Vector3 rotation_axis = normalize(cross_prod(uncorr, cam_to_earth_center_unit));
  corr_rot = vw::Quaternion<double>(rotation_axis, delta_alpha);
  vw::Vector3 corrected = corr_rot.rotate(uncorr);

  return corrected;
}

} // End namespace camera
} // End namespace vw
