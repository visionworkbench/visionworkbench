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


/// \file OrbitalCorrections.h
///
/// This file contains algorithms for velocity aberration and atmospheric
/// refraction correction. 

/// These currently only work for Earth.

#ifndef __VW_CAMERA_ORBITALCORRECTIONS_H__
#define __VW_CAMERA_ORBITALCORRECTIONS_H__

#include <vw/Math/Vector.h>
#include <vw/Math/Quaternion.h>

namespace vw {

const double DEFAULT_EARTH_RADIUS      = 6371000.0;  // In meters.
const double DEFAULT_SURFACE_ELEVATION = 0.0;

namespace camera {

  /// Returns the velocity corrected to account for the planetary rotation.
  /// - For efficiency, requires the uncorrected look vector at this location.
  Vector3 get_rotation_corrected_velocity(Vector3 const& camera_center,
                                          Vector3 const& camera_velocity,
                                          double         mean_earth_radius,
                                          Vector3 const& uncorrected_vector);

  /// Adjust a pixel vector to account for velocity aberration.
  Vector3 apply_velocity_aberration_correction(Vector3 const& camera_center,
                                               Vector3 const& camera_velocity,
                                               double         mean_earth_radius,
                                               Vector3 const& uncorrected_vector,
                                               vw::Quaternion<double> & corr_rot);

  /// Simple atmospheric atmospheric correction method.
  double saastamoinen_atmosphere_correction(double camera_alt, double ground_alt, double alpha);

  /// Account for atmospheric refraction.
  Vector3 apply_atmospheric_refraction_correction(Vector3 const& camera_center,
                                                  double         mean_earth_radius,
                                                  double         mean_surface_elevation,
                                                  Vector3 const& uncorrected_vector,
                                                  vw::Quaternion<double> & corr_rot);

}} // namespace vw::camera

#endif // __VW_CAMERA_ORBITALCORRECTIONS_H__

