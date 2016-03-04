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

// This file is included in LinescanModel.h

template <class PositionFuncT, class VelocityFuncT, class PoseFuncT, class TimeFuncT>
Vector2 LinescanModel<PositionFuncT, VelocityFuncT, PoseFuncT, TimeFuncT>
::point_to_pixel(Vector3 const& point) const {
  return point_to_pixel(point, -1);
}

// Here we use an initial guess for the line number
template <class PositionFuncT, class VelocityFuncT, class PoseFuncT, class TimeFuncT>
Vector2 LinescanModel<PositionFuncT, VelocityFuncT, PoseFuncT, TimeFuncT>
::point_to_pixel(Vector3 const& point, double starty) const {

  if (!m_correct_velocity_aberration)
    return point_to_pixel_uncorrected(point, starty);
  return point_to_pixel_corrected(point, starty);
}

/// Constants used by the internal solvers
namespace linescan {
  const double ABS_TOL = 1e-16;
  const double REL_TOL = 1e-16;
  const int    MAX_ITERATIONS = 1e+5;
}

template <class PositionFuncT, class VelocityFuncT, class PoseFuncT, class TimeFuncT>
Vector2 LinescanModel<PositionFuncT, VelocityFuncT, PoseFuncT, TimeFuncT>
::point_to_pixel_uncorrected(Vector3 const& point, double starty) const {

  // Solve for the correct line number to use
  LinescanLMA model( this, point );
  int status;
  Vector<double> objective(1), start(1);
  start[0] = m_image_size.y()/2;

  // Use a refined guess, if available
  if (starty >= 0)
    start[0] = starty;

  Vector<double> solution = math::levenberg_marquardt(model, start, objective, status,
                                                      linescan::ABS_TOL,
                                                      linescan::REL_TOL,
                                                      linescan::MAX_ITERATIONS);

  VW_ASSERT( status > 0, camera::PointToPixelErr() << "Unable to project point into LinescanDG model" );

  // Solve for sample location
  double  t  = m_time_func( solution[0] );
  Vector3 pt = inverse( m_pose_func(t) ).rotate( point - m_position_func(t) );
  pt *= m_focal_length / pt.z();

  return Vector2(pt.x() - m_detector_origin[0], solution[0]);
}

template <class PositionFuncT, class VelocityFuncT, class PoseFuncT, class TimeFuncT>
Vector2 LinescanModel<PositionFuncT, VelocityFuncT, PoseFuncT, TimeFuncT>
::point_to_pixel_corrected(Vector3 const& point, double starty) const {

  LinescanCorrLMA model( this, point );
  int status;
  Vector2 start = point_to_pixel_uncorrected(point, starty);

  Vector3 objective(0, 0, 0);
  Vector2 solution = math::levenberg_marquardt(model, start, objective, status,
                                               linescan::ABS_TOL,
                                               linescan::REL_TOL,
                                               linescan::MAX_ITERATIONS);
  VW_ASSERT( status > 0,
	     camera::PointToPixelErr() << "Unable to project point into LinescanDG model" );

  return solution;
}

template <class PositionFuncT, class VelocityFuncT, class PoseFuncT, class TimeFuncT>
Vector3 LinescanModel<PositionFuncT, VelocityFuncT, PoseFuncT, TimeFuncT>
::pixel_to_vector(Vector2 const& pix) const {

  // Compute local vector from the pixel out of the sensor
  // - m_detector_origin and m_focal_length have been converted into units of pixels
  Vector3 local_vec(pix[0]+m_detector_origin[0], m_detector_origin[1], m_focal_length);
  // Put the local vector in world coordinates using the pose information.
  Vector3 pix_to_vec = normalize(camera_pose(pix).rotate(local_vec));

  if (!m_correct_velocity_aberration) return pix_to_vec;

  // Correct for velocity aberration

  // 1. Find the distance from the camera to the first
  // intersection of the current ray with the Earth surface.
  Vector3 cam_ctr          = camera_center(pix);
  double  earth_ctr_to_cam = norm_2(cam_ctr);
  double  cam_angle_cos    = dot_prod(pix_to_vec, -normalize(cam_ctr));
  double  len_cos          = earth_ctr_to_cam*cam_angle_cos;
  double  earth_rad        = 6371000.0; // TODO: Vary by location?
  double  cam_to_surface   = len_cos - sqrt(earth_rad*earth_rad
					    + len_cos*len_cos
					    - earth_ctr_to_cam*earth_ctr_to_cam);

  // 2. Correct the camera velocity due to the fact that the Earth
  // rotates around its axis.
  double seconds_in_day = 86164.0905;
  Vector3 earth_rotation_vec(0.0, 0.0, 2*M_PI/seconds_in_day);
  Vector3 cam_vel = camera_velocity(pix);
  Vector3 cam_vel_corr1 = cam_vel - cam_to_surface * cross_prod(earth_rotation_vec, pix_to_vec);

  // 3. Find the component of the camera velocity orthogonal to the
  // direction the camera is pointing to.
  Vector3 cam_vel_corr2 = cam_vel_corr1 - dot_prod(cam_vel_corr1, pix_to_vec) * pix_to_vec;

  // 4. Correct direction for velocity aberration due to the speed of light.
  double light_speed = 299792458.0;
  Vector3 corr_pix_to_vec = pix_to_vec - cam_vel_corr2/light_speed;
  return normalize(corr_pix_to_vec);
}


template <class PositionFuncT, class VelocityFuncT, class PoseFuncT, class TimeFuncT>
typename LinescanModel<PositionFuncT, VelocityFuncT, PoseFuncT, TimeFuncT>::LinescanLMA::result_type
LinescanModel<PositionFuncT, VelocityFuncT, PoseFuncT, TimeFuncT>::LinescanLMA
::operator()( domain_type const& y ) const {
  double t = m_model->m_time_func( y[0] );

  // Rotate the point into our camera's frame
  Vector3 pt = inverse( m_model->m_pose_func(t) ).rotate( m_point - m_model->m_position_func(t) );
  pt *= m_model->m_focal_length / pt.z(); // Rescale to pixel units
  result_type result(1);
  result[0] = pt.y() -
    m_model->m_detector_origin[1]; // Error against the location of the detector
  return result;
}



template <class PositionFuncT, class VelocityFuncT, class PoseFuncT, class TimeFuncT>
std::ostream& operator<<( std::ostream& os, LinescanModel<PositionFuncT, VelocityFuncT, PoseFuncT, TimeFuncT> const& camera_model) {
  os << "\n-------------------- Linescan Camera Model -------------------\n\n";
  os << " Camera center @ origin :   " << camera_model.camera_center(Vector2(0,0)) << "\n";
  os << " Number of Lines        :   " << camera_model.number_of_lines()        << "\n";
  os << " Samples per Line       :   " << camera_model.samples_per_line()       << "\n";
  os << " Sample Offset          :   " << camera_model.sample_offset()          << "\n";
  os << " Focal Length           :   " << camera_model.focal_length()           << "\n";
  os << " Across Scan Pixel Size :   " << camera_model.across_scan_pixel_size() << "\n";
  os << " Along Scan Pixel Size  :   " << camera_model.along_scan_pixel_size()  << "\n";
  os << "\n------------------------------------------------------------------------\n\n";
  return os;
}

