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

#include <vw/Camera/OrbitingPushbroomModel.h>
#include <vw/Math/LinearAlgebra.h>
#include <math.h>

namespace vw {
namespace camera {


  // Fit a 3d parametric curve (2nd degree polynomial)
  // to a serios of N points in the Nx3 matrix points.
  // 
  // The curve is specified by the coefficient
  // matrix, which has the form:
  // 
  //          | a b c |
  //   A  =   | d e f |
  //          | g h i |
  // 
  // and the curve is is evaluatied using the equation:
  // 
  //              |  1  |
  //   x(t) = A * |  t  |
  //              | t^2 |
  //
  static Matrix<double,3,3> fit_curve_3d(std::vector<Vector3> const& points, const double t0, const double dt) {
    Matrix<double> Z(points.size()*3, 9);
    //    fill(Z, 0.0);
  
    Vector<double> p(points.size() * 3);
    // Reshape the points matrix into a column vector 
    for (unsigned int i = 0; i < points.size(); i++) {
      p(3*i) = points[i][0];
      p(3*i+1) = points[i][1];
      p(3*i+2) = points[i][2];
    }

    Vector<double> t(points.size());
    for (unsigned int i = 0; i < t.size(); i++) {
      t(i) = t0 + dt*i;
    }

    // Populate the Z matrix 
    for (unsigned int i = 0; i < points.size(); i++) {
      Z(3*i  , 0) = 1.0;
      Z(3*i  , 1) = t(i);
      Z(3*i  , 2) = t(i)*t(i);
      Z(3*i+1, 3) = 1.0;
      Z(3*i+1, 4) = t(i);
      Z(3*i+1, 5) = t(i)*t(i);
      Z(3*i+2, 6) = 1.0;
      Z(3*i+2, 7) = t(i);
      Z(3*i+2, 8) = t(i)*t(i);
    }

    Vector<double> x = least_squares(Z,p);
    //    std::cout << Z << "    " << "\n" << p << "\n" << x << "\n";
    Matrix<double,3,3> coeff;
    coeff(0,0) = x(0);  coeff(0,1) = x(1); coeff(0,2) = x(2);
    coeff(1,0) = x(3);  coeff(1,1) = x(4); coeff(1,2) = x(5);
    coeff(2,0) = x(6);  coeff(2,1) = x(7); coeff(2,2) = x(8);
 
    return coeff;
  } 

  static Quaternion<double> SLERP(double alpha, Quaternion<double> const& a, Quaternion<double> const& b, int spin) {
    const double SLERP_EPSILON = 1.0E-6; 	        // a tiny number
    double beta;			// complementary interp parameter 
    double theta;			// angle between A and B 
    double sin_t, cos_t;		// sine, cosine of theta 
    double phi;			// theta plus spins 
    int bflip;			// use negation of B? 

    // cosine theta = dot product of A and B 
    cos_t = a(1)*b(1) + a(2)*b(2) + a(3)*b(3) + a(0)*b(0);
  
    /* if B is on opposite hemisphere from A, use -B instead */
    if (cos_t < 0.0) {
      cos_t = -cos_t;
      bflip = true;
    } else {
      bflip = false;
    }

    // if B is (within precision limits) the same as A,
    // just linear interpolate between A and B.
    // Can't do spins, since we don't know what direction to spin.
    if (1.0 - cos_t < SLERP_EPSILON) {
      beta = 1.0 - alpha;
    } else {				/* normal case */
      theta = acos(cos_t);
      phi = theta + spin * M_PI;
      sin_t = sin(theta);
      beta = sin(theta - alpha*phi) / sin_t;
      alpha = sin(alpha*phi) / sin_t;
    }
  
    if (bflip)
      alpha = -alpha;
  
    // interpolate 
    Quaternion<double> q;
    q(1) = beta*a(1) + alpha*b(1);
    q(2) = beta*a(2) + alpha*b(2);
    q(3) = beta*a(3) + alpha*b(3);
    q(0) = beta*a(0) + alpha*b(0);
    return q;
  }

  // Evaluate the curve at the supplied time t.
  // 
  // The curve is evaluated using the spherical linear interpolation
  // routine above.
  static Vector3 eval_curve_3d(Matrix<double,3,3> const& coeff, double const t) {

    Vector3 T;
    T(0) = 1;
    T(1) = t;
    T(2) = t*t;

    return coeff * T;
  }

  // Given a sequence of evenly spaced quaternion keyframes seperated by dt
  // and starting at t0, extrapolate the pointing of the quaternion for an 
  // arbitrary time value.  The contiuous valu eof the quaternion is the result 
  // of a sphereical linear interpolation (SLERP) algorithm.
  static Quaternion<double> eval_quat_3d(std::vector<Quaternion<double> > const& quaternions, double t0, double dt, double t) {
 
    // Make sure that t lies within the range [t0, t0+dt*length(points)] 
    if ((t < t0) || (t > t0+dt*quaternions.size())) {
      std::cout << "Time: " << t << "   min: " << t0 << "   max: " << (t0+dt*quaternions.size()) <<"\n";
      throw ArgumentErr() << "Cannot extrapolate point for time given time. Out of valid range.";
    }

    unsigned int low_ind = (unsigned int)floor( (t-t0) / dt );
    unsigned int high_ind = (unsigned int)ceil( (t-t0) / dt );

    // If there are not enough ephemeris points to interpolate at the end,
    // we will limit the high_ind here.
    if (high_ind > quaternions.size()) {
      throw ArgumentErr() << "Attempted to interpolate a quaternion past the last available control point.";
    } else if (high_ind == quaternions.size()) {
      high_ind = quaternions.size() - 1;
    }
  
    double low_t =  t0 + dt * low_ind; 
    double norm_t = (t - low_t)/dt;
  
    Quaternion<double> a = quaternions[low_ind];
    Quaternion<double> b = quaternions[high_ind];
  
    return SLERP(norm_t,a,b,0);
  }
  // ------------------------------------------------
  //           Constructors / Destructors            
  // ------------------------------------------------
  OrbitingPushbroomModel::OrbitingPushbroomModel( double scan_duration,
                                                  double number_of_lines, 
                                                  double samples_per_line, 
                                                  double sample_offset, 
                                                  double focal_length, 
                                                  double along_scan_pixel_size,
                                                  double across_scan_pixel_size,
                                                  double t0_camera_pose,
                                                  double dt_camera_pose,
                                                  double t0_position,
                                                  double dt_position,
                                                  std::vector<Quaternion<double> > const& camera_poses,
                                                  std::vector<Vector3> const& positions,
                                                  std::vector<Vector3> const& velocity_vectors) {
  
    m_scan_duration = scan_duration;
    m_number_of_lines = number_of_lines;
    m_samples_per_line = samples_per_line;
    m_sample_offset = sample_offset;
    m_focal_length = focal_length;
    m_along_scan_pixel_size = along_scan_pixel_size;
    m_across_scan_pixel_size = across_scan_pixel_size;
  
    // Orbiting pushbroom specific
    set_positions(positions, t0_position, dt_position);
    set_velocity_vectors(velocity_vectors, t0_position, dt_position);
    set_camera_poses(camera_poses, t0_camera_pose, dt_camera_pose);
  }
  
  // ------------------------------------------------
  //           interface
  // ------------------------------------------------
  Vector2 OrbitingPushbroomModel::point_to_pixel(Vector3 const& vec) const {
    throw vw::NoImplErr() << "OrbitingPushbroomModel::point_to_pixel is not yet implemented.";
  }
  
  Vector3 OrbitingPushbroomModel::pixel_to_vector(Vector2 const& pix) const {
    // u is the crosstrack (i.e., perspective) direction
    // v is the downtrack (i.e., orthographic) direction
    double u = pix[0], v = pix[1];

    // Compute the time for the given pixel
    double time_per_line = m_scan_duration / m_number_of_lines;
    double pixel_time = v * time_per_line;

    // Compute the pose of the camera at pixel_time.
    Quaternion<double> pose = eval_quat_3d(m_camera_poses, 
                                           m_t0_camera_pose, 
                                           m_dt_camera_pose,
                                           pixel_time);
    
    // The view_matrix takes vectors from the camera (extrinsic)
    // coordinate system to the world frame
    //
    // The position and veloctiy are not actually needed, since we are
    // purely interested in returning the direction of the ray at this
    // point and not its origin.
    Matrix<double,3,3> rotation_matrix = transpose(pose.rotation_matrix());
    
    // The viewplane is the y-z plane of the camera coordinate system.
    // Assuming the origin of the coordinate system is at the center
    // of projection, the image plane is z = +f, and the pixel
    // position in camera coordinates is:
    double pixel_pos_y = (u + m_sample_offset) * m_across_scan_pixel_size;
    double f = m_focal_length;
    Vector<double, 3> pixel_pos(0.0, pixel_pos_y, f);

    // Transform to world coordinates using the rigid rotation
    Vector<double, 3> direction_vec = rotation_matrix * pixel_pos;
    return normalize(Vector3 ( direction_vec[0],
                               direction_vec[1], 
                               direction_vec[2] ));
  }

  
  Vector3 OrbitingPushbroomModel::camera_center(Vector2 const& pix ) const {
    double t = pix[1] / m_number_of_lines * m_scan_duration;
    return eval_curve_3d(m_position_coeff, t);
  }

  Vector3 OrbitingPushbroomModel::camera_position (double t) const {
    return eval_curve_3d(m_position_coeff, t);
  }

  Quaternion<double> OrbitingPushbroomModel::camera_pose (double t) const {
    return eval_quat_3d(m_camera_poses, m_t0_camera_pose, 
                        m_dt_camera_pose, t);
  }
  
  

  // ------------------------------------------------
  //                 Public Methods                   
  // ------------------------------------------------]
  const Matrix<double>& OrbitingPushbroomModel::get_camera_matrix(Vector2 const& pix ) {
    double t = pix.x() / m_number_of_lines * m_scan_duration;

    VW_ASSERT( m_camera_poses.size() != 0 && 
               m_camera_poses.size() == m_positions.size() &&
               m_camera_poses.size() == m_velocity_vectors.size(),
               LogicErr() << "OrbitingPushbroomModel: Improperly initialized -- vectors for pose, position and velocity must be non-zero and the same size.");

    m_camera_pose = eval_quat_3d(m_camera_poses,
                                 m_t0_camera_pose, 
                                 m_dt_camera_pose, 
                                 t);
    m_initial_position = eval_curve_3d(m_position_coeff, t);
    m_velocity_vector = eval_curve_3d(m_velocity_vector_coeff, t);
    
    LinearPushbroomModel::recompute_camera_matrix();
    return LinearPushbroomModel::get_camera_matrix(pix);

  }

  void OrbitingPushbroomModel::get_decomposed_camera_matrix(Matrix<double> &projection,
                                                            Matrix<double> &velocity_matrix,
                                                            Matrix<double> &similarity, 
                                                            Vector2 const& pix ) {
    double t = pix.x() / m_number_of_lines * m_scan_duration;

    VW_ASSERT( m_camera_poses.size() != 0 && 
               m_camera_poses.size() == m_positions.size() &&
               m_camera_poses.size() == m_velocity_vectors.size(),
               LogicErr() << "OrbitingPushbroomModel: Improperly initialized -- vectors for pose, position and velocity must be non-zero and the same size.");

    m_camera_pose = eval_quat_3d(m_camera_poses,
                                 m_t0_camera_pose, 
                                 m_dt_camera_pose, 
                                 t);
    m_initial_position = eval_curve_3d(m_position_coeff, t);
    m_velocity_vector = eval_curve_3d(m_velocity_vector_coeff, t);
    
    LinearPushbroomModel::recompute_camera_matrix();
    LinearPushbroomModel::get_decomposed_camera_matrix(projection, 
                                                       velocity_matrix, 
                                                       similarity, 
                                                       pix);
  }

  void OrbitingPushbroomModel::set_positions( std::vector<Vector3> const& val, double t0, double dt) { 
    m_t0_position = t0;
    m_dt_position = dt;
    m_positions = val; 
    m_position_coeff = fit_curve_3d(m_positions, m_t0_position, m_dt_position);
  }
  
  void OrbitingPushbroomModel::set_velocity_vectors( std::vector<Vector3> const& val, double t0, double dt) { 
    m_t0_position = t0;
    m_dt_position = dt;
    m_velocity_vectors = val;
    m_velocity_vector_coeff = fit_curve_3d(m_velocity_vectors, m_t0_position, m_dt_position);
  }

  std::ostream& operator<<( std::ostream& os, OrbitingPushbroomModel const& camera_model) {
    os << "\n-------------------- Linear Pushbroom Camera Model -------------------\n\n";
    os << " Camera Pose      :   " << camera_model.camera_pose() << "\n";
    os << " Initial Position :   " << camera_model.initial_position() << "\n";
    os << " Velocity Vector  :   " << camera_model.velocity_vector() << "\n";
    os << " Scan Duration    :   " << camera_model.scan_duration() << "\n";
    os << " Number of Lines  :   " << camera_model.number_of_lines() << "\n";
    os << " Samples per Line :   " << camera_model.samples_per_line() << "\n";
    os << " Focal Length     :   " << camera_model.focal_length() << "\n";
    os << " Sample Offset    :   " << camera_model.sample_offset() << "\n";
    os << " Across Scan Pixel Size :   " << camera_model.across_scan_pixel_size() << "\n";
    os << " Along Scan Pixel Size  :   " << camera_model.along_scan_pixel_size() << "\n";
    os << "\n------------------------------------------------------------------------\n\n";
  }

}} // namespace vw::camera

