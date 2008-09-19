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

/// \file EulerAngles.h
///
/// Provides a handy set of utilities for working with Euler angles.
///
#ifndef __VW_MATH_EULER_ANGLES_H__
#define __VW_MATH_EULER_ANGLES_H__

#include <vw/Math/Quaternion.h>

namespace vw {
namespace math {

  // Returns: A Vector3 containing the euler angles [phi, omega, kappa] inline
  inline Vector3 rotation_matrix_to_euler_xyz(const Matrix<double,3,3> rotation_matrix) {
    double omega = asin(-rotation_matrix(2,0));
    double kappa = atan2(rotation_matrix(1,0),rotation_matrix(0,0));
    double phi = atan2(rotation_matrix(2,1),rotation_matrix(2,2));

    return Vector3(phi, omega, kappa);
  }

  // Returns: A Vector3 containing the euler angles [phi, omega, kappa]
  inline Vector3 rotation_matrix_to_euler_yxz(const Matrix<double,3,3> rotation_matrix) {
    double cos_phi = sqrt (1 - rotation_matrix(2,1) * rotation_matrix(2,1));
    double phi = atan2(rotation_matrix(2,1), cos_phi);
    double omega = atan2(-rotation_matrix(2,0), rotation_matrix(2,2));
    double kappa = atan2(-rotation_matrix(0,1), rotation_matrix(1,1));
    return Vector3(omega, phi, kappa);
  }
  
  // Returns: A Vector3 containing the euler angles [phi, omega, kappa]
  inline Vector3 rotation_matrix_to_euler_zxy(const Matrix<double,3,3> rotation_matrix) {
    double sin_phi = -rotation_matrix(1,2);
    double cos_phi = sqrt (1 - sin_phi*sin_phi);
    double phi = atan2(sin_phi, cos_phi);
    double omega = atan2(rotation_matrix(0,2), rotation_matrix(2,2));
    double kappa = atan2(rotation_matrix(1,0), rotation_matrix(1,1));
    return Vector3(kappa, phi, omega);
  }

  // Return the rotation matrix for the rotation about the x-axis
  inline vw::Matrix<double,3,3> rotation_x_axis(double theta) {
    vw::Matrix<double,3,3> e;
    e.set_identity();
    e(1,1) = cos(theta);
    e(1,2) = -sin(theta);
    e(2,1) = sin(theta);
    e(2,2) = cos(theta);  
    return e;
  }
  
  // Return the rotation matrix for the rotation about the x-axis
  inline vw::Matrix<double,3,3> rotation_y_axis(double theta) {
    vw::Matrix<double,3,3> e;
    e.set_identity();
    e(0,0) = cos(theta);
    e(0,2) = sin(theta);
    e(2,0) = -sin(theta);
    e(2,2) = cos(theta);
    return e;
  }

  // Return the rotation matrix for the rotation about the x-axis
  inline vw::Matrix<double,3,3> rotation_z_axis(double theta) {
    vw::Matrix<double,3,3> e;
    e.set_identity();
    e(0,0) = cos(theta);
    e(0,1) = -sin(theta);
    e(1,0) = sin(theta);
    e(1,1) = cos(theta);
    return e;
  }
  

  /// \cond INTERNAL
  inline vw::Matrix<double,3,3> euler_rotation_helper(double theta, const char axis) {
    if (axis == 'X' || axis == 'x') 
      return rotation_x_axis(theta);
    else if (axis == 'Y' || axis == 'y') 
      return rotation_y_axis(theta);
    else if (axis == 'Z' || axis == 'z') 
      return rotation_z_axis(theta);
    else 
      vw_throw(vw::ArgumentErr() << "euler_to_quaternion(): unknown axis \"" << axis << "\"\n");
    
    // Should never reach this point.
    return vw::Matrix<double,3,3>();
  }
  /// \endcond


  /// Creates a rotation matrix that represents the same rotation
  /// through the sequence of euler angles, [phi, theta, kappa].  The
  /// euler angles are defined according to the convention specified
  /// in variable 'sequence'.  Sequence is a string containing any
  /// combination of the characters 'x', 'y', and 'z' (thoug the
  /// sequence must always be three characters long) that defines the
  /// axes of rotation for phi, theta, and kappa respectively.  For
  /// example, the sequence "XYZ" would create a rotation of phi
  /// degrees around the X axis, then omega degrees around the new Y
  /// axis, and finally kappa degrees around the new Z axis.
  ///
  /// In matrix notation, this would rotate a vector x_initial as
  /// follows:
  ///
  /// x_final = [ kappa ]_z * [ omega ]_y * [ phi ]_x * x_initial
  ///
  inline vw::Matrix<double,3,3> euler_to_rotation_matrix(double phi, double omega, double kappa, std::string const& sequence) {
    
    VW_ASSERT(sequence.size() == 3,
              vw::ArgumentErr() << "euler_to_rotation_matrix: rotation sequence must be a three character sequence composed of \'x\', \'y\', and \'z\'.");
    
    vw::Matrix<double,3,3> e_phi = euler_rotation_helper(phi, sequence[0]);
    vw::Matrix<double,3,3> e_omega = euler_rotation_helper(omega, sequence[1]);
    vw::Matrix<double,3,3> e_kappa = euler_rotation_helper(kappa, sequence[2]);

    return e_kappa*e_omega*e_phi;
  }

  /// Quaternion variant of euler_to_rotation_matrix()
  inline vw::Quaternion<double> euler_to_quaternion(double phi, double omega, double kappa, std::string const& sequence) {
    return vw::Quaternion<double>(euler_to_rotation_matrix(phi, omega, kappa, sequence));
  }

}} // namespace vw::math

#endif // __VW_MATH_EULER_ANGLES_H__
