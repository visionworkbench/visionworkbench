// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <cmath>
#include <gtest/gtest.h>
#include <vw/Math/EulerAngles.h>
#include <vw/tests/config_test.h>

using namespace vw;
using namespace math;

TEST(Euler, basic_euler_rotation) {
  // Should rotate the vector [0 1 0] 90 degrees around the x-axis
  // to [0 0 1] and then 90 degrees around the NEW y-axis to [1 0 0].
  Vector3 test_vec(0,1,0);
  Matrix3x3 rot1 = euler_to_rotation_matrix(M_PI/2, 0, 0, "xyz");
  Matrix3x3 rot2 = euler_to_rotation_matrix(M_PI/2, M_PI/2, 0, "xyz");
  Vector3 result1 = rot1 * test_vec;
  Vector3 result2 = rot2 * test_vec;

  EXPECT_VECTOR_NEAR( Vector3(0,0,1), result1, 1e-16 );
  EXPECT_VECTOR_NEAR( Vector3(1,0,0), result2, 1e-16 );
}

TEST(Euler, euler_angle_xyz) {
  for (unsigned p = 0; p < 360; p+=60){
    for (unsigned o = 0; o < 360; o+=60){
      for (unsigned k = 0; k < 360; k+=60){
        double phi   = p * M_PI / 180;
        double omega = o * M_PI / 180;
        double kappa = k * M_PI / 180;

        vw::Matrix<double,3,3> rot_matrix = euler_to_rotation_matrix( phi, omega, kappa, "xyz");
        Vector3 euler = rotation_matrix_to_euler_xyz( rot_matrix );

        // PASS 1
        vw::Matrix<double,3,3> rot_matrix2 = euler_to_rotation_matrix( euler[0], euler[1], euler[2], "xyz");
        EXPECT_MATRIX_NEAR( rot_matrix2, rot_matrix, 1e-15 );
      }
    }
  }
}

TEST(Euler, euler_angle_yxz) {
  for (unsigned p = 0; p < 360; p+=60){
    for (unsigned o = 0; o < 360; o+=60){
      for (unsigned k = 0; k < 360; k+=60){
        double phi = p * M_PI / 180;
        double omega = o * M_PI / 180;
        double kappa = k * M_PI / 180;

        vw::Matrix<double,3,3> rot_matrix = euler_to_rotation_matrix( phi, omega, kappa, "yxz");
        Vector3 euler = rotation_matrix_to_euler_yxz( rot_matrix );

        // PASS 1
        vw::Matrix<double,3,3> rot_matrix2 = euler_to_rotation_matrix( euler[0], euler[1], euler[2], "yxz");
        EXPECT_MATRIX_NEAR( rot_matrix2, rot_matrix, 1e-15 );
      }
    }
  }
}

TEST(Euler, euler_angle_zxy) {
  for (unsigned p = 0; p < 360; p+=60){
    for (unsigned o = 0; o < 360; o+=60){
      for (unsigned k = 0; k < 360; k+=60){
        double phi = p * M_PI / 180;
        double omega = o * M_PI / 180;
        double kappa = k * M_PI / 180;

        vw::Matrix<double,3,3> rot_matrix = euler_to_rotation_matrix( phi, omega, kappa, "zxy");
        Vector3 euler = rotation_matrix_to_euler_zxy( rot_matrix );

        // PASS 1
        vw::Matrix<double,3,3> rot_matrix2 = euler_to_rotation_matrix( euler[0], euler[1], euler[2], "zxy");
        EXPECT_MATRIX_NEAR( rot_matrix2, rot_matrix, 1e-15 );
      }
    }
  }
}
