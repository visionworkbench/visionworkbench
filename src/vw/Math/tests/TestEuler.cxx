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


#include <cmath>
#include <vw/Math/EulerAngles.h>
#include <test/Helpers.h>

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

TEST(Euler, euler_angle_zyx) {
  for (unsigned p = 0; p < 360; p+=60){
    for (unsigned o = 0; o < 360; o+=60){
      for (unsigned k = 0; k < 360; k+=60){
        double phi   = p * M_PI / 180;
        double omega = o * M_PI / 180;
        double kappa = k * M_PI / 180;

        vw::Matrix<double,3,3> rot_matrix = euler_to_rotation_matrix( phi, omega, kappa, "zyx");
        Vector3 euler = rotation_matrix_to_euler_zyx( rot_matrix );

        // PASS 1
        vw::Matrix<double,3,3> rot_matrix2 = euler_to_rotation_matrix( euler[0], euler[1], euler[2], "zyx");
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

TEST(Euler, quat_euler_xyz) {
  for ( unsigned x = 0; x < 360; x+= 60 ) {
    for ( unsigned y = 0; y < 360; y+= 60 ) {
      for ( unsigned z = 0; z < 360; z+= 60 ) {
        Vector3 euler_start( x,y,z );
        euler_start *= M_PI/180.0;
        Quaternion<double> q = euler_xyz_to_quaternion( euler_start );
        Quaternion<double> qalt = euler_to_quaternion( euler_start[0],
                                                       euler_start[1],
                                                       euler_start[2],
                                                       "xyz" );
        EXPECT_NEAR( fabs(q[0]), fabs(qalt[0]), 1e-8 );
        EXPECT_NEAR( fabs(q[1]), fabs(qalt[1]), 1e-8 );
        EXPECT_NEAR( fabs(q[2]), fabs(qalt[2]), 1e-8 );
        EXPECT_NEAR( fabs(q[3]), fabs(qalt[3]), 1e-8 );

        Vector3 euler_fin = quaternion_to_euler_xyz( q );

        Matrix3x3 rot_start = euler_to_rotation_matrix( euler_start[0],
                                                 euler_start[1],
                                                 euler_start[2], "xyz" );
        Matrix3x3 rot_end = euler_to_rotation_matrix( euler_fin[0],
                                               euler_fin[1],
                                               euler_fin[2], "xyz" );
        EXPECT_MATRIX_NEAR( rot_start, rot_end, 1e-8 );
      }
    }
  }
}
