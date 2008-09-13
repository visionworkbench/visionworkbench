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

#include <math.h>

#include <cxxtest/TestSuite.h>
#include <vw/Math/EulerAngles.h>

using namespace vw;
using namespace math;

#include <iostream>

class TestEuler : public CxxTest::TestSuite
{
public:

  void test_euler_angle_xyz() {
    std::cout << "\nTesting XYZ Euler Combination\n";
    for (unsigned p = 0; p < 360; p+=60){
      for (unsigned o = 0; o < 360; o+=60){
	for (unsigned k = 0; k < 360; k+=60){
	  double phi = p * M_PI / 180;
	  double omega = o * M_PI / 180;
	  double kappa = k * M_PI / 180;

	  vw::Matrix<double,3,3> rot_matrix = euler_to_rotation_matrix( phi, omega, kappa, "xyz");
	  Vector3 euler = rotation_matrix_to_euler_xyz( rot_matrix );

	  // PASS 1
	  vw::Matrix<double,3,3> rot_matrix2 = euler_to_rotation_matrix( euler[0], euler[1], euler[2], "xyz");
	  rot_matrix2 = rot_matrix2 - rot_matrix;
	  for (unsigned u = 0; u < 3; ++u){
	    for (unsigned v = 0; v < 3; ++v){
	      TS_ASSERT_DELTA( rot_matrix2(u,v), 0, 1e-4);
	    }
	  }
	}
      }
    }
    std::cout << "Finished.\n";
  }

  void test_euler_angle_yxz() {
    std::cout << "\nTesting YXZ Euler Combination\n";
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
	  rot_matrix2 = rot_matrix2 - rot_matrix;
	  for (unsigned u = 0; u < 3; ++u){
	    for (unsigned v = 0; v < 3; ++v){
	      TS_ASSERT_DELTA( rot_matrix2(u,v), 0, 1e-4);
	    }
	  }
	}
      }
    }
    std::cout << "Finished.\n";
  }

  void test_euler_angle_zxy() {
    std::cout << "\nTesting ZXY Euler Combination\n";
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
	  rot_matrix2 = rot_matrix2 - rot_matrix;
	  for (unsigned u = 0; u < 3; ++u){
	    for (unsigned v = 0; v < 3; ++v){
	      TS_ASSERT_DELTA( rot_matrix2(u,v), 0, 1e-4);
	    }
	  }
	}
      }
    }
    std::cout << "Finished.\n";
  }
}; // class TestEuler
