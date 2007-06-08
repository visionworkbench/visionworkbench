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

// TestMatrix.h
#include <cxxtest/TestSuite.h>
#include <vw/Math/Quaternion.h>

using namespace vw;

class TestQuaternion : public CxxTest::TestSuite
{
public:
  void test_axis_angle() {
    Vector<double,3> x_axis(1,0,0);
    Vector<double,3> y_axis(0,1,0);
    Vector<double,3> z_axis(0,0,1);

    // No rotation;  Quaternion should be (1,0,0,0)
    Quaternion<double> q1(x_axis, 0);
    TS_ASSERT_DELTA(math::norm_2(q1 - Quaternion<double>(1,0,0,0)), 0, 1e-10);

    // 90-degree rotation about X axis;  Quaternion should be (.707,.707,0,0)
    Quaternion<double> q2(x_axis, .5*M_PI);
    TS_ASSERT_DELTA(math::norm_2(q2 - Quaternion<double>(1./sqrt(2.),1./sqrt(2.),0,0)), 0, 1e-10);
    
    // 180-degree rotation about X axis;  Quaternion should be (0,1,0,0)
    Quaternion<double> q3(x_axis, M_PI);
    TS_ASSERT_DELTA(math::norm_2(q3 - Quaternion<double>(0,1,0,0)), 0, 1e-10);
    
    // 180-degree rotation about Y axis;  Quaternion should be (0,0,1,0)
    Quaternion<double> q4(y_axis, M_PI);
    TS_ASSERT_DELTA(math::norm_2(q4 - Quaternion<double>(0,0,1,0)), 0, 1e-10);
    
    // 180-degree rotation about Z axis;  Quaternion should be (0,0,0,1)
    Quaternion<double> q5(z_axis, M_PI);
    TS_ASSERT_DELTA(math::norm_2(q5 - Quaternion<double>(0,0,0,1)), 0, 1e-10);

    // Arbitrary axis
    Vector<double,3> axis(1.2, -2.3, 3.4);
    axis /= math::norm_2(axis);

    for (double angle = -4*M_PI; angle < 4*M_PI; angle += .1) {
      Quaternion<double> q(axis, angle);
      Vector<double,3> test_axis;
      double test_angle;
      q.axis_angle(test_axis, test_angle);
      //std::cout << "\nin: " << axis << "," << angle << std::endl;
      //std::cout << "quaternion: " << q << std::endl;
      //std::cout << "out: " << test_axis << "," << test_angle << std::endl;

      if (math::dot_prod(test_axis, axis) < 0) {
        test_axis = -test_axis;
        test_angle = -test_angle;
      }

      TS_ASSERT_DELTA(fmod(angle - test_angle, M_PI*2), 0, 1e-10);
      if (fabs(fmod(angle, M_PI*2)) > 1e-10) {
        TS_ASSERT_DELTA(math::norm_2(test_axis - axis), 0, 1e-10);
      }
    }
  }

}; // class TestQuaternion
