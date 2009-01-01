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
      //TS_TRACE(stringify("in: "        ) + stringify(axis) + "," + stringify(angle));
      //TS_TRACE(stringify("quaternion: ") + stringify(q));
      //TS_TRACE(stringify("out: "       ) + stringify(test_axis) + "," + stringify(test_angle))

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

  void test_rotation()
  {
    // Create a test vector to rotate
    Vector3 test_vec(0.1, 0.84, 0.23);

    double theta = 32*M_PI/180;  // 32 degree rotation
    Matrix3x3 rotation = math::identity_matrix<3>();
    rotation(0,0) = cos(theta);
    rotation(0,1) = -sin(theta);
    rotation(1,0) = sin(theta);
    rotation(1,1) = cos(theta);

    Quaternion<double> quat(rotation);
    Matrix3x3 rotation2 = quat.rotation_matrix();

    Vector3 sol1 = rotation * test_vec;
    Vector3 sol2 = rotation2 * test_vec;
    Vector3 sol3 = quat.rotate(test_vec);
    Vector3 sol4 = -quat.rotate(-test_vec);

    TS_ASSERT_DELTA(sol1(0),sol2(0),1e-9);
    TS_ASSERT_DELTA(sol1(1),sol2(1),1e-9);
    TS_ASSERT_DELTA(sol1(2),sol2(2),1e-9);
    TS_ASSERT_DELTA(sol1(0),sol3(0),1e-9);
    TS_ASSERT_DELTA(sol1(1),sol3(1),1e-9);
    TS_ASSERT_DELTA(sol1(2),sol3(2),1e-9);
    TS_ASSERT_DELTA(sol1(0),sol4(0),1e-9);
    TS_ASSERT_DELTA(sol1(1),sol4(1),1e-9);
    TS_ASSERT_DELTA(sol1(2),sol4(2),1e-9);
  }


  void test_rotationmatrix() {
    vw::Matrix3x3 I;
    I.set_identity();

    vw::Matrix3x3 T, dQ;


    // test ww
    double R_data_ww[9] = { 0.36, 0.48, -0.8, -0.8, 0.6, 0, 0.48, 0.64, 0.6 };
    vw::Matrix3x3 R_ww(R_data_ww);

    // check if R*transpose(R) is identity
    T = R_ww * transpose(R_ww) - I;
    for (Matrix3x3::const_iterator it = T.begin(); it != T.end(); ++it)
      TS_ASSERT_DELTA(*it, 0, 1e-14);

    // converting the matrix to quaternion and back should not change it
    dQ = R_ww - vw::Quaternion<double>(R_ww).rotation_matrix();
    for (Matrix3x3::const_iterator it = dQ.begin(); it != dQ.end(); ++it)
      TS_ASSERT_DELTA(*it, 0, 1e-14);


    // test xx
    double R_data_xx[9] = { 0.11276615807984,   -0.0161365698014631,  0.993490515660293,
			    -0.633815522454304, -0.771198883306014,   0.0594151990953122,
			    0.765220018744864, -0.636389733970151,  -0.097192743605217,
    };
    vw::Matrix3x3 R_xx(R_data_xx);

    // check if is identity
    T = R_xx * transpose(R_xx) - I;
    for (Matrix3x3::const_iterator it = T.begin(); it != T.end(); ++it)
      TS_ASSERT_DELTA(*it, 0, 1e-14);

    // converting the matrix to quaternion and back should not change it
    dQ = R_xx - vw::Quaternion<double>(R_xx).rotation_matrix();
    for (Matrix3x3::const_iterator it = dQ.begin(); it != dQ.end(); ++it)
      TS_ASSERT_DELTA(*it, 0, 1e-14);


    // test yy
    double R_data_yy[9] = { -0.771198883306014, -0.633815522454304,  0.0594151990953122,
			    -0.0161365698014631, 0.11276615807984,   0.993490515660293,
			    -0.636389733970151,  0.765220018744864, -0.097192743605217,
    };
    vw::Matrix3x3 R_yy(R_data_yy);

    // check if is identity
    T = R_yy * transpose(R_yy) - I;
    for (Matrix3x3::const_iterator it = T.begin(); it != T.end(); ++it)
      TS_ASSERT_DELTA(*it, 0, 1e-14);

    // converting the matrix to quaternion and back should not change it
    dQ = R_yy - vw::Quaternion<double>(R_yy).rotation_matrix();
    for (Matrix3x3::const_iterator it = dQ.begin(); it != dQ.end(); ++it)
      TS_ASSERT_DELTA(*it, 0, 1e-14);


    // test zz
    double R_data_zz[9] = { -0.771198883306014,  0.0594151990953122, -0.633815522454304,
			    -0.636389733970151, -0.097192743605217,   0.765220018744864,
			    -0.0161365698014631, 0.993490515660293,   0.11276615807984 
    };
    vw::Matrix3x3 R_zz(R_data_zz);

    T = R_zz * transpose(R_zz) - I;
    for (Matrix3x3::const_iterator it = T.begin(); it != T.end(); ++it)
      TS_ASSERT_DELTA(*it, 0, 1e-14);

    // converting the matrix to quaternion and back should not change it
    dQ = R_zz - vw::Quaternion<double>(R_zz).rotation_matrix();
    for (Matrix3x3::const_iterator it = dQ.begin(); it != dQ.end(); ++it)
      TS_ASSERT_DELTA(*it, 0, 1e-14);
  }
}; // class TestQuaternion
