// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestCorrelator.h
#include <gtest/gtest.h>
#include <test/Helpers.h>

#include <vw/Core.h>
#include <vw/Camera/PinholeModel.h>
#include <vw/Camera/CameraModel.h>
#include <vw/Stereo/StereoModel.h>
#include <vw/Math/EulerAngles.h>

using namespace vw;
using namespace vw::stereo;
using namespace vw::camera;
using namespace vw::math;

TEST( StereoModel, PinholeStereo ) {
  Vector3 pos1, pos2;
  pos2 = Vector3(1,0,0);

  Matrix3x3 pose1, pose2;
  pose1.set_identity();
  pose2.set_identity();
  camera::PinholeModel pin1( pos1, pose1,
                             1, 1, 0, 0);

  camera::PinholeModel pin2( pos2, pose2,
                             1, 1, 0, 0);

  Vector3 point(2,0,1);
  Vector2 px1 = pin1.point_to_pixel(point);
  Vector2 px2 = pin2.point_to_pixel(point);

  StereoModel st(&pin1,&pin2);
  double error;
  Vector3 pt2 = st(px1, px2, error);

  EXPECT_VECTOR_DOUBLE_EQ( point, pt2 );
}

TEST( StereoModel, AdjustedStereo ) {

  Vector3 pos1, pos2;
  pos2 = Vector3(1,0,0);

  Matrix3x3 pose1, pose2;
  pose1.set_identity();
  pose2.set_identity();

  boost::shared_ptr<CameraModel> pin1 = boost::shared_ptr<CameraModel>(new camera::PinholeModel( pos1, pose1, 1, 1, 0, 0));
  boost::shared_ptr<CameraModel> pin2 = boost::shared_ptr<CameraModel>(new camera::PinholeModel( pos2, pose2, 1, 1, 0, 0));

  Vector3 point(2,0,1);

  camera::AdjustedCameraModel adj1(pin1);
  camera::AdjustedCameraModel adj2(pin2);

  // Set a "random" rotation
  Quaternion<double> q1 = euler_to_quaternion(M_PI/8, M_PI/12, M_PI/15, "xyz");
  Quaternion<double> q2 = euler_to_quaternion(M_PI/9, M_PI/13, M_PI/9.9, "xyz");
  adj1.set_rotation(q1);
  adj2.set_rotation(q2);

  // and a "random" translation
  adj1.set_translation(Vector3(0.2, 0.14, 0.033));
  adj2.set_translation(Vector3(0.1, 0.04, 0.123));

  Vector2 px1 = adj1.point_to_pixel(point);
  Vector2 px2 = adj2.point_to_pixel(point);

  StereoModel st(&adj1,&adj2);
  double error;
  Vector3 pt2 = st(px1, px2, error);

  EXPECT_VECTOR_NEAR( point, pt2, 1e-6 );
}


