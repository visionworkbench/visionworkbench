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


// TestCorrelator.h
#include <gtest/gtest_VW.h>
#include <test/Helpers.h>

#include <vw/Image/ImageView.h>
#include <vw/Math/EulerAngles.h>
#include <vw/Camera/PinholeModel.h>
#include <vw/Camera/CameraModel.h>
#include <vw/Stereo/StereoModel.h>
#include <vw/Stereo/StereoView.h>

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

  {
    StereoModel st(&pin1,&pin2);
    double error;
    Vector3 pt2 = st(px1, px2, error);
    EXPECT_VECTOR_NEAR( point, pt2, 1e-6 );
  }
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

  {
    StereoModel st(&adj1,&adj2);
    double error;
    Vector3 pt2 = st(px1, px2, error);

    EXPECT_VECTOR_NEAR( point, pt2, 1e-6 );
  }
}

TEST(StereoView, PixelMaskVec2) {
  Vector3 pos1, pos2;
  pos2 = Vector3(1,0,0);

  Matrix3x3 pose1, pose2;
  pose1.set_identity();
  pose2.set_identity();

  camera::PinholeModel pin1( Vector3(),
                             math::identity_matrix<3>(),
                             1, 1, 1, 0);

  camera::PinholeModel pin2( Vector3(1,0,0),
                             math::identity_matrix<3>(),
                             1, 1, 1, 0);

  ImageView<PixelMask<Vector2f> > disparity(3,1);
  disparity(1,0) = PixelMask<Vector2f>( Vector2f(-1,0) );
  disparity(2,0) = PixelMask<Vector2f>( Vector2f(-1.3,0) );

  ImageView<Vector3> pc = stereo_triangulate( disparity, &pin1, &pin2 );
  EXPECT_VECTOR_DOUBLE_EQ( pc(0,0), Vector3() );
  EXPECT_VECTOR_NEAR( pc(1,0), Vector3(0,0,1), 1e-2 );
  EXPECT_VECTOR_NEAR( pc(2,0), Vector3(0.769,0,0.769), 1e-2 );
}

TEST( StereoView, Vec2 ) {
  Vector3 pos1, pos2;
  pos2 = Vector3(1,0,0);

  Matrix3x3 pose1, pose2;
  pose1.set_identity();
  pose2.set_identity();

  camera::PinholeModel pin1( Vector3(),
                             math::identity_matrix<3>(),
                             1, 1, 1, 0);

  camera::PinholeModel pin2( Vector3(1,0,0),
                             math::identity_matrix<3>(),
                             1, 1, 1, 0);

  ImageView<Vector2f> disparity(3,1);
  disparity(0,0) = Vector2f(-1.5,0);
  disparity(1,0) = Vector2f(-1,0);
  disparity(2,0) = Vector2f(-1.3,0);

  ImageView<Vector3> pc = stereo_triangulate( disparity, &pin1, &pin2 );
  EXPECT_VECTOR_NEAR( pc(0,0), Vector3(-0.666,0,0.666), 1e-2 );
  EXPECT_VECTOR_NEAR( pc(1,0), Vector3(0,0,1), 1e-2 );
  EXPECT_VECTOR_NEAR( pc(2,0), Vector3(0.769,0,0.769), 1e-2 );
}

TEST( StereoView, PixelMaskFloat ) {
  Vector3 pos1, pos2;
  pos2 = Vector3(1,0,0);

  Matrix3x3 pose1, pose2;
  pose1.set_identity();
  pose2.set_identity();

  camera::PinholeModel pin1( Vector3(),
                             math::identity_matrix<3>(),
                             1, 1, 1, 0);

  camera::PinholeModel pin2( Vector3(1,0,0),
                             math::identity_matrix<3>(),
                             1, 1, 1, 0);

  ImageView<PixelMask<float> > disparity(3,1);
  disparity(1,0) = PixelMask<float>(-1);
  disparity(2,0) = PixelMask<float>(-1.3);

  ImageView<Vector3> pc = stereo_triangulate( disparity, &pin1, &pin2 );
  EXPECT_VECTOR_DOUBLE_EQ( pc(0,0), Vector3() );
  EXPECT_VECTOR_NEAR( pc(1,0), Vector3(0,0,1), 1e-2 );
  EXPECT_VECTOR_NEAR( pc(2,0), Vector3(0.769,0,0.769), 1e-2 );
}

TEST( StereoView, Float ) {
  Vector3 pos1, pos2;
  pos2 = Vector3(1,0,0);

  Matrix3x3 pose1, pose2;
  pose1.set_identity();
  pose2.set_identity();

  camera::PinholeModel pin1( Vector3(),
                             math::identity_matrix<3>(),
                             1, 1, 1, 0);

  camera::PinholeModel pin2( Vector3(1,0,0),
                             math::identity_matrix<3>(),
                             1, 1, 1, 0);

  ImageView<float> disparity(3,1);
  disparity(0,0) = -1.5;
  disparity(1,0) = -1;
  disparity(2,0) = -1.3;

  ImageView<Vector3> pc = stereo_triangulate( disparity, &pin1, &pin2 );
  EXPECT_VECTOR_NEAR( pc(0,0), Vector3(-0.666,0,0.666), 1e-2 );
  EXPECT_VECTOR_NEAR( pc(1,0), Vector3(0,0,1), 1e-2 );
  EXPECT_VECTOR_NEAR( pc(2,0), Vector3(0.769,0,0.769), 1e-2 );
}
