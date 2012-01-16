// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestPinholeModel.h
#include <gtest/gtest_VW.h>
#include <test/Helpers.h>
#include <vw/Math/EulerAngles.h>

#include <vw/Camera/PinholeModel.h>

using namespace vw;
using namespace vw::camera;
using namespace vw::test;

TEST( AdjustedCameraModel, StandardConstruct ) {
  Matrix<double,3,3> pose;
  pose.set_identity();

  // Create an imaginary 1000x1000 pixel imager
  boost::shared_ptr<CameraModel> pinhole(
      new PinholeModel( Vector3(0,0,0), // camera center
                        pose,           // camera pose
                        500,500,        // fx, fy
                        500,500,        // cx, cy
                        NullLensDistortion()) );

  AdjustedCameraModel adjcam( pinhole );
  EXPECT_STREQ( "Adjusted", adjcam.type().c_str() );

  // No mods, so should have equal results
  EXPECT_VECTOR_NEAR( pinhole->camera_center(Vector2()),
                      adjcam.camera_center(Vector2()), 1e-5 );
  EXPECT_VECTOR_NEAR( pinhole->pixel_to_vector( Vector2(750,750) ),
                      adjcam.pixel_to_vector(   Vector2(750,750) ),
                      1e-5 );
  EXPECT_VECTOR_NEAR( pinhole->pixel_to_vector( Vector2(55,677) ),
                      adjcam.pixel_to_vector(   Vector2(55,677) ),
                      1e-5 );
  EXPECT_NEAR( pinhole->camera_pose(Vector2())[0],
               adjcam.camera_pose(Vector2())[0], 1e-5 );
  EXPECT_NEAR( pinhole->camera_pose(Vector2())[1],
               adjcam.camera_pose(Vector2())[1], 1e-5 );
  EXPECT_NEAR( pinhole->camera_pose(Vector2())[2],
               adjcam.camera_pose(Vector2())[2], 1e-5 );
}

TEST( AdjustedCameraModel, AdjustedConstruct ) {

  Matrix<double,3,3> pose = math::euler_to_rotation_matrix(1.3,2.0,-.7,"xyz");

  // Create an imaginary 1000x1000 pixel imager
  boost::shared_ptr<CameraModel> pinhole(
      new PinholeModel( Vector3(0,0,0), // camera center
                        pose,           // camera pose
                        500,500,        // fx, fy
                        500,500,        // cx, cy
                        NullLensDistortion()) );

  AdjustedCameraModel adjcam( pinhole, Vector3(1,0,0), Quat(1,0,0,0) );
  EXPECT_VECTOR_NEAR( Vector3(1,0,0), adjcam.camera_center(Vector2()) -
                      pinhole->camera_center(Vector2()), 1e-5 );
  EXPECT_VECTOR_NEAR( pinhole->pixel_to_vector( Vector2(750,750) ),
                      adjcam.pixel_to_vector( Vector2(750,750) ), 1e-5 );
  EXPECT_VECTOR_DOUBLE_EQ( Vector3(1,0,0), adjcam.translation() );

  UnlinkName adjustment("adjust.txt");
  adjcam.write( adjustment );

  // Read back in adjustment
  AdjustedCameraModel adjcam2( pinhole );
  adjcam2.read( adjustment );
  EXPECT_VECTOR_NEAR( adjcam.translation(),
                      adjcam2.translation(), 1e-5 );
  EXPECT_MATRIX_NEAR( adjcam.rotation_matrix(),
                      adjcam2.rotation_matrix(), 1e-5 );
  EXPECT_VECTOR_NEAR( adjcam.axis_angle_rotation(),
                      adjcam2.axis_angle_rotation(), 1e-5 );

  // check enforcement that pose returns the rotation from camera
  // frame to world frame.
  Vector2 center_pixel(500,500);
  Quaternion<double> center_pose =
    adjcam2.camera_pose(center_pixel);
  double angle_from_z =
    acos(dot_prod(Vector3(0,0,1),inverse(center_pose).rotate(adjcam2.pixel_to_vector(center_pixel))));
  EXPECT_LT( angle_from_z, 0.5 );
}
