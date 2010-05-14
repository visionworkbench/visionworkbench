// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestPinholeModel.h
#include <gtest/gtest.h>
#include <test/Helpers.h>

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

  Matrix<double,3,3> pose;
  pose.set_identity();

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
}
