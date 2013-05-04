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


#include <gtest/gtest_VW.h>

#include <vw/Math/Vector.h>
#include <vw/Camera/LinearPushbroomModel.h>
#include <test/Helpers.h>

using namespace vw;
using namespace camera;

TEST( LinearPushbroom, PixelToVector ) {

  Quaternion<double> pose(0,0,0,1);
  Vector3 position(0,0,1);
  Vector3 velocity(1,0,0);

  // create a simplistic orbiting pushbroom camera model
  LinearPushbroomModel cam(10.0, // scan duration
                           1000, // number of lines
                           1024, // samples per line
                           -512, // sample offset
                           1.0,  // focal length
                           0.01, // along_scan_pixel_size
                           0.01, // across_scan_pixel_size
                           Vector3(0,0,1),  // looakat vector
                           Vector3(0,1,0),  // horizontal pixel vector
                           pose,
                           position,
                           velocity);
  EXPECT_STREQ( "LinearPushbroom", cam.type().c_str() );

  Vector3 camera_center, pointing_vector;
  pointing_vector = cam.pixel_to_vector(Vector2(0,0));
  camera_center = cam.camera_center(Vector2(0,0));
  EXPECT_VECTOR_DOUBLE_EQ( camera_center,
                           Vector3(0,0,1) );
  EXPECT_EQ( pointing_vector[0], 0 );
  EXPECT_VECTOR_NEAR( subvector(pointing_vector,1,2),
                      Vector2( 0.981455, 0.191691), 1e-4 );

  pointing_vector = cam.pixel_to_vector(Vector2(512,0));
  camera_center = cam.camera_center(Vector2(512,0));
  EXPECT_VECTOR_DOUBLE_EQ( camera_center,
                           Vector3(0,0,1) );
  EXPECT_VECTOR_DOUBLE_EQ( pointing_vector,
                           Vector3(0,0,1) );

  pointing_vector = cam.pixel_to_vector(Vector2(0,512));
  camera_center = cam.camera_center(Vector2(0,512));
  EXPECT_VECTOR_NEAR( camera_center,
                      Vector3(5.12,0,1),
                      1e-3 );
  EXPECT_EQ( pointing_vector[0], 0 );
  EXPECT_VECTOR_NEAR( subvector(pointing_vector,1,2),
                      Vector2( 0.981455, 0.191691), 1e-4 );

  // check enforcement that pose returns the rotation from camera
  // frame to world frame.
  Vector2 center_pixel(512,500);
  Quaternion<double> center_pose = cam.camera_pose(center_pixel);
  double angle_from_z =
    acos(dot_prod(Vector3(0,0,1),inverse(center_pose).rotate(cam.pixel_to_vector(center_pixel))));
  EXPECT_LT( angle_from_z, 0.5 );
}
