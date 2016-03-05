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
#include <vw/Camera/LinescanModel.h>
#include <vw/Camera/Extrinsics.h>
#include <test/Helpers.h>

using namespace vw;
using namespace camera;

// TODO: Is there a good test we can do here?
//       - Not as important since we are not using templates anymore.

/*
// A very simple Linescan model for testing purposes
class LinearPushbroomModel {
  public:

    LinearPushbroomModel(LinearPositionInterpolation const& position,
		                     LinearPositionInterpolation const& velocity,
	                       ConstantPoseInterpolation   const& pose,
	                       LinearTimeInterpolation     const& time,
	                       vw::Vector2i  const& image_size,
	                       vw::Vector2   const& detector_origin,
	                       double        const  focal_length
		    ) : vw::camera::LinescanModel(image_size, true), // Always correct velocity aberration
		        m_position_func(position), m_velocity_func(velocity),
            m_pose_func(pose),         m_time_func(time),
            m_detector_origin(detector_origin),
            m_focal_length(focal_length) {} 
    virtual ~LinearPushbroomModel() {}
    virtual std::string type() const { return "LinearPushbroomModel"; }
    
    

}; // End class LinearPushbroomModel
*/
TEST( LinearPushbroom, PixelToVector ) {
  EXPECT_TRUE(true);

/*
  typedef LinescanModel<LinearPositionInterpolation, // Position
                        LinearPositionInterpolation, // Velocity
                        ConstantPoseInterpolation,   // Pose
                        LinearTimeInterpolation // Time
                       > LinearPushbroomModel;

  






  Quaternion<double> pose(0,0,0,1);
  Vector3 position(0,0,1);
  Vector3 velocity(1,0,0);

  // create a simplistic orbiting pushbroom camera model
  LinearPushbroomModel cam(LinearPositionInterpolation(position, velocity),       // Linear position
                           LinearPositionInterpolation(velocity, Vector3(0,0,0)), // Constant velocity
                           ConstantPoseInterpolation(pose),   // Pose
                           LinearTimeInterpolation(0, 0.01),  // Start time, time per line
                           vw::Vector2i(1024, 1000), // Num samples, num lines
                           vw::Vector2(-512.0, 0),   // Detector origin in pixels
                           100.0, // Focal length in pixels
                           false  // Correct velocity aberration
                          );

  const double EPS = 1e-4;
  Vector3 camera_center, pointing_vector;
  pointing_vector = cam.pixel_to_vector(Vector2(0,0));
  camera_center   = cam.camera_center  (Vector2(0,0));
  EXPECT_VECTOR_DOUBLE_EQ( camera_center,
                           Vector3(0,0,1) );
  EXPECT_VECTOR_NEAR( pointing_vector, Vector3(0.981455, 0.0, 0.191691), EPS );

  pointing_vector = cam.pixel_to_vector(Vector2(512,0));
  camera_center   = cam.camera_center  (Vector2(512,0));
  EXPECT_VECTOR_DOUBLE_EQ( camera_center,   Vector3(0,0,1) );
  EXPECT_VECTOR_DOUBLE_EQ( pointing_vector, Vector3(0,0,1) );

  pointing_vector = cam.pixel_to_vector(Vector2(1024,512));
  camera_center   = cam.camera_center  (Vector2(1024,512));
  EXPECT_VECTOR_NEAR( camera_center, Vector3(5.12,0,1), EPS );
  EXPECT_VECTOR_NEAR( pointing_vector, Vector3(-0.981455, 0.0, 0.191691), EPS );

  // check enforcement that pose returns the rotation from camera
  // frame to world frame.
  Vector2 center_pixel(512,500);
  Quaternion<double> center_pose = cam.camera_pose(center_pixel);
  double angle_from_z =
    acos(dot_prod(Vector3(0,0,1),inverse(center_pose).rotate(cam.pixel_to_vector(center_pixel))));
  EXPECT_LT( angle_from_z, 0.5 );
  */
}

