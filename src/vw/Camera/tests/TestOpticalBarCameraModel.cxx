// __BEGIN_LICENSE__
//  Copyright (c) 2009-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NGT platform is licensed under the Apache License, Version 2.0 (the
//  "License"); you may not use this file except in compliance with the
//  License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#include <vw/Camera/OpticalBarModel.h>

#include <boost/scoped_ptr.hpp>
#include <test/Helpers.h>

using namespace vw;
using namespace vw::test;

// Build an optical bar camera from parameters, round-trip it through a temp
// file, and verify project-out-and-back self-consistency.
TEST(OpticalBarModel, CreateCamera) {

  Vector3   gcc(-2470899.9105873918, 5542443.1851118281, 2502370.2339997198);
  Matrix3x3 rot_mat(-0.057502980758053046, -0.74594045857037883,  0.66352561327483939,
                    -0.16417103320903326,  -0.64851223071215924, -0.74328982131589649,
                     0.98475442576259453,  -0.15167306578477252, -0.085170429471915887);
  vw::Quat q(rot_mat);
  Vector3  initial_position(gcc);
  Vector3  initial_pose = q.axis_angle();
  
  double   pixel_size =  0.000112;
  Vector2i image_size(3910, 2290);
  Vector2  center_loc_pixels = Vector2(1957.2358579744291, 1112.512227705716);
  double   focal_length       = 1.9600000381469727;
  double   scan_time = 0.5;
  //double   scan_angle_radians = 1.2915399999999999;
  //double   scan_rate_radians  = 3.0125446909089169;
  bool     scan_left_to_right = true;
  double   forward_tilt_radians = 0;
  double   velocity = 5253.0854975797192;
  bool     use_motion_comp = true;


  OpticalBarModel* raw_ptr = new OpticalBarModel(image_size, center_loc_pixels, pixel_size,
                                                 focal_length, scan_time, scan_left_to_right,
                                                 forward_tilt_radians,
                                                 initial_position, initial_pose,
                                                 velocity, use_motion_comp);

  // Basic file I/O
  OpticalBarModel cpy;
  UnlinkName tsaiFile("test.tsai");
  raw_ptr->write(tsaiFile);
  cpy.read (tsaiFile);
  EXPECT_NO_THROW(raw_ptr->write(tsaiFile));
  EXPECT_NO_THROW(cpy.read (tsaiFile));

  EXPECT_VECTOR_EQ(raw_ptr->get_image_size(), cpy.get_image_size());

  vw::CamPtr cam1;
  cam1 = vw::CamPtr(raw_ptr);

  ASSERT_TRUE(cam1.get() != 0);


  // A more accurate test is just to project out and back into the same camera
  for (size_t i = 0; i < 3000; i += 200) {
    for (size_t j = 0; j < 2400; j += 200) {
      Vector2 pixel(i,j);
      Vector3 center = cam1->camera_center(pixel);
      EXPECT_VECTOR_NEAR(pixel,
                         cam1->point_to_pixel(center + 2e4*cam1->pixel_to_vector(pixel)),
                         1e-2 /*pixels*/);
    }
  }
}
