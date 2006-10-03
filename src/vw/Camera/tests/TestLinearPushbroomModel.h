// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
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

// TestLinearPushbroomModel.h
#include <cxxtest/TestSuite.h>

#include <vw/Math/Vector.h>
#include <vw/Camera/LinearPushbroomModel.h>

#include <vector>

using namespace std;
using namespace vw;

class TestLinearPushbroomModel : public CxxTest::TestSuite
{
public:

  void testVectorToPixel()
  {
    // create a simplistic orbiting pushbroom camera model
    vw::camera::LinearPushbroomModel cam;
    
    cam.set_scan_duration(10.0);
    cam.set_number_of_lines(1000);
    cam.set_samples_per_line(1024);
    cam.set_focal_length(1.0);
    cam.set_across_scan_pixel_size(0.01);
    cam.set_along_scan_pixel_size(0.01);
    cam.set_sample_offset(-512); 
    
    Quaternion<double> quat(0,0,0,1);
    Vector3 position(0,0,1);
    Vector3 velocity(1,0,0);

    cam.set_camera_pose(quat);
    cam.set_initial_position(position);
    cam.set_velocity_vector(velocity);

    std::cout << "\n\n";

    std::cout << "[0,0]: " << cam.pixel_to_vector(Vector2(0,0)) << "\n";
    std::cout << "       " << cam.camera_center(Vector2(0,0)) << "\n\n";

    std::cout << "[512,0]: " << cam.pixel_to_vector(Vector2(512,0)) << "\n";
    std::cout << "       " << cam.camera_center(Vector2(512,0)) << "\n\n";

    std::cout << "[1024,0]: " << cam.pixel_to_vector(Vector2(1024,0)) << "\n";
    std::cout << "       " << cam.camera_center(Vector2(1024,0)) << "\n\n";

    std::cout << "[0,500]: " << cam.pixel_to_vector(Vector2(0,512)) << "\n";
    std::cout << "         " << cam.camera_center(Vector2(0,512)) << "\n\n";
    
    Vector3 camera_center, pointing_vector;
    pointing_vector = cam.pixel_to_vector(Vector2(0,0));
    camera_center = cam.camera_center(Vector2(0,0));
    TS_ASSERT_EQUALS(camera_center[0], 0);
    TS_ASSERT_EQUALS(camera_center[1], 0);
    TS_ASSERT_EQUALS(camera_center[2], 1);
    TS_ASSERT_EQUALS(pointing_vector[0], 0);
    TS_ASSERT_DELTA(pointing_vector[1], 0.981455, 0.0001);
    TS_ASSERT_DELTA(pointing_vector[2], 0.191691, 0.0001);

    pointing_vector = cam.pixel_to_vector(Vector2(512,0));
    camera_center = cam.camera_center(Vector2(512,0));
    TS_ASSERT_EQUALS(camera_center[0], 0);
    TS_ASSERT_EQUALS(camera_center[1], 0);
    TS_ASSERT_EQUALS(camera_center[2], 1);
    TS_ASSERT_EQUALS(pointing_vector[0], 0);
    TS_ASSERT_EQUALS(pointing_vector[1], 0);
    TS_ASSERT_EQUALS(pointing_vector[2], 1);

    pointing_vector = cam.pixel_to_vector(Vector2(0,512));
    camera_center = cam.camera_center(Vector2(0,512));
    TS_ASSERT_EQUALS(camera_center[0], 5.12);
    TS_ASSERT_EQUALS(camera_center[1], 0);
    TS_ASSERT_EQUALS(camera_center[2], 1);
    TS_ASSERT_EQUALS(pointing_vector[0], 0);
    TS_ASSERT_DELTA(pointing_vector[1], 0.981455, 0.0001);
    TS_ASSERT_DELTA(pointing_vector[2], 0.191691, 0.0001);
    
  }

};


