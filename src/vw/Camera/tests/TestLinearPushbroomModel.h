// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
    Quaternion<double> pose(0,0,0,1);
    Vector3 position(0,0,1);
    Vector3 velocity(1,0,0);

    // create a simplistic orbiting pushbroom camera model
    vw::camera::LinearPushbroomModel cam(10.0, // scan duration
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

//     TS_TRACE( stringify("[0,0]: ") + stringify(cam.pixel_to_vector(Vector2(0,0))));
//     TS_TRACE( stringify("       ") + stringify(cam.camera_center(Vector2(0,0))));

//     TS_TRACE( stringify("[512,0]: ") + stringify(cam.pixel_to_vector(Vector2(512,0))));
//     TS_TRACE( stringify("         ") + stringify(cam.camera_center(Vector2(512,0))));

//     TS_TRACE( stringify("[1024,0]: ") + stringify(cam.pixel_to_vector(Vector2(1024,0))));
//     TS_TRACE( stringify("          ") + stringify(cam.camera_center(Vector2(1024,0))));

//     TS_TRACE( stringify("[0,500]: ") + stringify(cam.pixel_to_vector(Vector2(0,512))));
//     TS_TRACE( stringify("         ") + stringify(cam.camera_center(Vector2(0,512))));

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
    TS_ASSERT_DELTA(camera_center[0], 5.12, 0.001);
    TS_ASSERT_EQUALS(camera_center[1], 0);
    TS_ASSERT_EQUALS(camera_center[2], 1);
    TS_ASSERT_EQUALS(pointing_vector[0], 0);
    TS_ASSERT_DELTA(pointing_vector[1], 0.981455, 0.0001);
    TS_ASSERT_DELTA(pointing_vector[2], 0.191691, 0.0001);

  }

};


