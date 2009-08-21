// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestCorrelator.h
#include <cxxtest/TestSuite.h>

#include <vw/Core.h>
#include <vw/Camera/PinholeModel.h>
#include <vw/Camera/CameraModel.h>
#include <vw/Stereo/StereoModel.h>
#include <vw/Math/EulerAngles.h>

using namespace vw;
using namespace vw::stereo;
using namespace vw::camera;
using namespace vw::math;

class TestStereoModel : public CxxTest::TestSuite
{
public:

  void test_pinhole_stereo_model() {

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
    //    TS_TRACE(stringify("Original Point: ") + stringify(point));
    Vector2 px1 = pin1.point_to_pixel(point);
    Vector2 px2 = pin2.point_to_pixel(point);
    //     TS_TRACE(stringify("PX 1: ") + stringify(px1));
    //     TS_TRACE(stringify("PX 2: ") + stringify(px2));
    
    StereoModel st(pin1,pin2);
    double error;
    Vector3 pt2 = st(px1, px2, error);
    
    //     TS_TRACE(stringify(pin1.camera_center(px1)));
    //     TS_TRACE(stringify(pin2.camera_center(px2)));
    //     TS_TRACE(stringify(pin1.pixel_to_vector(px1)));
    //     TS_TRACE(stringify(pin2.pixel_to_vector(px2)));
    
    //     TS_TRACE(stringify("Reconstructed Point: ") + stringify(pt2));
    
    TS_ASSERT_EQUALS(point[0], pt2[0]);
    TS_ASSERT_EQUALS(point[1], pt2[1]);
    TS_ASSERT_EQUALS(point[2], pt2[2]);
  }

  void test_adjusted_stereo_model() {

    Vector3 pos1, pos2;
    pos2 = Vector3(1,0,0);

    Matrix3x3 pose1, pose2;
    pose1.set_identity();
    pose2.set_identity();

    boost::shared_ptr<CameraModel> pin1 = boost::shared_ptr<CameraModel>(new camera::PinholeModel( pos1, pose1,
                                                                                                   1, 1, 0, 0));
    boost::shared_ptr<CameraModel> pin2 = boost::shared_ptr<CameraModel>(new camera::PinholeModel( pos2, pose2,
                                                                                                   1, 1, 0, 0));

    Vector3 point(2,0,1);
//     TS_TRACE(stringify("Original Point: ") + stringify(point));

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
//     TS_TRACE(stringify("PX 1: ") + stringify(px1));
//     TS_TRACE(stringify("PX 2: ") + stringify(px2));

    StereoModel st(adj1,adj2);
    double error;
    Vector3 pt2 = st(px1, px2, error);

//     TS_TRACE(stringify(adj1.camera_center(px1)));
//     TS_TRACE(stringify(adj2.camera_center(px2)));

//     TS_TRACE(stringify(adj1.pixel_to_vector(px1)));
//     TS_TRACE(stringify(adj2.pixel_to_vector(px2)));

//     TS_TRACE(stringify("Reconstructed Point: ") + stringify(pt2));

    TS_ASSERT_DELTA(point[0], pt2[0], 1e-6);
    TS_ASSERT_DELTA(point[1], pt2[1], 1e-6);
    TS_ASSERT_DELTA(point[2], pt2[2], 1e-6);
  }


};
