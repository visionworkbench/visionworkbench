// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestLinearPushbroomModel.h
#include <cxxtest/TestSuite.h>

#include <vw/Camera/PinholeModel.h>
#include <vw/Math/Vector.h>

#include <vw/FileIO.h>
#include <vw/Camera/CameraTransform.h>
#include <vw/Math/EulerAngles.h>

using namespace vw;
using namespace vw::camera;

class TestPinholeModel : public CxxTest::TestSuite
{
public:

  void test_basic_pinhole_model()
  {
    Matrix<double,3,3> pose;
    pose.set_identity();

    // Create an imaginary 1000x1000 pixel imager
    PinholeModel pinhole( Vector3(0,0,0), // camera center
                          pose,           // camera pose
                          500,500,        // fx, fy
                          500,500,
                          NullLensDistortion());       // cx, cy

    // TS_TRACE(stringify(pinhole));
    TS_ASSERT_EQUALS(pinhole.point_to_pixel(Vector3(0,0,10)),Vector2(500,500));
    TS_ASSERT_EQUALS(pinhole.point_to_pixel(Vector3(-10,0,10)),Vector2(0,500));
    TS_ASSERT_EQUALS(pinhole.point_to_pixel(Vector3(10,0,10)),Vector2(1000,500));
    TS_ASSERT_EQUALS(pinhole.point_to_pixel(Vector3(0,-10,10)),Vector2(500,1000));
    TS_ASSERT_EQUALS(pinhole.point_to_pixel(Vector3(0,10,10)),Vector2(500,0));

  }


  void test_coordinate_frames()
  {
    Matrix<double,3,3> pose;
    pose.set_identity();

    // Create an imaginary 1000x1000 pixel imager, where the camera
    // coordinate system is mapped as follows:
    //
    // +u : along the camera +Y axis
    // +v : along the camera +X axis
    // +w : along the camera -Z axis
    PinholeModel pinhole( Vector3(0,0,0), // camera center
                          pose,           // camera pose
                          500,500,        // fx, fy
                          500,500,
                          Vector3(0, 1, 0),
                          Vector3(1, 0, 0),
                          Vector3(0, 0, -1),

                          NullLensDistortion());       // cx, cy

    // TS_TRACE(stringify(pinhole));
    TS_ASSERT_EQUALS(pinhole.point_to_pixel(Vector3(0,0,-10)),Vector2(500,500));
    TS_ASSERT_EQUALS(pinhole.point_to_pixel(Vector3(-10,0,-10)),Vector2(500,0));
    TS_ASSERT_EQUALS(pinhole.point_to_pixel(Vector3(10,0,-10)),Vector2(500,1000));
    TS_ASSERT_EQUALS(pinhole.point_to_pixel(Vector3(0,-10,-10)),Vector2(0,500));
    TS_ASSERT_EQUALS(pinhole.point_to_pixel(Vector3(0,10,-10)),Vector2(1000,500));

    TS_TRACE(stringify(pinhole.point_to_pixel(Vector3(-10,0,-10))));
    TS_TRACE(stringify(pinhole.point_to_pixel(Vector3(10,0,-10))));
    TS_TRACE(stringify(pinhole.point_to_pixel(Vector3(0,-10,-10))));
    TS_TRACE(stringify(pinhole.point_to_pixel(Vector3(0,10,-10))));
  }

  void test_pixel_to_vector()
  {
    Matrix<double,3,3> pose;
    pose.set_identity();

    // Create an imaginary 1000x1000 pixel imager
    PinholeModel pinhole( Vector3(0,0,0), // camera center
                          pose,           // camera pose
                          500,500,        // fx, fy
                          500,500,
                          NullLensDistortion());       // cx, cy


    PinholeModel pinhole2( Vector3(10,10,10), // camera center
                          pose,           // camera pose
                          500,500,        // fx, fy
                          500,500,
                          NullLensDistortion());       // cx, cy

    Matrix<double,3,3> rot = vw::math::euler_to_quaternion(1.15, 0.0, -1.57, "xyz").rotation_matrix();
    PinholeModel pinhole3(Vector3(-0.329, 0.065, -0.82),
                         rot,
                         605.320556640625,
                         606.3638305664062,
                         518.89208984375,
                         387.5555114746094);

    PinholeModel pinhole4(Vector3(-0.329, 0.065, -0.82),
                          rot,
                          605.320556640625,
                          606.3638305664062,
                          518.89208984375,
                          387.5555114746094,
                          TsaiLensDistortion(Vector4(-0.2796604335308075,
                                                     0.1031486615538597,
                                                     -0.0007824968779459596,
                                                     0.0009675505571067333)));

    Vector2 result1 = pinhole.point_to_pixel(pinhole.pixel_to_vector(Vector2(0,0))+pinhole.camera_center(Vector2(0,0)));
    Vector2 result2 = pinhole2.point_to_pixel(pinhole2.pixel_to_vector(Vector2(0,0))+pinhole2.camera_center(Vector2(0,0)));
    Vector2 result3 = pinhole3.point_to_pixel(pinhole3.pixel_to_vector(Vector2(0,0))+pinhole3.camera_center(Vector2(0,0)));
#if defined(VW_HAVE_PKG_LAPACK) && VW_HAVE_PKG_LAPACK==1
    Vector2 result4 = pinhole4.point_to_pixel(pinhole4.pixel_to_vector(Vector2(0,0))+pinhole4.camera_center(Vector2(0,0)));
#endif
    TS_ASSERT_EQUALS(result1, Vector2(0,0));
    TS_ASSERT_EQUALS(result2, Vector2(0,0));
    TS_ASSERT_DELTA(result3[0], Vector2(0,0)[0], 1e-8);
    TS_ASSERT_DELTA(result3[1], Vector2(0,0)[1], 1e-8);
#if defined(VW_HAVE_PKG_LAPACK) && VW_HAVE_PKG_LAPACK==1
    TS_ASSERT_DELTA(result4[0], Vector4(0,0)[0], 1e-3);
    TS_ASSERT_DELTA(result4[1], Vector4(0,0)[1], 1e-3);
#endif
  }

  void test_tsai_distortion()
  {
    // Create an imaginary 1000x1000 pixel imager
    PinholeModel pinhole( Vector3(0,0,0),                 // camera center
                          math::identity_matrix<3>(),     // camera pose
                          500,500,                        // fx, fy
                          500,500,                        // cx, cy
                          TsaiLensDistortion(Vector4(-0.2805362343788147,
                                                     0.1062035113573074,
                                                     -0.0001422458299202845,
                                                     0.00116333004552871)));
    const LensDistortion* distortion = pinhole.lens_distortion();

#if defined(VW_HAVE_PKG_LAPACK) && VW_HAVE_PKG_LAPACK==1
    Vector2 distorted_pix = distortion->distorted_coordinates(pinhole, Vector2(200,200));
    Vector2 undistorted_pix = distortion->undistorted_coordinates(pinhole, distorted_pix);

    TS_ASSERT_DELTA(distorted_pix[0], 244.865, 0.1);
    TS_ASSERT_DELTA(distorted_pix[1], 244.395, 0.1);
    TS_ASSERT_DELTA(undistorted_pix[0], 200, 0.1);
    TS_ASSERT_DELTA(undistorted_pix[1], 200, 0.1);
#endif
  }

  void test_scale_pinhole()
  {
    Matrix<double,3,3> rot = vw::math::euler_to_quaternion(1.15, 0.0, -1.57, "xyz").rotation_matrix();
    PinholeModel pinhole4(Vector3(-0.329, 0.065, -0.82),
                          rot,
                          605.320556640625,
                          606.3638305664062,
                          518.89208984375,
                          387.5555114746094,
                          TsaiLensDistortion(Vector4(-0.2796604335308075,
                                                     0.1031486615538597,
                                                     -0.0007824968779459596,
                                                     0.0009675505571067333)));
    PinholeModel scaled = scale_camera(pinhole4, .1);

    Vector3 point = Vector3(2,-1,1) +5*pinhole4.pixel_to_vector(Vector2(500,500))+pinhole4.camera_center();

    Vector2 o_return = pinhole4.point_to_pixel(point);
    Vector2 s_return = scaled.point_to_pixel(point);

    TS_ASSERT_DELTA(o_return[0],s_return[0]*10,5); // Lens distortion doesn't
    TS_ASSERT_DELTA(o_return[1],s_return[1]*10,5); // seem to be too accurate
  }

};
