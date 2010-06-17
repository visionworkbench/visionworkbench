// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestPinholeModel.h
#include <gtest/gtest.h>

#include <boost/random.hpp>

#include <vw/Math/Vector.h>
#include <vw/Math/EulerAngles.h>
#include <vw/Math/LinearAlgebra.h>
#include <test/Helpers.h>

#include <vw/FileIO.h>

#include <vw/Camera/PinholeModel.h>
#include <vw/Camera/CameraTransform.h>
#include <vw/Camera/CameraGeometry.h>

using namespace vw;
using namespace vw::camera;
using namespace vw::test;

TEST( PinholeModel, StandardConstruct ) {
  Matrix<double,3,3> pose;
  pose.set_identity();

  // Create an imaginary 1000x1000 pixel imager
  PinholeModel pinhole( Vector3(0,0,0), // camera center
                        pose,           // camera pose
                        500,500,        // fx, fy
                        500,500,        // cx, cy
                        NullLensDistortion());

  EXPECT_VECTOR_NEAR(pinhole.point_to_pixel(Vector3(0,0,10)),
                     Vector2(500,500), 1e-6);
  EXPECT_VECTOR_NEAR(pinhole.point_to_pixel(Vector3(-10,0,10)),
                     Vector2(0,500), 1e-6);
  EXPECT_VECTOR_NEAR(pinhole.point_to_pixel(Vector3(10,0,10)),
                     Vector2(1000,500), 1e-6);
  EXPECT_VECTOR_NEAR(pinhole.point_to_pixel(Vector3(0,-10,10)),
                     Vector2(500,1000), 1e-6);
  EXPECT_VECTOR_NEAR(pinhole.point_to_pixel(Vector3(0,10,10)),
                     Vector2(500,0), 1e-6);
  EXPECT_STREQ( "Pinhole", pinhole.type().c_str() );
}

TEST( PinholeModel, CoordinateFrame ) {
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

  EXPECT_VECTOR_NEAR(pinhole.point_to_pixel(Vector3(0,0,-10)),
                     Vector2(500,500), 1e-5);
  EXPECT_VECTOR_NEAR(pinhole.point_to_pixel(Vector3(-10,0,-10)),
                     Vector2(500,0), 1e-5);
  EXPECT_VECTOR_NEAR(pinhole.point_to_pixel(Vector3(10,0,-10)),
                     Vector2(500,1000), 1e-5);
  EXPECT_VECTOR_NEAR(pinhole.point_to_pixel(Vector3(0,-10,-10)),
                     Vector2(0,500), 1e-5);
  EXPECT_VECTOR_NEAR(pinhole.point_to_pixel(Vector3(0,10,-10)),
                     Vector2(1000,500), 1e-5);
}

TEST( PinholeModel, PixelToVector ) {
  Matrix<double,3,3> pose;
  pose.set_identity();

  // Create an imaginary 1000x1000 pixel imager
  PinholeModel pinhole( Vector3(0,0,0), // camera center
                        pose,           // camera pose
                        500,500,        // fx, fy
                        500,500,        // cx, cy
                        NullLensDistortion());

  PinholeModel pinhole2( Vector3(10,10,10), // camera center
                         pose,              // camera pose
                         500,500,           // fx, fy
                         500,500,           // cx, cy
                         NullLensDistortion());

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

  Vector2 zero_vector;
  EXPECT_VECTOR_NEAR( result1, zero_vector, 1e-8 );
  EXPECT_VECTOR_NEAR( result2, zero_vector, 1e-8 );
  EXPECT_VECTOR_NEAR( result3, zero_vector, 1e-8 );

#if defined(VW_HAVE_PKG_LAPACK) && VW_HAVE_PKG_LAPACK==1
  EXPECT_VECTOR_NEAR( result4, zero_vector, 1e-3 );
#endif
}

TEST( PinholeModel, TsaiDistortion ) {
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

  EXPECT_VECTOR_NEAR( distorted_pix,
                      Vector2(244.865,244.395),
                      1e-1 );
  EXPECT_VECTOR_NEAR( undistorted_pix,
                      Vector2(200,200),
                      1e-1 );
#endif
}

TEST( PinholeModel, ScalePinhole ) {
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

  Vector3 point = Vector3(2,-1,1) +
    5*pinhole4.pixel_to_vector(Vector2(500,500)) +
    pinhole4.camera_center();

  Vector2 o_return = pinhole4.point_to_pixel(point);
  Vector2 s_return = scaled.point_to_pixel(point);

  EXPECT_NEAR(o_return[0],s_return[0]*10,5); // Lens distortion doesn't
  EXPECT_NEAR(o_return[1],s_return[1]*10,5); // seem to be too accurate
}

TEST( PinholeModel, ProjectiveMatrix ) {
  // First set control camera
  Matrix<double,3,3> pose = math::euler_to_rotation_matrix(1.3,2.0,-.7,"xyz");

  // Create an imaginary 1000x1000 pixel imager
  PinholeModel control_pinhole( Vector3(0,4,-10),
                                pose, 600, 700,
                                500, 500,
                                NullLensDistortion() );

  // Make solve control camera (w/ random input for now)
  pose.set_identity();
  PinholeModel solved_pinhole( Vector3(-5,-2,5),
                               pose, 833, 544,
                               400, 700,
                               NullLensDistortion() );

  // Create Measurements used to solve for camera matrix
  std::vector<Vector<double> > world_m, image_m;
  boost::minstd_rand random_gen(42u);
  boost::normal_distribution<double> normal(0,20);
  boost::variate_generator<boost::minstd_rand&,
    boost::normal_distribution<double> > generator( random_gen, normal );
  for ( uint8 i = 0; i < 6; i++ ) {
    Vector3 point( generator(), generator(), generator()+60.0 );
    world_m.push_back( Vector4(point[0],point[1],point[2],1) );
    Vector2 pixel = control_pinhole.point_to_pixel(point);
    image_m.push_back( Vector3(pixel[0],pixel[1],1) );
  }

  // Building Camera Matrix
  CameraMatrixFittingFunctor fitfunc;
  Matrix<double> P = fitfunc(world_m,image_m);

  solved_pinhole.set_camera_matrix( P );

  // Compare camera matrices
  {
    Vector3 solved = solved_pinhole.camera_center( Vector2() );
    Vector3 control = control_pinhole.camera_center( Vector2() );
    EXPECT_VECTOR_NEAR( solved,
                        control, 1e-8 );
  }
  {
    Quaternion<double> solved = solved_pinhole.camera_pose( Vector2() );
    Quaternion<double> control = control_pinhole.camera_pose( Vector2() );
    for ( uint i = 0; i < 4; i ++ )
      EXPECT_NEAR( solved[i], control[i], 1e-8 );
  }
  {
    double s_fu, s_fv, s_cu, s_cv;
    double c_fu, c_fv, c_cu, c_cv;
    solved_pinhole.intrinsic_parameters( s_fu, s_fv, s_cu, s_cv );
    control_pinhole.intrinsic_parameters( c_fu, c_fv, c_cu, c_cv );
    EXPECT_NEAR( s_fv/s_fu, c_fv/c_fu, 1e-2 ); // Only accurate to 1%
    EXPECT_NEAR( s_cu/s_fu, c_cu/c_fu, 1e-2 );
    EXPECT_NEAR( s_cv/s_fu, c_cv/c_fu, 1e-2 );
  }
}

TEST( PinholeModel, BrownConradyDistortion ) {
  // First set control camera
  Matrix<double,3,3> pose = math::euler_to_rotation_matrix(1.3,2.0,-.7,"xyz");

  // Create an imaginary 1000x1000 pixel imager
  BrownConradyDistortion distortion( Vector2(-0.006,-0.002),
                                     Vector3(-.13361854e-5,
                                             0.52261757e-9,
                                             -.50728336e-13),
                                     Vector2(-.54958195e-6,
                                             -.46089420e-10),
                                     2.9659070 );
  PinholeModel control_pinhole( Vector3(0,4,-10), pose, 76.054, 76.054,
                                65, 65, distortion );
  for ( float i = 0; i < 100; i+=4 ) {
    for ( float j = 0; j < 100; j+=4 ) {
      Vector2 starting_pixel( i+7, j+7 );
      Vector3 point = control_pinhole.pixel_to_vector( starting_pixel );
      point = 50*point + control_pinhole.camera_center();
      Vector2 loop = control_pinhole.point_to_pixel( point );
      EXPECT_VECTOR_NEAR(loop, starting_pixel, 2e-6 );
    }
  }
}

#if defined(VW_HAVE_PKG_PROTOBUF) && VW_HAVE_PKG_PROTOBUF==1

// WriteRead TEST
class WriteReadTest : public ::testing::Test {
protected:
  WriteReadTest() {}

  virtual void SetUp() {
    // First set control camera
    Matrix<double,3,3> pose = math::euler_to_rotation_matrix(1.3,2.0,-.7,"xyz");
    control_pinhole = PinholeModel( Vector3(0,4,-10),
                                    pose, 600, 700,
                                    511, 533,
                                    NullLensDistortion() );
     control_R =
       control_pinhole.camera_pose().rotation_matrix();
  }

  void TestReadBack( std::string const& file ) {
    control_pinhole.write_file( file );
    PinholeModel read_back;
    read_back.read_file( file );
    EXPECT_VECTOR_DOUBLE_EQ( read_back.camera_center(),
                             control_pinhole.camera_center() );
    Matrix<double> read_back_R =
      read_back.camera_pose().rotation_matrix();
    EXPECT_MATRIX_DOUBLE_EQ( read_back_R, control_R );
    EXPECT_VECTOR_DOUBLE_EQ( read_back.coordinate_frame_u_direction(),
                             control_pinhole.coordinate_frame_u_direction() );
    EXPECT_VECTOR_DOUBLE_EQ( read_back.coordinate_frame_v_direction(),
                             control_pinhole.coordinate_frame_v_direction() );
    EXPECT_VECTOR_DOUBLE_EQ( read_back.coordinate_frame_w_direction(),
                             control_pinhole.coordinate_frame_w_direction() );
    EXPECT_STREQ( read_back.lens_distortion()->name().c_str(),
                  control_pinhole.lens_distortion()->name().c_str() );
    EXPECT_VECTOR_DOUBLE_EQ( read_back.lens_distortion()->distortion_parameters(),
                             control_pinhole.lens_distortion()->distortion_parameters() );
  }

  PinholeModel control_pinhole;
  Matrix<double> control_R;
};

TEST_F( WriteReadTest, Loop ) {
  { // NULL Distortion
    UnlinkName file("NullCam.tsai");
    TestReadBack( file );
  }

  control_pinhole.set_lens_distortion( TsaiLensDistortion( Vector4(-1,2,.3,.4) ) );
  { // TSAI Distortion
    UnlinkName file("TsaiCam.tsai");
    TestReadBack( file );
  }

  control_pinhole.set_lens_distortion( BrownConradyDistortion( Vector2(-1,3), Vector3(.05,1.052e-3,9.72564e-12), Vector2(-3.23415e-5,7.5271e-6), 2.345 ) );
  { // BROWN CONRADY Distortion
    UnlinkName file("BrownConrady.tsai");
    TestReadBack( file );
  }
}

#endif
