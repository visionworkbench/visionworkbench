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


// TestPinholeModel.h
#include <gtest/gtest_VW.h>
#include <vw/Math/Vector.h>
#include <vw/Math/EulerAngles.h>
#include <vw/Math/LinearAlgebra.h>
#include <vw/Camera/CameraGeometry.h>
#include <vw/Camera/CameraTransform.h>
#include <vw/Camera/LensDistortion.h>
#include <vw/Camera/PinholeModel.h>
#include <test/Helpers.h>

#include <boost/random.hpp>

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
                        500,500);       // cx, cy

  EXPECT_VECTOR_NEAR(pinhole.point_to_pixel(Vector3(0,0,10)),
                     Vector2(500,500), 1e-6);
  EXPECT_VECTOR_NEAR(pinhole.point_to_pixel(Vector3(-10,0,10)),
                     Vector2(0,500), 1e-6);
  EXPECT_VECTOR_NEAR(pinhole.point_to_pixel(Vector3(10,0,10)),
                     Vector2(1000,500), 1e-6);
  EXPECT_VECTOR_NEAR(pinhole.point_to_pixel(Vector3(0,-10,10)),
                     Vector2(500,0), 1e-6);
  EXPECT_VECTOR_NEAR(pinhole.point_to_pixel(Vector3(0,10,10)),
                     Vector2(500,1000), 1e-6);
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
                        500,500,        // cx, cy
                        Vector3(0, 1, 0),
                        Vector3(1, 0, 0),
                        Vector3(0, 0, -1));

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

TEST( PinholeModel, TsaiDistortion ) {
  double distortion_arr[] = {-0.2805362343788147, 0.1062035113573074,
                             -0.0001422458299202845, 0.00116333004552871};
  Vector<double> distortion_vec(sizeof(distortion_arr)/sizeof(double), distortion_arr);
  TsaiLensDistortion lens(distortion_vec);
  // Create an imaginary 1000x1000 pixel imager
  PinholeModel pinhole( Vector3(0,0,0),                 // camera center
                        math::identity_matrix<3>(),     // camera pose
                        500,500,                        // fx, fy
                        500,500,                        // cx, cy
                        &lens);
  const LensDistortion* distortion = pinhole.lens_distortion();

#if defined(VW_HAVE_PKG_LAPACK)
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
  double distortion_arr[] = {-0.2796604335308075, 0.1031486615538597,
                             -0.0007824968779459596, 0.0009675505571067333};
  Vector<double> distortion_vec(sizeof(distortion_arr)/sizeof(double), distortion_arr);
  TsaiLensDistortion lens(distortion_vec);
  PinholeModel pinhole4(Vector3(-0.329, 0.065, -0.82),
                        rot,
                        605.320556640625,
                        606.3638305664062,
                        518.89208984375,
                        387.5555114746094,
                        Vector3(1, 0, 0),
                        Vector3(0, -1, 0),
                        Vector3(0, 0, 1),
                        &lens);
  PinholeModel scaled = scale_camera(pinhole4, .1);

  Vector3 point = Vector3(2,-1,1) +
    5*pinhole4.pixel_to_vector(Vector2(500,500)) +
    pinhole4.camera_center();

  Vector2 o_return = pinhole4.point_to_pixel(point);
  Vector2 s_return = scaled.point_to_pixel(point);

  EXPECT_NEAR(o_return[0],s_return[0]*10,1);
  EXPECT_NEAR(o_return[1],s_return[1]*10,1);
}

TEST( PinholeModel, ProjectiveMatrix ) {
  // First set control camera
  Matrix<double,3,3> pose = math::euler_to_rotation_matrix(1.3,2.0,-.7,"xyz");

  // Create an imaginary 1000x1000 pixel imager
  PinholeModel control_pinhole( Vector3(0,4,-10),
                                pose, 600, 700,
                                500, 500);

  // Make solve control camera (w/ random input for now)
  pose.set_identity();
  PinholeModel solved_pinhole( Vector3(-5,-2,5),
                               pose, 833, 544,
                               400, 700);

  // Create Measurements used to solve for camera matrix
  std::vector<Vector<double> > world_m, image_m;
  boost::minstd_rand random_gen(42u);
  boost::normal_distribution<double> normal(0,20);
  boost::variate_generator<boost::minstd_rand&,
    boost::normal_distribution<double> > generator( random_gen, normal );
  while(world_m.size() < 6) {
    try {
      Vector3 point( generator(), generator(), generator()+60.0 );
      Vector2 pixel = control_pinhole.point_to_pixel(point);

      world_m.push_back( Vector4(point[0],point[1],point[2],1) );
      image_m.push_back( Vector3(pixel[0],pixel[1],1) );
    } catch(vw::camera::PointToPixelErr) {}
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
    for ( uint32 i = 0; i < 4; i ++ )
      EXPECT_NEAR( solved[i], control[i], 1e-8 );
  }
  {
    Vector2 s_focal = solved_pinhole.focal_length();
    Vector2 c_focal = control_pinhole.focal_length();
    EXPECT_NEAR( s_focal[0]/s_focal[1],c_focal[0]/c_focal[1], 1e-2 );
    EXPECT_VECTOR_NEAR( elem_quot(solved_pinhole.point_offset(),s_focal),
                        elem_quot(control_pinhole.point_offset(),c_focal), 1e-2 );
  }
}

/// Helper class for testing PinholeModel objects
class PinholeTest : public ::testing::Test {
protected:
  PinholeTest() {}

  /// Init with a set of constants and no lens distortion
  virtual void SetUp() {
    expect_rot = vw::math::euler_to_rotation_matrix(1.15, 0.0, -1.57, "xyz");
    pinhole = PinholeModel(Vector3(-0.329, 0.065, -0.82),
                           expect_rot,
                           605.320556640625,
                           606.3638305664062,
                           518.89208984375,
                           387.5555114746094);
  }

  /// Iterate through a grid of pixels and make sure that the pixels project
  ///  out and then back in to the same pixel.
  void projection_test(double tolerance=1e-6) {
    Vector2 image_size = pinhole.point_offset();
    image_size *= 2;
    for ( unsigned x = 10; x < image_size.x(); x+=80 ) {
      for ( unsigned y = 10; y < image_size.y(); y+=80 ) {
          EXPECT_VECTOR_NEAR( Vector2(x,y),
                              pinhole.point_to_pixel(pinhole.pixel_to_vector(Vector2(x,y)) +
                                                     pinhole.camera_center(Vector2(x,y))),
                            tolerance );
      }
    }
  }

  /// Write a .tsai file, then read it back in and make sure nothing has changed.
  void readback_test(std::string const& file) {
    pinhole.write( file );
    PinholeModel read_back;
    read_back.read( file );
    EXPECT_VECTOR_DOUBLE_EQ( read_back.camera_center(),
                             pinhole.camera_center() );
    Matrix<double> read_back_R =
      read_back.camera_pose().rotation_matrix();
    EXPECT_MATRIX_NEAR( read_back_R, expect_rot,1e-8 );
    EXPECT_VECTOR_DOUBLE_EQ( read_back.coordinate_frame_u_direction(),
                             pinhole.coordinate_frame_u_direction() );
    EXPECT_VECTOR_DOUBLE_EQ( read_back.coordinate_frame_v_direction(),
                             pinhole.coordinate_frame_v_direction() );
    EXPECT_VECTOR_DOUBLE_EQ( read_back.coordinate_frame_w_direction(),
                             pinhole.coordinate_frame_w_direction() );
    EXPECT_STREQ( read_back.lens_distortion()->name().c_str(),
                  pinhole.lens_distortion()->name().c_str() );
    EXPECT_VECTOR_DOUBLE_EQ( read_back.lens_distortion()->distortion_parameters(),
                             pinhole.lens_distortion()->distortion_parameters() );
  }

  PinholeModel pinhole;
  Matrix<double> expect_rot;
};

TEST_F( PinholeTest, NullLensDistortion ) {
  projection_test();
  UnlinkName file("NullCam.tsai");
  readback_test( file );

  // check enforcement that pose returns the rotation from camera
  // frame to world frame.
  Vector2 center_pixel = pinhole.point_offset();
  Quaternion<double> center_pose =
    pinhole.camera_pose(center_pixel);
  double angle_from_z =
    acos(dot_prod(Vector3(0,0,1),inverse(center_pose).rotate(pinhole.pixel_to_vector(center_pixel))));
  EXPECT_LT( angle_from_z, 0.5 );
}

TEST_F( PinholeTest, TsaiLensDistortion ) {
  double distortion_arr[] = {-0.2796604335308075, 0.1031486615538597,
                             -0.0007824968779459596, 0.0009675505571067333};
  Vector<double> distortion_vec(sizeof(distortion_arr)/sizeof(double), distortion_arr);
  TsaiLensDistortion lens(distortion_vec);
  pinhole.set_lens_distortion(&lens);
#if defined(VW_HAVE_PKG_LAPACK) && VW_HAVE_PKG_LAPACK==1
  projection_test(1e-4);
#endif
  UnlinkName file("TsaiCam.tsai");
  readback_test( file );
}

TEST_F( PinholeTest, BrownConradyDistortion ) {
  BrownConradyDistortion lens(Vector2(-0.6,-0.2),
                              Vector3(.1336185e-8,
                                      -0.5226175e-12,
                                      0),
                              Vector2(.5495819e-9,
                                      0),
                              0.201);
  pinhole.set_lens_distortion(&lens);
#if defined(VW_HAVE_PKG_LAPACK)
  projection_test(1e-4);
#endif
  UnlinkName file("BrownConrady.tsai");
  readback_test( file );
}

TEST_F( PinholeTest, AdjustableTsaiDistortion ) {
  Vector<double> distort_coeff(6);
  distort_coeff[0] = 0.007646500298509824;   // k1
  distort_coeff[1] = -0.01743067138801845;   // k2
  distort_coeff[2] = 0.00980946292640812;    // k3
  distort_coeff[3] = -2.98092556225311e-05;  // p1
  distort_coeff[4] = -1.339089765674149e-05; // p2
  distort_coeff[5] = -1.221974557659228e-05; // alpha = skew
  AdjustableTsaiLensDistortion lens(distort_coeff);
  pinhole.set_lens_distortion(&lens);
#if defined(VW_HAVE_PKG_LAPACK)
  projection_test(1e-4);
#endif
  UnlinkName file("AdjustedTsai.tsai");
  readback_test( file );
}

TEST_F( PinholeTest, OldFormatReadTest ) {
  UnlinkName filename("monkey.tsai");
  std::ofstream filestream( filename.c_str() );
  filestream << "fu = 54.6\n";
  filestream << "fv = 45.3\n";
  filestream << "cu = 3\n";
  filestream << "cv = 5\n";
  filestream << "u_direction = 1 0 0\n";
  filestream << "v_direction = 0 1 0\n";
  filestream << "w_direction = 0 0 1\n";
  filestream << "C = 18.6 14.4 13.3\n";
  filestream << "R = 1 0 0 0 1 0 0 0 1\n";
  filestream << "k1 = 0.001\n";
  filestream << "k2 = 0.001\n";
  filestream << "p1 = 0.010\n";
  filestream << "p2 = 1\n";
  filestream.close();

  PinholeModel monkey;
  ASSERT_NO_THROW( monkey.read( filename ) );
  EXPECT_VECTOR_NEAR( Vector2(54.6,45.3),
                      monkey.focal_length(), 1e-5 );
  EXPECT_VECTOR_NEAR( Vector2(3,5),
                      monkey.point_offset(), 1e-5 );
  EXPECT_VECTOR_NEAR( Vector3(18.6,14.4,13.3),
                      monkey.camera_center(Vector2()), 1e-5 );
}
