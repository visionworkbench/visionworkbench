// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
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

#include <vw/Camera/PinholeModel.h>
#include <vw/Math/Vector.h>

#include <vw/FileIO.h>
#include <vw/Camera/CameraTransform.h>
#include <vw/Math/EulerAngles.h>

using namespace std;
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

    //    std::cout << "\n" << pinhole << "\n";
    TS_ASSERT_EQUALS(pinhole.point_to_pixel(Vector3(0,0,10)),Vector2(500,500));
    TS_ASSERT_EQUALS(pinhole.point_to_pixel(Vector3(-10,0,10)),Vector2(0,500));
    TS_ASSERT_EQUALS(pinhole.point_to_pixel(Vector3(10,0,10)),Vector2(1000,500));
    TS_ASSERT_EQUALS(pinhole.point_to_pixel(Vector3(0,-10,10)),Vector2(500,1000));
    TS_ASSERT_EQUALS(pinhole.point_to_pixel(Vector3(0,10,10)),Vector2(500,0));

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
    
    //    std::cout << "\n" << pinhole << "\n";
//     std::cout << "\n\nPinhole 1\n";
//     std::cout << "Pixel to vector: " << pinhole.pixel_to_vector(Vector2(0,0)) << "\n";
//     std::cout << "Camera center: " << pinhole.camera_center(Vector2(0,0)) << "\n";
//     std::cout << "Point to pixel: " << pinhole.point_to_pixel(pinhole.pixel_to_vector(Vector2(0,0))+pinhole.camera_center(Vector2(0,0))) << "\n";

//     std::cout << "\n\nPinhole 2\n";
//     std::cout << "Pixel to vector: " << pinhole2.pixel_to_vector(Vector2(0,0)) << "\n";
//     std::cout << "Camera center: " << pinhole2.camera_center(Vector2(0,0)) << "\n";
//     std::cout << "Point to pixel: " << pinhole2.point_to_pixel(pinhole2.pixel_to_vector(Vector2(0,0))+pinhole2.camera_center(Vector2(0,0))) << "\n";

//     std::cout << "\n\nPinhole 3\n";
//     std::cout << "Pixel to vector: " << pinhole3.pixel_to_vector(Vector2(0,0)) << "\n";
//     std::cout << "Camera center: " << pinhole3.camera_center(Vector2(0,0)) << "\n";
//     std::cout << "Point to pixel: " << pinhole3.point_to_pixel(pinhole3.pixel_to_vector(Vector2(0,0))+pinhole3.camera_center(Vector2(0,0))) << "\n";

//     std::cout << "\n\nPinhole 4\n";
//     std::cout << "Pixel to vector: " << pinhole4.pixel_to_vector(Vector2(0,0)) << "\n";
//     std::cout << "Camera center: " << pinhole4.camera_center(Vector2(0,0)) << "\n";
//     std::cout << "Point to pixel: " << pinhole4.point_to_pixel(pinhole4.pixel_to_vector(Vector2(0,0))+pinhole4.camera_center(Vector2(0,0))) << "\n";

    Vector2 result1 = pinhole.point_to_pixel(pinhole.pixel_to_vector(Vector2(0,0))+pinhole.camera_center(Vector2(0,0)));
    Vector2 result2 = pinhole2.point_to_pixel(pinhole2.pixel_to_vector(Vector2(0,0))+pinhole2.camera_center(Vector2(0,0)));
    Vector2 result3 = pinhole3.point_to_pixel(pinhole3.pixel_to_vector(Vector2(0,0))+pinhole3.camera_center(Vector2(0,0)));
    Vector2 result4 = pinhole4.point_to_pixel(pinhole4.pixel_to_vector(Vector2(0,0))+pinhole4.camera_center(Vector2(0,0)));

    TS_ASSERT_EQUALS(result1, Vector2(0,0));
    TS_ASSERT_EQUALS(result2, Vector2(0,0));
    TS_ASSERT_DELTA(result3[0], Vector2(0,0)[0], 1e-8);
    TS_ASSERT_DELTA(result3[1], Vector2(0,0)[1], 1e-8);
    TS_ASSERT_DELTA(result4[0], Vector4(0,0)[0], 1e-3);
    TS_ASSERT_DELTA(result4[1], Vector4(0,0)[1], 1e-3);

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
    boost::shared_ptr<LensDistortion> distortion = pinhole.lens_distortion();
    
    Vector2 distorted_pix = distortion->get_distorted_coordinates(Vector2(200,200));
    Vector2 undistorted_pix = distortion->get_undistorted_coordinates(distorted_pix);
    
    TS_ASSERT_DELTA(distorted_pix[0], 244.968, 0.001);
    TS_ASSERT_DELTA(distorted_pix[1], 244.600, 0.001);
    TS_ASSERT_DELTA(undistorted_pix[0], 200, 0.001); 
    TS_ASSERT_DELTA(undistorted_pix[1], 200, 0.001); 
  }

//   void test_foo_distortion()
//   {
//         // Create an imaginary 1000x1000 pixel imager
//     PinholeModel distorted_pinhole_camera( Vector3(0,0,0),                 // camera center
//                                     math::identity_matrix<3>(),     // camera pose
//                                     500,500,                        // fx, fy
//                                     500,500,                        // cx, cy
//                                     TsaiLensDistortion(Vector4(-0.2805362343788147,
//                                                                0.1062035113573074,
//                                                                -0.0001422458299202845,
//                                                                0.00116333004552871))); 

//     ImageView<PixelRGB<uint8> > src_image;
//     read_image(src_image, "DSC_0264.jpg");

//     PinholeModel undistorted_pinhole_camera = linearize_camera(distorted_pinhole_camera);
//     ImageView<PixelRGB<uint8> > dst_image = transform(src_image,
//                                                       CameraTransform<PinholeModel,PinholeModel>(distorted_pinhole_camera,
//                                                                                                  undistorted_pinhole_camera));
//     write_image("test.jpg", dst_image);

//   }


};
