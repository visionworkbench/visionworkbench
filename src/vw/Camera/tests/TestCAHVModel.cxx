// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>

#include <vw/Camera/CAHVModel.h>
#include <vw/Camera/PinholeModel.h>
#include <vw/Camera/CameraTransform.h>
#include <test/Helpers.h>

#include <cstdlib>
#include <ctime>

using namespace vw;
using namespace vw::camera;

TEST( CAHVModel, PinholeConversion ) {
  Matrix3x3 pose;
  pose(0,0) = 0.0665748562545205;
  pose(0,1) = -0.997208917590965;
  pose(0,2) = 0.0337959049552638;
  pose(1,0) = 0.992575027692545;
  pose(1,1) = 0.0627338358137274;
  pose(1,2) = -0.104207870361314;
  pose(2,0) = 0.101796870854825;
  pose(2,1) = 0.0404825952867597;
  pose(2,2) = 0.993981165094699;

  // Create fictitious Canon Rebel XT
  PinholeModel pinhole( Vector3( -57.315, 13.2, -106.947 ),
                        pose, 4210, 4210, 0, 0 );

  CAHVModel test_a(pinhole); // One routine
  CAHVModel test_b;
  test_b = pinhole;          // The other one

  for ( int x = 0; x < 10; x+=2 ){
    for ( int y = 0; y < 20; y+=4 ) {
      for ( int z = -5; z < 5; z+=3 ) {
        Vector3 test_point(x,y,z);
        Vector2 pin_px, cahv_a_px, cahv_b_px;
        pin_px = pinhole.point_to_pixel( test_point );
        cahv_a_px = test_a.point_to_pixel( test_point );
        cahv_b_px = test_b.point_to_pixel( test_point );

        EXPECT_VECTOR_NEAR( pin_px, cahv_a_px, 1e-4 );
        EXPECT_VECTOR_NEAR( pin_px, cahv_b_px, 1e-4 );
      }
    }
  }

  // Running another pinhole to CAHV conversion
  pose.set_identity();
  PinholeModel pinhold( Vector3(10,13,12),
                        pose, 2000, 2000, 3, -1 );

  CAHVModel test_c(pinhole);
  EXPECT_STREQ( "CAHV", test_c.type().c_str() );
  CAHVModel test_d;
  test_d = pinhole;

  for ( int x = 0; x < 10; x+=2 ){
    for ( int y = 0; y < 20; y+=4 ) {
      for ( int z = -5; z < 5; z+=3 ) {
        Vector3 test_point(x,y,z);
        Vector2 pin_px, cahv_c_px, cahv_d_px;
        pin_px = pinhole.point_to_pixel( test_point );
        cahv_c_px = test_c.point_to_pixel( test_point );
        cahv_d_px = test_d.point_to_pixel( test_point );

        EXPECT_VECTOR_NEAR( pin_px, cahv_c_px, 1e-4 );
        EXPECT_VECTOR_NEAR( pin_px, cahv_d_px, 1e-4 );
      }
    }
  }
}

// This is not really epipolar rectification but it does some of the
// work to force one image into the perspective of the other.
TEST( CAHVModel, FakeEpipolarConversion ) {

  // Building fake flat terrain
  std::vector<Vector3> points(100);
  for ( unsigned i = 0; i <100; i++ ){
    for ( unsigned s = 0; s < 3; s++ ) {
      if ( s == 2 ) {
        points[i][s] = 0;
      } else {
        points[i][s] = double(rand())/double(RAND_MAX) * 30 - 15;
      }
    }
  }

  // Building cameras
  CAHVModel model_a(2,Vector2(.1,.1),0,0,Vector3(30,0,30),
                    normalize(Vector3(-1,0,-1)),
                    Vector3(0,1,0),
                    cross_prod(normalize(Vector3(-1,0,-1)),Vector3(0,1,0)));
  CAHVModel model_b(2,Vector2(.1,.1),0,0,Vector3(-30,0,30),
                    normalize(Vector3(1,0,-1)),
                    Vector3(0,-1,0),
                    cross_prod(normalize(Vector3(1,0,-1)),Vector3(0,-1,0)));

  // Creating images
  std::vector<Vector2> image_a(100);
  std::vector<Vector2> image_b(100);
  for ( unsigned i = 0; i <100; i++ ) {
    image_a[i] = model_a.point_to_pixel(points[i]);
    image_b[i] = model_b.point_to_pixel(points[i]);
  }

  CAHVModel epimodel_a, epimodel_b;
  epipolar( model_a, model_b, epimodel_a, epimodel_b );

  // Building transform images
  std::vector<Vector2> new_image_a(100);
  std::vector<Vector2> new_image_b(100);
  for ( unsigned i = 0; i < 100; i++ ) {
    TransformRef trans_a( CameraTransform<CAHVModel, CAHVModel>(model_a, epimodel_a) );
    TransformRef trans_b( CameraTransform<CAHVModel, CAHVModel>(model_b, epimodel_b) );
    new_image_a[i] = trans_a.forward( image_a[i] );
    new_image_b[i] = trans_b.forward( image_b[i] );
  }

  // Calculate the average error now between images
  double sum = 0;
  for ( unsigned i = 0; i < 100; i++ ) {
    Vector2 temp;
    temp = new_image_a[i] - new_image_b[i];
    sum += norm_2(temp);

    EXPECT_NEAR( temp[0], -40, .1 ); // doesn't correct for translation
    EXPECT_NEAR( temp[1],   0, .1 ); // but it did correct rotation.
  }
}
