// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


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

// TestCAHVModel.h
#include <cxxtest/TestSuite.h>

#include <vw/Camera/PinholeModel.h>
#include <vw/Camera/CAHVModel.h>
#include <vw/Math.h>
#include <vw/Image/Transform.h>
#include <vw/Camera/CameraTransform.h>
#include <vw/FileIO.h>

#include <stdlib.h>
#include <time.h>

using namespace vw;
using namespace vw::camera;

class TestCAHVModel : public CxxTest::TestSuite
{
 public:
  void test_pinhole_conversion ()
  {
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

	  TS_ASSERT_DELTA( pin_px[0], cahv_a_px[0], 0.0001 );
	  TS_ASSERT_DELTA( pin_px[1], cahv_a_px[1], 0.0001 );
	  TS_ASSERT_DELTA( pin_px[0], cahv_b_px[0], 0.0001 );
	  TS_ASSERT_DELTA( pin_px[1], cahv_b_px[1], 0.0001 );
	}
      }
    }

    // Running another pinhole to CAHV conversion
    pose.set_identity();
    PinholeModel pinhold( Vector3(10,13,12),
			  pose, 2000, 2000, 3, -1 );

    CAHVModel test_c(pinhole);
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

	  TS_ASSERT_DELTA( pin_px[0], cahv_c_px[0], 0.0001 );
	  TS_ASSERT_DELTA( pin_px[1], cahv_c_px[1], 0.0001 );
	  TS_ASSERT_DELTA( pin_px[0], cahv_d_px[0], 0.0001 );
	  TS_ASSERT_DELTA( pin_px[1], cahv_d_px[1], 0.0001 );
	}
      }
    }
  }

  void test_CAHV_pinhole() {
    // Building fake flat terrain
    srandom((unsigned int) clock());
    std::vector<Vector3> points(100);
    for ( unsigned i = 0; i <100; i++ ){
      for ( unsigned s = 0; s < 3; s++ ) {
	if ( s == 2 ) {
	  points[i][s] = 0;
	} else {
	  points[i][s] = double(random())/double(pow(2,31)-1) * 30 - 15;
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
    
    // Writing debug images
    {
      /*
      BBox2 bb_image_a;
      BBox2 bb_image_b;
      for ( unsigned i = 0; i < 100; i++ ) {
	bb_image_a.grow(image_a[i]);
	bb_image_b.grow(image_b[i]);
      }
      
      Vector2 size_a = bb_image_a.size() + Vector2(10,10);
      Vector2 size_b = bb_image_b.size() + Vector2(10,10);
      ImageView<PixelRGB<uint8> > px_image_a(4*size_a.x(),4*size_a.y()), px_image_b(4*size_b.x(),4*size_b.y());
      
      for ( unsigned i = 0; i < 100; i++ ) {
	Vector2i location_a = 4*image_a[i] + 2*size_a;
	Vector2i location_b = 4*image_b[i] + 2*size_b;
	std::cout << "p_a:" << location_a << " p_b:" << location_b << std::endl;
	
	if ( i < 50 ) {
	  px_image_a(location_a.x(), location_a.y()) = PixelRGB<uint8>(255-4*i,
								       255,4*i);
	  px_image_b(location_b.x(), location_b.y()) = PixelRGB<uint8>(255-4*i,
								       255,4*i);
	} else {
	  px_image_a(location_a.x(), location_a.y()) = PixelRGB<uint8>(255-4*(i-50),
								       0,4*(i-50));
	  px_image_b(location_b.x(), location_b.y()) = PixelRGB<uint8>(255-4*(i-50),
								       0,4*(i-50));
	}
      }
      write_image("image_a.png", px_image_a);
      write_image("image_b.png", px_image_b);
      */
    }

    // Creating epipolar cameras
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

    // Writing debug images of the epipolar stuff
    {
      /*
      BBox2 bb_image_a;
      BBox2 bb_image_b;
      for ( unsigned i = 0; i < 100; i++ ) {
	bb_image_a.grow(new_image_a[i]);
	bb_image_b.grow(new_image_b[i]);
      }
      
      Vector2 size_a = bb_image_a.size() + Vector2(10,10);
      Vector2 size_b = bb_image_b.size() + Vector2(10,10);
      ImageView<PixelRGB<uint8> > px_image_a(4*size_a.x(),4*size_a.y()), px_image_b(4*size_b.x(),4*size_b.y());
      
      for ( unsigned i = 0; i < 100; i++ ) {
	Vector2i location_a = 4*new_image_a[i] + 2*size_a;
	Vector2i location_b = 4*new_image_b[i] + 2*size_b;
	std::cout << "p_a:" << location_a << " p_b:" << location_b << std::endl;
	
	if ( i < 50 ) {
	  px_image_a(location_a.x(), location_a.y()) = PixelRGB<uint8>(255-4*i,
								       255,4*i);
	  px_image_b(location_b.x(), location_b.y()) = PixelRGB<uint8>(255-4*i,
								       255,4*i);
	} else {
	  px_image_a(location_a.x(), location_a.y()) = PixelRGB<uint8>(255-4*(i-50),
								       0,4*(i-50));
	  px_image_b(location_b.x(), location_b.y()) = PixelRGB<uint8>(255-4*(i-50),
								       0,4*(i-50));
	}
      }
      write_image("epiimage_a.png", px_image_a);
      write_image("epiimage_b.png", px_image_b);
      */
    }

    // Calculate the average error now between images
    double sum = 0;
    for ( unsigned i = 0; i < 100; i++ ) {
      Vector2 temp;
      temp = new_image_a[i] - new_image_b[i];
      sum += norm_2(temp);

      TS_ASSERT_DELTA( temp[0], -40, .1 ); // Epipolar doesn't correct for translation
      TS_ASSERT_DELTA( temp[1], 0, .1 );   // but it did correct rotation.
    }
    sum /= 100;
    TS_TRACE("Average error out: " + stringify(sum) + "\n");
  }
};
