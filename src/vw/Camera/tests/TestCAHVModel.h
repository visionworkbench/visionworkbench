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
  }
};
