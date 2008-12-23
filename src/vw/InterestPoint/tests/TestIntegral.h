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

#include <cxxtest/TestSuite.h>
#include <vw/InterestPoint/IntegralImage.h>
#include <vw/Image.h>
#include <vw/FileIO.h>

using namespace vw;
using namespace vw::ip;

class TestIntegral : public CxxTest::TestSuite
{
 public:

  void test_integral_summing() {
    ImageView<float> graffiti;
    read_image( graffiti, "sub.png" );
    ImageView<double> integral;
    integral = IntegralImage( graffiti );
    
    for (unsigned size = 10; size <= 70; size+= 10) {
      float sum = IntegralBlock( integral,
				 Vector2i(0,0),
				 Vector2i(size,size) );
      float actual_sum = 0;
      for (unsigned i = 0; i < size; i++ ) {
	for (unsigned j = 0; j < size; j++ ) {
	  actual_sum+=graffiti(i,j);
	}
      }

      TS_ASSERT_DELTA(actual_sum, sum, .01);
    }
  }

  void test_haar_filters() { 
    ImageView<float> graffiti;
    read_image( graffiti, "sub.png" );
    ImageView<double> integral;
    integral = IntegralImage( graffiti );
    
    float first_H = HHaarWavelet( integral,
				  49.5, 49.5, 20 );
    float first_V = VHaarWavelet( integral,
				  49.5, 49.5, 20 );
    
    ImageView<float> rotated;
    rotated = rotate_180(graffiti);
    integral = IntegralImage( rotated );
    
    TS_ASSERT_DELTA( HHaarWavelet( integral,
				   49.5, 49.5, 20 ),
		     -first_H, .1 );
    TS_ASSERT_DELTA( VHaarWavelet( integral,
				   49.5, 49.5, 20 ),
		     -first_V, .1 );

    float hand_response = 0;
    for ( unsigned i = 0; i < 2; i++ ) {
      for ( unsigned j = 0; j < 2; j++ ) {
	if ( i == 0 ) {
	  hand_response -= rotated(i,j);
	} else {
	  hand_response += rotated(i,j);
	}
      }
    }
    TS_ASSERT_DELTA( hand_response,
		     HHaarWavelet( integral,
				   1,1,2 ),
		     .01 );
    
    hand_response = 0;
    for ( unsigned i = 0; i < 2; i++ ) {
      for ( unsigned j = 0; j < 2; j++ ) {
	if ( j == 0 ) {
	  hand_response -= rotated(i,j);
	} else {
	  hand_response += rotated(i,j);
	}
      }
    }
    TS_ASSERT_DELTA( hand_response,
		     VHaarWavelet( integral,
				   1, 1, 2 ),
		     .01 );
  
  }
};

