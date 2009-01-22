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

// TestTransform.h
#include <cxxtest/TestSuite.h>

#include <vw/Image.h>
#include <vw/GPU.h>

#define DELTA_PRECISION 1e-4

using namespace vw;
using namespace GPU;

class TestTransform : public CxxTest::TestSuite
{
 public:

  void testTranslateTransform() // NOT DONE
    {
      /*
	GPUImage<float> im(2,3); im.pixel(0,0)=1; im.pixel(1,0)=2; im.pixel(0,1)=3; im.pixel(1,1)=4; im.pixel(0,2)=5; im.pixel(1,2)=6;
	GPUImage<float> im2 = translate(im, TranslateTransform(1,1));
	TransformView<InterpolationView<EdgeExtendView<GPUImage<float>, ZeroEdgeExtend>, BilinearInterpolation>, TranslateTransform> im2 = transform(im, TranslateTransform(1,1));
	TS_ASSERT_DELTA( im2.cols(), 2 );
	TS_ASSERT_DELTA( im2.rows(), 3 );
	TS_ASSERT_DELTA( im2(1,1), 1 );
	TS_ASSERT_DELTA( im2(0,0), 0 );
	TS_ASSERT_DELTA( im2(1,2), 3 );
	TS_ASSERT_DELTA( im2(-1,-1), 0 );
	TS_ASSERT_DELTA( im2(-100,-100), 0 );

	// Test accessor
	TS_ASSERT_DELTA( *(im2.origin().advance(1,1)), 1 );
	TS_ASSERT_DELTA( *(im2.origin().advance(-1,-1)), 0 );

	// Test the traits
	TS_ASSERT( bool_trait<IsFloatingPointIndexable>( im2 ) );
	TS_ASSERT( bool_trait<IsFloatingPointIndexable>( transform(im, TranslateTransform(1,1)) ) );
	TS_ASSERT( !bool_trait<IsMultiplyAccessible>( transform(im, TranslateTransform(1,1)) ) );
      */
    }

  void testHomographyTransform() {

  }

#define print_image(im) {	\
	for(int y=0; y < im.rows(); y++) {	\
		for(int x=0; x < im.cols(); x++) {	\
			printf("%.10f  ", im(x, y));	\
		}	\
		printf("\n");	\
	}	\
}

  void testResample_BilinearInterpolation()
    {
      gpu_init();
      GPUImage<float> im(2,3); im.pixel(0,0)=1; im.pixel(1,0)=2; im.pixel(0,1)=3; im.pixel(1,1)=4; im.pixel(0,2)=5; im.pixel(1,2)=6;
      GPUImage<float> im2 = resample(im, 2, 2, 4, 6, BilinearInterpolation(), ZeroEdgeExtend());
      TS_ASSERT_DELTA( im2.cols(), 4, DELTA_PRECISION);
      TS_ASSERT_DELTA( im2.rows(), 6, DELTA_PRECISION);
      TS_ASSERT_DELTA( im2(0,0), 1, DELTA_PRECISION);
      TS_ASSERT_DELTA( im2(0,2), 3, DELTA_PRECISION);
      TS_ASSERT_DELTA( im2(2,2), 4, DELTA_PRECISION);
      TS_ASSERT_DELTA( im2(1,1), 2.5, DELTA_PRECISION);
      TS_ASSERT_DELTA( im2(1,2), 3.5, DELTA_PRECISION);
    }
};
