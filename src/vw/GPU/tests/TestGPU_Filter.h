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

// TestFilter.h
#include <cxxtest/TestSuite.h>

#include <vw/GPU.h>

#include <vector>

using namespace std;
using namespace vw;
using namespace GPU;

class TestFilter : public CxxTest::TestSuite
{
public:

  void test_gaussian_filter() {
	gpu_init(true, true);
    GPUImage<float> src(2,2); src.pixel(0,0)=1; src.pixel(1,0)=2; src.pixel(0,1)=3; src.pixel(1,1)=4;
    GPUImage<float> dst = gaussian_filter( src, 1.0, 0.0, 5.0, 0.0);
	TS_ASSERT_DELTA( (float) dst(0,0), 0.3877403988*1+0.2447702197*2, 1e-3 );
    TS_ASSERT_DELTA( (float) dst(1,0), 0.3877403988*2+0.2447702197*1, 1e-3 );
    TS_ASSERT_DELTA( (float) dst(0,1), 0.3877403988*3+0.2447702197*4, 1e-3 );
    TS_ASSERT_DELTA( (float) dst(1,1), 0.3877403988*4+0.2447702197*3, 1e-3 );
  }

   void test_derivative_filter() {
    std::vector<float> x_kernel, y_kernel;
    generate_derivative_kernel(x_kernel, 1, 0);
    generate_derivative_kernel(y_kernel, 0, 0);
	printf("x_kernel: ");
	for(int i=0; i < x_kernel.size(); i++)
		printf("%f ", x_kernel[i]);
	printf("\n");
	printf("y_kernel: ");
	for(int i=0; i < y_kernel.size(); i++)
		printf("%f ", y_kernel[i]);
	printf("\n");
   
	gpu_init(true, true);
    GPUImage<float> src(2,2); 
	src.pixel(0,0)=1; src.pixel(1,0)=2; src.pixel(0,1)=3; src.pixel(1,1)=4;
    GPUImage<float> dst = derivative_filter(src, 1, 0 ); // Implicit ConstantEdgeExtend, which the results are based on
    TS_ASSERT_EQUALS( (float) dst(0,0), 0.5 );
    TS_ASSERT_EQUALS( (float) dst(1,0), 0.5 );
    TS_ASSERT_EQUALS( (float) dst(0,1), 0.5 );
    TS_ASSERT_EQUALS( (float) dst(1,1), 0.5 );
  }

  void test_laplacian_filter() {
	// Problem: I think long shaders don't run properly on the mac. It gives annoying straightforward looking garbage values.
	gpu_init(true, true);
    GPUImage<float> src(2,2); src.pixel(0,0)=1; src.pixel(1,0)=2; src.pixel(0,1)=3; src.pixel(1,1)=4;
    GPUImage<float> dst = laplacian_filter(src);
    TS_ASSERT_EQUALS( (float) dst(0,0), 1 );
    TS_ASSERT_EQUALS( (float) dst(1,0), -3 );
    TS_ASSERT_EQUALS( (float) dst(0,1), -7 );
    TS_ASSERT_EQUALS( (float) dst(1,1), -11 );
  }

}; // class TestFilter
