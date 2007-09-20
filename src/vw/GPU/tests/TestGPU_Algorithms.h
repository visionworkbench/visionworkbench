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

// Image/tests/TestAlgorithms.h
#include <cxxtest/TestSuite.h>

#include <vw/GPU.h>

using namespace std;
using namespace vw;
using namespace GPU;

#define DELTA_PRECISION 1e-5

class TestAlgorithms : public CxxTest::TestSuite
{
public:

  void test_image_algo_fill()
  {
	gpu_init();
    PixelRGB<float> fillval(1,2,3);
    GPUImage<PixelRGB<float> > fl(2,1);
    fill(fl, fillval);
    TS_ASSERT_DELTA( fl(0,0).r(), fillval.r(), DELTA_PRECISION);
    TS_ASSERT_DELTA( fl(0,0).g(), fillval.g(), DELTA_PRECISION);
    TS_ASSERT_DELTA( fl(0,0).b(), fillval.b(), DELTA_PRECISION);
    TS_ASSERT_DELTA( fl(1,0).r(), fillval.r(), DELTA_PRECISION);
    TS_ASSERT_DELTA( fl(1,0).g(), fillval.g(), DELTA_PRECISION);
    TS_ASSERT_DELTA( fl(1,0).b(), fillval.b(), DELTA_PRECISION);
  }

  void test_image_algo_clamp()
  {
    gpu_init();
    GPUImage<PixelRGB<float> > im(2,1);
    im.pixel(0,0) = PixelRGB<float>(0.5,0.9,1.5);
    im.pixel(1,0) = PixelRGB<float>(-1,0,1);

    GPUImage<PixelRGB<float> > c1 = clamp(im,-0.2,0.2);
    TS_ASSERT_DELTA( c1(0,0).r(), 0.2, DELTA_PRECISION );
    TS_ASSERT_DELTA( c1(0,0).g(), 0.2, DELTA_PRECISION );
    TS_ASSERT_DELTA( c1(0,0).b(), 0.2, DELTA_PRECISION );
    TS_ASSERT_DELTA( c1(1,0).r(), -0.2, DELTA_PRECISION );
    TS_ASSERT_DELTA( c1(1,0).g(), 0.0, DELTA_PRECISION );
    TS_ASSERT_DELTA( c1(1,0).b(), 0.2, DELTA_PRECISION );

    GPUImage<PixelRGB<float> > c2 = clamp(im,0.8);
    TS_ASSERT_DELTA( c2(0,0).r(), 0.5, DELTA_PRECISION );
    TS_ASSERT_DELTA( c2(0,0).g(), 0.8, DELTA_PRECISION );
    TS_ASSERT_DELTA( c2(0,0).b(), 0.8, DELTA_PRECISION );
    TS_ASSERT_DELTA( c2(1,0).r(), 0.0, DELTA_PRECISION );
    TS_ASSERT_DELTA( c2(1,0).g(), 0.0, DELTA_PRECISION );
    TS_ASSERT_DELTA( c2(1,0).b(), 0.8, DELTA_PRECISION );

    GPUImage<PixelRGB<float> > c3 = clamp(im);
    TS_ASSERT_DELTA( c3(0,0).r(), 0.5, DELTA_PRECISION );
    TS_ASSERT_DELTA( c3(0,0).g(), 0.9, DELTA_PRECISION );
    TS_ASSERT_DELTA( c3(0,0).b(), 1.0, DELTA_PRECISION );
    TS_ASSERT_DELTA( c3(1,0).r(), 0.0, DELTA_PRECISION );
    TS_ASSERT_DELTA( c3(1,0).g(), 0.0, DELTA_PRECISION );
    TS_ASSERT_DELTA( c3(1,0).b(), 1.0, DELTA_PRECISION );
  }

  void test_image_algo_normalize()
  {
    gpu_init();
    GPUImage<PixelRGB<float> > im(2,1);
    im.pixel(0,0) = PixelRGB<float>(0.5,1.5,2.5);
    im.pixel(1,0) = PixelRGB<float>(-1.5,0,1);

    GPUImage<PixelRGB<float> > n1 = normalize(im,-0.2,0.2);
    TS_ASSERT_DELTA( n1(0,0).r(), 0.0, DELTA_PRECISION );
    TS_ASSERT_DELTA( n1(0,0).g(), 0.1, DELTA_PRECISION );
    TS_ASSERT_DELTA( n1(0,0).b(), 0.2, DELTA_PRECISION );
    TS_ASSERT_DELTA( n1(1,0).r(), -0.2, DELTA_PRECISION );
    TS_ASSERT_DELTA( n1(1,0).g(), -0.05, DELTA_PRECISION );
    TS_ASSERT_DELTA( n1(1,0).b(), 0.05, DELTA_PRECISION );

    GPUImage<PixelRGB<float> > n2 = normalize(im, 0.8);
    TS_ASSERT_DELTA( n2(0,0).r(), 0.4, DELTA_PRECISION );
    TS_ASSERT_DELTA( n2(0,0).g(), 0.6, DELTA_PRECISION );
    TS_ASSERT_DELTA( n2(0,0).b(), 0.8, DELTA_PRECISION );
    TS_ASSERT_DELTA( n2(1,0).r(), 0.0, DELTA_PRECISION );
    TS_ASSERT_DELTA( n2(1,0).g(), 0.3, DELTA_PRECISION );
    TS_ASSERT_DELTA( n2(1,0).b(), 0.5, DELTA_PRECISION );

    GPUImage<PixelRGB<float> > n3 = normalize(im);
    TS_ASSERT_DELTA( n3(0,0).r(), 0.5, DELTA_PRECISION );
    TS_ASSERT_DELTA( n3(0,0).g(), 0.75, DELTA_PRECISION );
    TS_ASSERT_DELTA( n3(0,0).b(), 1.0, DELTA_PRECISION );
    TS_ASSERT_DELTA( n3(1,0).r(), 0.0, DELTA_PRECISION );
    TS_ASSERT_DELTA( n3(1,0).g(), 0.375, DELTA_PRECISION );
    TS_ASSERT_DELTA( n3(1,0).b(), 0.625, DELTA_PRECISION );
  }

  void test_image_algo_threshold()
  {
    gpu_init();
    GPUImage<PixelRGB<float> > im(2,1);
    im.pixel(0,0) = PixelRGB<float>(0.5,1.5,2.5);
    im.pixel(1,0) = PixelRGB<float>(-1.5,0,1);

    GPUImage<PixelRGB<float> > t1 = threshold(im,0.5,-0.2,0.2);
    TS_ASSERT_DELTA( t1(0,0).r(), -0.2, DELTA_PRECISION );
    TS_ASSERT_DELTA( t1(0,0).g(), 0.2, DELTA_PRECISION );
    TS_ASSERT_DELTA( t1(0,0).b(), 0.2, DELTA_PRECISION );
    TS_ASSERT_DELTA( t1(1,0).r(), -0.2, DELTA_PRECISION );
    TS_ASSERT_DELTA( t1(1,0).g(), -0.2, DELTA_PRECISION );
    TS_ASSERT_DELTA( t1(1,0).b(), 0.2, DELTA_PRECISION );

    GPUImage<PixelRGB<float> > t2 = threshold(im,0.6,0.8);
    TS_ASSERT_DELTA( t2(0,0).r(), 0.0, DELTA_PRECISION );
    TS_ASSERT_DELTA( t2(0,0).g(), 0.8, DELTA_PRECISION );
    TS_ASSERT_DELTA( t2(0,0).b(), 0.8, DELTA_PRECISION );
    TS_ASSERT_DELTA( t2(1,0).r(), 0.0, DELTA_PRECISION );
    TS_ASSERT_DELTA( t2(1,0).g(), 0.0, DELTA_PRECISION );
    TS_ASSERT_DELTA( t2(1,0).b(), 0.8, DELTA_PRECISION );

    GPUImage<PixelRGB<float> > t3 = threshold(im,0.6);
    TS_ASSERT_DELTA( t3(0,0).r(), 0.0, DELTA_PRECISION );
    TS_ASSERT_DELTA( t3(0,0).g(), 1.0, DELTA_PRECISION );
    TS_ASSERT_DELTA( t3(0,0).b(), 1.0, DELTA_PRECISION );
    TS_ASSERT_DELTA( t3(1,0).r(), 0.0, DELTA_PRECISION );
    TS_ASSERT_DELTA( t3(1,0).g(), 0.0, DELTA_PRECISION );
    TS_ASSERT_DELTA( t3(1,0).b(), 1.0, DELTA_PRECISION );

    GPUImage<PixelRGB<float> > t4 = threshold(im);
    TS_ASSERT_DELTA( t4(0,0).r(), 1.0, DELTA_PRECISION );
    TS_ASSERT_DELTA( t4(0,0).g(), 1.0, DELTA_PRECISION );
    TS_ASSERT_DELTA( t4(0,0).b(), 1.0, DELTA_PRECISION );
    TS_ASSERT_DELTA( t4(1,0).r(), 0.0, DELTA_PRECISION );
    TS_ASSERT_DELTA( t4(1,0).g(), 0.0, DELTA_PRECISION );
    TS_ASSERT_DELTA( t4(1,0).b(), 1.0, DELTA_PRECISION );
  }

}; // class TestAlgorithms
