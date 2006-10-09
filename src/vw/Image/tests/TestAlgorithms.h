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

#include <vw/Image/Algorithms.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>

using namespace std;
using namespace vw;

class TestAlgorithms : public CxxTest::TestSuite
{
public:

  void test_image_algo_fill()
  {
    PixelRGB<double> fillval(1,2,3);
    ImageView<PixelRGB<double> > fl(2,1);
    fill(fl, fillval);
    TS_ASSERT_EQUALS( fl(0,0).r(), fillval.r());
    TS_ASSERT_EQUALS( fl(0,0).g(), fillval.g());
    TS_ASSERT_EQUALS( fl(0,0).b(), fillval.b());
    TS_ASSERT_EQUALS( fl(1,0).r(), fillval.r());
    TS_ASSERT_EQUALS( fl(1,0).g(), fillval.g());
    TS_ASSERT_EQUALS( fl(1,0).b(), fillval.b());
  }

  void test_image_algo_clamp()
  {
    ImageView<PixelRGB<double> > im(2,1);
    im(0,0) = PixelRGB<double>(0.5,0.9,1.5);
    im(1,0) = PixelRGB<double>(-1,0,1);

    ImageView<PixelRGB<double> > c1 = clamp(im,-0.2,0.2);
    TS_ASSERT_EQUALS( c1(0,0).r(), 0.2 );
    TS_ASSERT_EQUALS( c1(0,0).g(), 0.2 );
    TS_ASSERT_EQUALS( c1(0,0).b(), 0.2 );
    TS_ASSERT_EQUALS( c1(1,0).r(), -0.2 );
    TS_ASSERT_EQUALS( c1(1,0).g(), 0.0 );
    TS_ASSERT_EQUALS( c1(1,0).b(), 0.2 );

    ImageView<PixelRGB<double> > c2 = clamp(im,0.8);
    TS_ASSERT_EQUALS( c2(0,0).r(), 0.5 );
    TS_ASSERT_EQUALS( c2(0,0).g(), 0.8 );
    TS_ASSERT_EQUALS( c2(0,0).b(), 0.8 );
    TS_ASSERT_EQUALS( c2(1,0).r(), 0.0 );
    TS_ASSERT_EQUALS( c2(1,0).g(), 0.0 );
    TS_ASSERT_EQUALS( c2(1,0).b(), 0.8 );

    ImageView<PixelRGB<double> > c3 = clamp(im);
    TS_ASSERT_EQUALS( c3(0,0).r(), 0.5 );
    TS_ASSERT_EQUALS( c3(0,0).g(), 0.9 );
    TS_ASSERT_EQUALS( c3(0,0).b(), 1.0 );
    TS_ASSERT_EQUALS( c3(1,0).r(), 0.0 );
    TS_ASSERT_EQUALS( c3(1,0).g(), 0.0 );
    TS_ASSERT_EQUALS( c3(1,0).b(), 1.0 );
  }

  void test_image_algo_normalize()
  {
    ImageView<PixelRGB<double> > im(2,1);
    im(0,0) = PixelRGB<double>(0.5,1.5,2.5);
    im(1,0) = PixelRGB<double>(-1.5,0,1);

    ImageView<PixelRGB<double> > n1 = normalize(im,-0.2,0.2);
    TS_ASSERT_DELTA( n1(0,0).r(), 0.0, 1e-7 );
    TS_ASSERT_DELTA( n1(0,0).g(), 0.1, 1e-7 );
    TS_ASSERT_DELTA( n1(0,0).b(), 0.2, 1e-7 );
    TS_ASSERT_DELTA( n1(1,0).r(), -0.2, 1e-7 );
    TS_ASSERT_DELTA( n1(1,0).g(), -0.05, 1e-7 );
    TS_ASSERT_DELTA( n1(1,0).b(), 0.05, 1e-7 );

    ImageView<PixelRGB<double> > n2 = normalize(im,0.8);
    TS_ASSERT_DELTA( n2(0,0).r(), 0.4, 1e-7 );
    TS_ASSERT_DELTA( n2(0,0).g(), 0.6, 1e-7 );
    TS_ASSERT_DELTA( n2(0,0).b(), 0.8, 1e-7 );
    TS_ASSERT_DELTA( n2(1,0).r(), 0.0, 1e-7 );
    TS_ASSERT_DELTA( n2(1,0).g(), 0.3, 1e-7 );
    TS_ASSERT_DELTA( n2(1,0).b(), 0.5, 1e-7 );

    ImageView<PixelRGB<double> > n3 = normalize(im);
    TS_ASSERT_DELTA( n3(0,0).r(), 0.5, 1e-7 );
    TS_ASSERT_DELTA( n3(0,0).g(), 0.75, 1e-7 );
    TS_ASSERT_DELTA( n3(0,0).b(), 1.0, 1e-7 );
    TS_ASSERT_DELTA( n3(1,0).r(), 0.0, 1e-7 );
    TS_ASSERT_DELTA( n3(1,0).g(), 0.375, 1e-7 );
    TS_ASSERT_DELTA( n3(1,0).b(), 0.625, 1e-7 );
  }

  void test_image_algo_threshold()
  {
    ImageView<PixelRGB<double> > im(2,1);
    im(0,0) = PixelRGB<double>(0.5,1.5,2.5);
    im(1,0) = PixelRGB<double>(-1.5,0,1);

    ImageView<PixelRGB<double> > t1 = threshold(im,0.5,-0.2,0.2);
    TS_ASSERT_EQUALS( t1(0,0).r(), -0.2 );
    TS_ASSERT_EQUALS( t1(0,0).g(), 0.2 );
    TS_ASSERT_EQUALS( t1(0,0).b(), 0.2 );
    TS_ASSERT_EQUALS( t1(1,0).r(), -0.2 );
    TS_ASSERT_EQUALS( t1(1,0).g(), -0.2 );
    TS_ASSERT_EQUALS( t1(1,0).b(), 0.2 );

    ImageView<PixelRGB<double> > t2 = threshold(im,0.6,0.8);
    TS_ASSERT_EQUALS( t2(0,0).r(), 0.0 );
    TS_ASSERT_EQUALS( t2(0,0).g(), 0.8 );
    TS_ASSERT_EQUALS( t2(0,0).b(), 0.8 );
    TS_ASSERT_EQUALS( t2(1,0).r(), 0.0 );
    TS_ASSERT_EQUALS( t2(1,0).g(), 0.0 );
    TS_ASSERT_EQUALS( t2(1,0).b(), 0.8 );

    ImageView<PixelRGB<double> > t3 = threshold(im,0.6);
    TS_ASSERT_EQUALS( t3(0,0).r(), 0.0 );
    TS_ASSERT_EQUALS( t3(0,0).g(), 1.0 );
    TS_ASSERT_EQUALS( t3(0,0).b(), 1.0 );
    TS_ASSERT_EQUALS( t3(1,0).r(), 0.0 );
    TS_ASSERT_EQUALS( t3(1,0).g(), 0.0 );
    TS_ASSERT_EQUALS( t3(1,0).b(), 1.0 );

    ImageView<PixelRGB<double> > t4 = threshold(im);
    TS_ASSERT_EQUALS( t4(0,0).r(), 1.0 );
    TS_ASSERT_EQUALS( t4(0,0).g(), 1.0 );
    TS_ASSERT_EQUALS( t4(0,0).b(), 1.0 );
    TS_ASSERT_EQUALS( t4(1,0).r(), 0.0 );
    TS_ASSERT_EQUALS( t4(1,0).g(), 0.0 );
    TS_ASSERT_EQUALS( t4(1,0).b(), 1.0 );
  }

}; // class TestAlgorithms
