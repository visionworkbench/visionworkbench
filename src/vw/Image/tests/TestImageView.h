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

// TestImageView.h
#include <cxxtest/TestSuite.h>

#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>

using namespace std;
using namespace vw;

class TestImageView : public CxxTest::TestSuite
{
public:

  void testDefaultConstructor()
  {
    ImageView<double> test_double;
    TS_ASSERT_EQUALS( (bool)test_double, false );
    TS_ASSERT_EQUALS( !test_double, true );
    TS_ASSERT_EQUALS(test_double.cols(), 0);
    TS_ASSERT_EQUALS(test_double.rows(), 0);
    TS_ASSERT_EQUALS(test_double.planes(), 0);
    TS_ASSERT_EQUALS(test_double.data(), (double*)0);

    ImageView<PixelRGBA<uint8> > test_rgba;
    TS_ASSERT_EQUALS( (bool)test_rgba, false );
    TS_ASSERT_EQUALS( !test_rgba, true );
    TS_ASSERT_EQUALS(test_rgba.cols(), 0);
    TS_ASSERT_EQUALS(test_rgba.rows(), 0);
    TS_ASSERT_EQUALS(test_rgba.planes(), 0);
    TS_ASSERT_EQUALS(test_rgba.data(), (PixelRGBA<uint8>*)0);
  }

  void testColsRowsConstructor()
  {
    ImageView<double> test_double(3,4);
    TS_ASSERT_EQUALS( (bool)test_double, true );
    TS_ASSERT_EQUALS( !test_double, false );
    TS_ASSERT_EQUALS(test_double.cols(), 3);
    TS_ASSERT_EQUALS(test_double.rows(), 4);
    TS_ASSERT_EQUALS(test_double.planes(), 1);
    TS_ASSERT_DIFFERS(test_double.data(), (double*)0);

    ImageView<PixelRGBA<uint8> > test_rgba(3,4);
    TS_ASSERT_EQUALS( (bool)test_rgba, true );
    TS_ASSERT_EQUALS( !test_rgba, false );
    TS_ASSERT_EQUALS(test_rgba.cols(), 3);
    TS_ASSERT_EQUALS(test_rgba.rows(), 4);
    TS_ASSERT_EQUALS(test_rgba.planes(), 1);
    TS_ASSERT_DIFFERS(test_rgba.data(), (PixelRGBA<uint8>*)0);
  }

  void testColsRowsPlanesConstructor()
  {
    ImageView<double> test_double(4,3,2);
    TS_ASSERT_EQUALS( (bool)test_double, true );
    TS_ASSERT_EQUALS( !test_double, false );
    TS_ASSERT_EQUALS(test_double.cols(), 4);
    TS_ASSERT_EQUALS(test_double.rows(), 3);
    TS_ASSERT_EQUALS(test_double.planes(), 2);
    TS_ASSERT_DIFFERS(test_double.data(), (double*)0);

    ImageView<PixelRGBA<uint8> > test_rgba(4,3,2);
    TS_ASSERT_EQUALS( (bool)test_rgba, true );
    TS_ASSERT_EQUALS( !test_rgba, false );
    TS_ASSERT_EQUALS(test_rgba.cols(), 4);
    TS_ASSERT_EQUALS(test_rgba.rows(), 3);
    TS_ASSERT_EQUALS(test_rgba.planes(), 2);
    TS_ASSERT_DIFFERS(test_rgba.data(), (PixelRGBA<uint8>*)0);
  }

  void testCopyConstructor()
  {
    ImageView<double> test_double(3,4);
    ImageView<double> test2_double( test_double );
    TS_ASSERT_EQUALS( (bool)test2_double, true );
    TS_ASSERT_EQUALS( !test2_double, false );
    TS_ASSERT_EQUALS(test2_double.cols(), 3);
    TS_ASSERT_EQUALS(test2_double.rows(), 4);
    TS_ASSERT_EQUALS(test2_double.planes(), 1);
    TS_ASSERT_EQUALS(test2_double.data(),test_double.data());

    ImageView<PixelRGBA<uint8> > test_rgba(3,4);
    ImageView<PixelRGBA<uint8> > test2_rgba( test_rgba );
    TS_ASSERT_EQUALS( (bool)test2_rgba, true );
    TS_ASSERT_EQUALS( !test2_rgba, false );
    TS_ASSERT_EQUALS(test2_rgba.cols(), 3);
    TS_ASSERT_EQUALS(test2_rgba.rows(), 4);
    TS_ASSERT_EQUALS(test2_rgba.planes(), 1);
    TS_ASSERT_EQUALS(test2_rgba.data(), test_rgba.data());
  }

  void testSetSize()
  {
    ImageView<double> test_double(3,4);
    TS_ASSERT_EQUALS( (bool)test_double, true );
    test_double.set_size(2,3);
    TS_ASSERT_EQUALS( (bool)test_double, true );
    TS_ASSERT_EQUALS(test_double.cols(), 2);
    TS_ASSERT_EQUALS(test_double.rows(), 3);
    TS_ASSERT_EQUALS(test_double.planes(), 1);
    test_double.set_size(0,0);
    TS_ASSERT_EQUALS( (bool)test_double, false );
    TS_ASSERT_EQUALS(test_double.cols(), 0);
    TS_ASSERT_EQUALS(test_double.rows(), 0);
    TS_ASSERT_EQUALS(test_double.data(), (double*)0);

    ImageView<PixelRGBA<uint8> > test_rgba(3,4);
    TS_ASSERT_EQUALS( (bool)test_rgba, true );
    test_rgba.set_size(2,3);
    TS_ASSERT_EQUALS( (bool)test_rgba, true );
    TS_ASSERT_EQUALS(test_rgba.cols(), 2);
    TS_ASSERT_EQUALS(test_rgba.rows(), 3);
    TS_ASSERT_EQUALS(test_rgba.planes(), 1);
    test_rgba.set_size(0,0);
    TS_ASSERT_EQUALS( (bool)test_rgba, false );
    TS_ASSERT_EQUALS(test_rgba.cols(), 0);
    TS_ASSERT_EQUALS(test_rgba.rows(), 0);
    TS_ASSERT_EQUALS(test_rgba.data(), (PixelRGBA<uint8>*)0);
  }

  void testReset()
  {
    ImageView<double> test_double(3,4);
    TS_ASSERT_EQUALS( (bool)test_double, true );
    test_double.reset();
    TS_ASSERT_EQUALS( (bool)test_double, false );
    TS_ASSERT_EQUALS(test_double.cols(), 0);
    TS_ASSERT_EQUALS(test_double.rows(), 0);
    TS_ASSERT_EQUALS(test_double.planes(), 0);
    TS_ASSERT_EQUALS(test_double.data(), (double*)0);

    ImageView<PixelRGBA<uint8> > test_rgba(3,4);
    TS_ASSERT_EQUALS( (bool)test_rgba, true );
    test_rgba.reset();
    TS_ASSERT_EQUALS( (bool)test_rgba, false );
    TS_ASSERT_EQUALS(test_rgba.cols(), 0);
    TS_ASSERT_EQUALS(test_rgba.rows(), 0);
    TS_ASSERT_EQUALS(test_rgba.planes(), 0);
    TS_ASSERT_EQUALS(test_rgba.data(), (PixelRGBA<uint8>*)0);
  }

  void testRasterization()
  {
    ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
    ImageView<double> im2(2,2);
    TS_ASSERT_THROWS_NOTHING( im1.rasterize( im2, BBox2i(0,0,2,2) ) );
    TS_ASSERT_EQUALS( im1(0,0), im2(0,0) );
    TS_ASSERT_EQUALS( im1(1,0), im2(1,0) );
    TS_ASSERT_EQUALS( im1(0,1), im2(0,1) );
    TS_ASSERT_EQUALS( im1(1,1), im2(1,1) );
    TS_ASSERT_DIFFERS( &im1(0,0), &im2(0,0) );
    TS_ASSERT_DIFFERS( &im1(1,0), &im2(1,0) );
    TS_ASSERT_DIFFERS( &im1(0,1), &im2(0,1) );
    TS_ASSERT_DIFFERS( &im1(1,1), &im2(1,1) );
  }

  void testDefaultRasterization1()
  {
    ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
    ImageView<double> im2(2,2);
    TS_ASSERT_THROWS_NOTHING( vw::rasterize( im1, im2, BBox2i(0,0,2,2) ) );
    TS_ASSERT_EQUALS( im1(0,0), im2(0,0) );
    TS_ASSERT_EQUALS( im1(1,0), im2(1,0) );
    TS_ASSERT_EQUALS( im1(0,1), im2(0,1) );
    TS_ASSERT_EQUALS( im1(1,1), im2(1,1) );
    TS_ASSERT_DIFFERS( &im1(0,0), &im2(0,0) );
    TS_ASSERT_DIFFERS( &im1(1,0), &im2(1,0) );
    TS_ASSERT_DIFFERS( &im1(0,1), &im2(0,1) );
    TS_ASSERT_DIFFERS( &im1(1,1), &im2(1,1) );
  }

  void testDefaultRasterization2()
  {
    ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
    ImageView<double> im2(2,2);
    TS_ASSERT_THROWS_NOTHING( vw::rasterize( im1, im2 ) );
    TS_ASSERT_EQUALS( im1(0,0), im2(0,0) );
    TS_ASSERT_EQUALS( im1(1,0), im2(1,0) );
    TS_ASSERT_EQUALS( im1(0,1), im2(0,1) );
    TS_ASSERT_EQUALS( im1(1,1), im2(1,1) );
    TS_ASSERT_DIFFERS( &im1(0,0), &im2(0,0) );
    TS_ASSERT_DIFFERS( &im1(1,0), &im2(1,0) );
    TS_ASSERT_DIFFERS( &im1(0,1), &im2(0,1) );
    TS_ASSERT_DIFFERS( &im1(1,1), &im2(1,1) );
  }

  void testDefaultRasterization3()
  {
    ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
    ImageView<double> im2;
    TS_ASSERT_THROWS_NOTHING( vw::rasterize( im1, im2, BBox2i(0,0,2,2) ) );
    TS_ASSERT_EQUALS( im1(0,0), im2(0,0) );
    TS_ASSERT_EQUALS( im1(1,0), im2(1,0) );
    TS_ASSERT_EQUALS( im1(0,1), im2(0,1) );
    TS_ASSERT_EQUALS( im1(1,1), im2(1,1) );
    TS_ASSERT_DIFFERS( &im1(0,0), &im2(0,0) );
    TS_ASSERT_DIFFERS( &im1(1,0), &im2(1,0) );
    TS_ASSERT_DIFFERS( &im1(0,1), &im2(0,1) );
    TS_ASSERT_DIFFERS( &im1(1,1), &im2(1,1) );
  }

  void testDefaultRasterization4()
  {
    ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
    ImageView<double> im2;
    TS_ASSERT_THROWS_NOTHING( vw::rasterize( im1, im2 ) );
    TS_ASSERT_EQUALS( im1(0,0), im2(0,0) );
    TS_ASSERT_EQUALS( im1(1,0), im2(1,0) );
    TS_ASSERT_EQUALS( im1(0,1), im2(0,1) );
    TS_ASSERT_EQUALS( im1(1,1), im2(1,1) );
    TS_ASSERT_DIFFERS( &im1(0,0), &im2(0,0) );
    TS_ASSERT_DIFFERS( &im1(1,0), &im2(1,0) );
    TS_ASSERT_DIFFERS( &im1(0,1), &im2(0,1) );
    TS_ASSERT_DIFFERS( &im1(1,1), &im2(1,1) );
  }

  void testPartialRasterization()
  {
    ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
    ImageView<double> im2(1,2);
    TS_ASSERT_THROWS_NOTHING( im1.rasterize( im2, BBox2i(1,0,1,2) ) );
    TS_ASSERT_EQUALS( im1(1,0), im2(0,0) );
    TS_ASSERT_EQUALS( im1(1,1), im2(0,1) );
    TS_ASSERT_DIFFERS( &im1(1,0), &im2(0,0) );
    TS_ASSERT_DIFFERS( &im1(1,1), &im2(0,1) );
    ImageView<double> im3(2,1);
    TS_ASSERT_THROWS_NOTHING( im1.rasterize( im3, BBox2i(0,1,2,1) ) );
    TS_ASSERT_EQUALS( im1(0,1), im3(0,0) );
    TS_ASSERT_EQUALS( im1(1,1), im3(1,0) );
    TS_ASSERT_DIFFERS( &im1(0,1), &im3(0,0) );
    TS_ASSERT_DIFFERS( &im1(1,1), &im3(1,0) );
  }

  void testIterator()
  {
    ImageView<double> im(2,2,2);
    im(0,0,0)=1; im(1,0,0)=2; im(0,1,0)=3; im(1,1,0)=4;
    im(0,0,1)=5; im(1,0,1)=6; im(0,1,1)=7; im(1,1,1)=8;
    ImageView<double>::iterator i = im.begin();
    TS_ASSERT_DIFFERS( i, im.end() );
    TS_ASSERT_EQUALS( &(*i), &(im(0,0,0)) );
    TS_ASSERT_EQUALS( i.col(), 0 );
    TS_ASSERT_EQUALS( i.row(), 0 );
    TS_ASSERT_EQUALS( i.plane(), 0 );
    ++i;
    TS_ASSERT_DIFFERS( i, im.end() );
    TS_ASSERT_EQUALS( &(*i), &(im(1,0,0)) );
    TS_ASSERT_EQUALS( i.col(), 1 );
    TS_ASSERT_EQUALS( i.row(), 0 );
    TS_ASSERT_EQUALS( i.plane(), 0 );
    ++i;
    TS_ASSERT_DIFFERS( i, im.end() );
    TS_ASSERT_EQUALS( &(*i), &(im(0,1,0)) );
    TS_ASSERT_EQUALS( i.col(), 0 );
    TS_ASSERT_EQUALS( i.row(), 1 );
    TS_ASSERT_EQUALS( i.plane(), 0 );
    ++i;
    TS_ASSERT_DIFFERS( i, im.end() );
    TS_ASSERT_EQUALS( &(*i), &(im(1,1,0)) );
    TS_ASSERT_EQUALS( i.col(), 1 );
    TS_ASSERT_EQUALS( i.row(), 1 );
    TS_ASSERT_EQUALS( i.plane(), 0 );
    ++i;
    TS_ASSERT_DIFFERS( i, im.end() );
    TS_ASSERT_EQUALS( &(*i), &(im(0,0,1)) );
    TS_ASSERT_EQUALS( i.col(), 0 );
    TS_ASSERT_EQUALS( i.row(), 0 );
    TS_ASSERT_EQUALS( i.plane(), 1 );
    ++i;
    TS_ASSERT_DIFFERS( i, im.end() );
    TS_ASSERT_EQUALS( &(*i), &(im(1,0,1)) );
    TS_ASSERT_EQUALS( i.col(), 1 );
    TS_ASSERT_EQUALS( i.row(), 0 );
    TS_ASSERT_EQUALS( i.plane(), 1 );
    ++i;
    TS_ASSERT_DIFFERS( i, im.end() );
    TS_ASSERT_EQUALS( &(*i), &(im(0,1,1)) );
    TS_ASSERT_EQUALS( i.col(), 0 );
    TS_ASSERT_EQUALS( i.row(), 1 );
    TS_ASSERT_EQUALS( i.plane(), 1 );
    ++i;
    TS_ASSERT_DIFFERS( i, im.end() );
    TS_ASSERT_EQUALS( &(*i), &(im(1,1,1)) );
    TS_ASSERT_EQUALS( i.col(), 1 );
    TS_ASSERT_EQUALS( i.row(), 1 );
    TS_ASSERT_EQUALS( i.plane(), 1 );
    ++i;
    TS_ASSERT_EQUALS( i, im.end() );
  }

#ifdef VW_IMAGE_BOUNDS_CHECK
  void testBoundsCheck()
  {
    ImageView<PixelRGB<uint8> > test_img(2,2);

    // First, test to make sure that the operator() bounds checking is
    // working properly.
    PixelRGB<uint8> test = test_img(0,0);
    TS_ASSERT_THROWS_NOTHING ( test = test_img(1,0) );
    TS_ASSERT_THROWS_NOTHING ( test = test_img(0,1) );
    TS_ASSERT_THROWS_NOTHING ( test = test_img(1,1) );
    TS_ASSERT_THROWS_ANYTHING ( test = test_img(0,2) );
    TS_ASSERT_THROWS_ANYTHING ( test = test_img(2,0) );
    TS_ASSERT_THROWS_ANYTHING ( test = test_img(-1,0) );
    TS_ASSERT_THROWS_ANYTHING ( test = test_img(0,-1) );
    
    // Next, test to make sure that the MemoryStridingPixelAccessor
    // bounds checking is working properly.
    ImageView<PixelRGB<uint8> >::pixel_accessor acc = test_img.origin();

    TS_ASSERT_THROWS_NOTHING( test = *acc );

    // First, check to make sure that negative indices cause an error.
    acc.prev_col();
    TS_ASSERT_THROWS_ANYTHING ( test = *acc );
    
    // Proceed through colums
    acc = test_img.origin();
    acc.next_col();
    TS_ASSERT_THROWS_NOTHING ( test = *acc );
    acc.next_col();
    TS_ASSERT_THROWS_NOTHING ( test = *acc );
    acc.next_col();
    TS_ASSERT_THROWS_NOTHING ( test = *acc );
    acc.next_col();
    TS_ASSERT_THROWS_ANYTHING ( test = *acc );

    // Proceed through rows
    acc = test_img.origin();
    acc.next_row();
    TS_ASSERT_THROWS_NOTHING ( test = *acc );
    acc.next_row();
    TS_ASSERT_THROWS_ANYTHING ( test = *acc );
  }
#endif

};
