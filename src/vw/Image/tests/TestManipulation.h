// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
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

// TestManipulation.h
#include <cxxtest/TestSuite.h>

#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/Manipulation.h>

using namespace vw;

class TestImageView : public CxxTest::TestSuite
{
public:

  void testCopy()
  {
    ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
    ImageView<double> im2 = copy(im);
    TS_ASSERT_EQUALS( im2.cols(), 2 );
    TS_ASSERT_EQUALS( im2.rows(), 3 );
    TS_ASSERT_EQUALS( im2(0,0), 1 );
    TS_ASSERT_EQUALS( im2(1,0), 2 );
    TS_ASSERT_EQUALS( im2(0,1), 3 );
    TS_ASSERT_EQUALS( im2(1,1), 4 );
    TS_ASSERT_EQUALS( im2(0,2), 5 );
    TS_ASSERT_EQUALS( im2(1,2), 6 );
    // Make sure it's really deep.
    TS_ASSERT_DIFFERS( im2.data(), im.data() );
    TS_ASSERT_EQUALS( copy(im)(1,0), im(1,0) );
    TS_ASSERT_DIFFERS( &(copy(im)(1,0)), &(im(1,0)) );
  }

  void testTranspose()
  {
    ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
    ImageView<double> im2 = transpose(im);
    TS_ASSERT_EQUALS( im2.cols(), 3 );
    TS_ASSERT_EQUALS( im2.rows(), 2 );
    TS_ASSERT_EQUALS( im2(0,0), 1 );
    TS_ASSERT_EQUALS( im2(1,0), 3 );
    TS_ASSERT_EQUALS( im2(2,0), 5 );
    TS_ASSERT_EQUALS( im2(0,1), 2 );
    TS_ASSERT_EQUALS( im2(1,1), 4 );
    TS_ASSERT_EQUALS( im2(2,1), 6 );
    // Make sure it's really shallow.
    TS_ASSERT_EQUALS( transpose(im)(1,0), im(0,1) );
    TS_ASSERT_EQUALS( &(transpose(im)(1,0)), &(im(0,1)) );
  }

  void testRotate180()
  {
    ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
    ImageView<double> im2 = rotate_180(im);
    TS_ASSERT_EQUALS( im2.cols(), 2 );
    TS_ASSERT_EQUALS( im2.rows(), 3 );
    TS_ASSERT_EQUALS( im2(0,0), 6 );
    TS_ASSERT_EQUALS( im2(1,0), 5 );
    TS_ASSERT_EQUALS( im2(0,1), 4 );
    TS_ASSERT_EQUALS( im2(1,1), 3 );
    TS_ASSERT_EQUALS( im2(0,2), 2 );
    TS_ASSERT_EQUALS( im2(1,2), 1 );
    // Make sure it's really shallow.
    TS_ASSERT_EQUALS( rotate_180(im)(1,0), im(0,2) );
    TS_ASSERT_EQUALS( &(rotate_180(im)(1,0)), &(im(0,2)) );
  }

  void testRotate90CW()
  {
    ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
    ImageView<double> im2 = rotate_90_cw(im);
    TS_ASSERT_EQUALS( im2.cols(), 3 );
    TS_ASSERT_EQUALS( im2.rows(), 2 );
    TS_ASSERT_EQUALS( im2(0,0), 5 );
    TS_ASSERT_EQUALS( im2(1,0), 3 );
    TS_ASSERT_EQUALS( im2(2,0), 1 );
    TS_ASSERT_EQUALS( im2(0,1), 6 );
    TS_ASSERT_EQUALS( im2(1,1), 4 );
    TS_ASSERT_EQUALS( im2(2,1), 2 );
    // Make sure it's really shallow.
    TS_ASSERT_EQUALS( rotate_90_cw(im)(2,1), im(1,0) );
    TS_ASSERT_EQUALS( &(rotate_90_cw(im)(2,1)), &(im(1,0)) );
  }

  void testRotate90CCW()
  {
    ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
    ImageView<double> im2 = rotate_90_ccw(im);
    TS_ASSERT_EQUALS( im2.cols(), 3 );
    TS_ASSERT_EQUALS( im2.rows(), 2 );
    TS_ASSERT_EQUALS( im2(0,0), 2 );
    TS_ASSERT_EQUALS( im2(1,0), 4 );
    TS_ASSERT_EQUALS( im2(2,0), 6 );
    TS_ASSERT_EQUALS( im2(0,1), 1 );
    TS_ASSERT_EQUALS( im2(1,1), 3 );
    TS_ASSERT_EQUALS( im2(2,1), 5 );
    // Make sure it's really shallow.
    TS_ASSERT_EQUALS( rotate_90_ccw(im)(0,0), im(1,0) );
    TS_ASSERT_EQUALS( &(rotate_90_ccw(im)(0,0)), &(im(1,0)) );
  }

  void testFlipVertical()
  {
    ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
    ImageView<double> im2 = flip_vertical(im);
    TS_ASSERT_EQUALS( im2.cols(), 2 );
    TS_ASSERT_EQUALS( im2.rows(), 3 );
    TS_ASSERT_EQUALS( im2(0,0), 5 );
    TS_ASSERT_EQUALS( im2(1,0), 6 );
    TS_ASSERT_EQUALS( im2(0,1), 3 );
    TS_ASSERT_EQUALS( im2(1,1), 4 );
    TS_ASSERT_EQUALS( im2(0,2), 1 );
    TS_ASSERT_EQUALS( im2(1,2), 2 );
    // Make sure it's really shallow.
    TS_ASSERT_EQUALS( flip_vertical(im)(1,0), im(1,2) );
    TS_ASSERT_EQUALS( &(flip_vertical(im)(1,0)), &(im(1,2)) );
  }
  
  void testFlipHorizontal()
  {
    ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
    ImageView<double> im2 = flip_horizontal(im);
    TS_ASSERT_EQUALS( im2.cols(), 2 );
    TS_ASSERT_EQUALS( im2.rows(), 3 );
    TS_ASSERT_EQUALS( im2(0,0), 2 );
    TS_ASSERT_EQUALS( im2(1,0), 1 );
    TS_ASSERT_EQUALS( im2(0,1), 4 );
    TS_ASSERT_EQUALS( im2(1,1), 3 );
    TS_ASSERT_EQUALS( im2(0,2), 6 );
    TS_ASSERT_EQUALS( im2(1,2), 5 );
    // Make sure it's really shallow.
    TS_ASSERT_EQUALS( flip_horizontal(im)(1,0), im(0,0) );
    TS_ASSERT_EQUALS( &(flip_horizontal(im)(1,0)), &(im(0,0)) );
  }
  
  void testCrop()
  {
    ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
    ImageView<double> im2 = crop(im,1,1,1,2);
    TS_ASSERT_EQUALS( im2.cols(), 1 );
    TS_ASSERT_EQUALS( im2.rows(), 2 );
    TS_ASSERT_EQUALS( im2(0,0), 4 );
    TS_ASSERT_EQUALS( im2(0,1), 6 );

    ImageView<double> im3 = crop(im,BBox2i(1,1,1,2));
    TS_ASSERT_EQUALS( im3.cols(), 1 );
    TS_ASSERT_EQUALS( im3.rows(), 2 );
    TS_ASSERT_EQUALS( im3(0,0), 4 );
    TS_ASSERT_EQUALS( im3(0,1), 6 );

    // Make sure it's really shallow.
    TS_ASSERT_EQUALS( crop(im,1,1,1,2)(0,0), im(1,1) );
    TS_ASSERT_EQUALS( &(crop(im,1,1,1,2)(0,0)), &(im(1,1)) );

  }

  void testSubsample()
  {
    ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
    ImageView<double> im2 = subsample(im,1,2);
    TS_ASSERT_EQUALS( im2.cols(), 2 );
    TS_ASSERT_EQUALS( im2.rows(), 2 );
    TS_ASSERT_EQUALS( im2(0,0), 1 );
    TS_ASSERT_EQUALS( im2(1,0), 2 );
    TS_ASSERT_EQUALS( im2(0,1), 5 );
    TS_ASSERT_EQUALS( im2(1,1), 6 );
    // Make sure it's really shallow.
    TS_ASSERT_EQUALS( subsample(im,2)(0,1), im(0,2) );
    TS_ASSERT_EQUALS( &(subsample(im,2)(0,1)), &(im(0,2)) );
  }

  void testSelectCol()
  {
    ImageView<double> im(2,2); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4;
    ImageView<double> im2 = select_col(im,1);
    TS_ASSERT_EQUALS( im2.cols(), 1 );
    TS_ASSERT_EQUALS( im2.rows(), 2 );
    TS_ASSERT_EQUALS( im2.planes(), 1 );
    TS_ASSERT_EQUALS( im2(0,0), 2 );
    TS_ASSERT_EQUALS( im2(0,1), 4 );
    // Make sure it's really shallow.
    TS_ASSERT_EQUALS( select_col(im,1)(0,1), im(1,1) );
    TS_ASSERT_EQUALS( &(select_col(im,1)(0,1)), &(im(1,1)) );
  }

  void testSelectRow()
  {
    ImageView<double> im(2,2); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4;
    ImageView<double> im2 = select_row(im,1);
    TS_ASSERT_EQUALS( im2.cols(), 2 );
    TS_ASSERT_EQUALS( im2.rows(), 1 );
    TS_ASSERT_EQUALS( im2.planes(), 1 );
    TS_ASSERT_EQUALS( im2(0,0), 3 );
    TS_ASSERT_EQUALS( im2(1,0), 4 );
    // Make sure it's really shallow.
    TS_ASSERT_EQUALS( select_row(im,1)(1,0), im(1,1) );
    TS_ASSERT_EQUALS( &(select_row(im,1)(1,0)), &(im(1,1)) );
  }

  void testSelectPlane()
  {
    ImageView<double> im(1,2,2); im(0,0,0)=1; im(0,1,0)=2; im(0,0,1)=3; im(0,1,1)=4;
    ImageView<double> im2 = select_plane(im,1);
    TS_ASSERT_EQUALS( im2.cols(), 1 );
    TS_ASSERT_EQUALS( im2.rows(), 2 );
    TS_ASSERT_EQUALS( im2.planes(), 1 );
    TS_ASSERT_EQUALS( im2(0,0), 3 );
    TS_ASSERT_EQUALS( im2(0,1), 4 );
    // Make sure it's really shallow.
    TS_ASSERT_EQUALS( select_plane(im,0)(0,1), im(0,1) );
    TS_ASSERT_EQUALS( &(select_plane(im,0)(0,1)), &(im(0,1)) );
  }

  void testChannelsToPlanes()
  {
    ImageView<PixelRGB<double> > im(1,2); im(0,0)=PixelRGB<double>(1,2,3); im(0,1)=PixelRGB<double>(4,5,6);
    ImageView<double> im2 = channels_to_planes(im);
    TS_ASSERT_EQUALS( im2.cols(), 1 );
    TS_ASSERT_EQUALS( im2.rows(), 2 );
    TS_ASSERT_EQUALS( im2.planes(), 3 );
    TS_ASSERT_EQUALS( im2(0,0,0), 1 );
    TS_ASSERT_EQUALS( im2(0,1,0), 4 );
    TS_ASSERT_EQUALS( im2(0,0,1), 2 );
    TS_ASSERT_EQUALS( im2(0,1,1), 5 );
    TS_ASSERT_EQUALS( im2(0,0,2), 3 );
    TS_ASSERT_EQUALS( im2(0,1,2), 6 );
    // Make sure it's really shallow.
    TS_ASSERT_EQUALS( channels_to_planes(im)(0,1,1), im(0,1)[1] );
    TS_ASSERT_EQUALS( &(channels_to_planes(im)(0,1,1)), &(im(0,1)[1]) );
  }

  void testPlanesToChannels()
  {
    ImageView<double> im(1,2,3); im(0,0,0)=1; im(0,0,1)=2; im(0,0,2)=3; im(0,1,0)=4; im(0,1,1)=5; im(0,1,2)=6;
    ImageView<PixelRGB<double> > im2 = planes_to_channels<PixelRGB<double> >(im);
    TS_ASSERT_EQUALS( im2.cols(), 1 );
    TS_ASSERT_EQUALS( im2.rows(), 2 );
    TS_ASSERT_EQUALS( im2.planes(), 1 );
    TS_ASSERT_EQUALS( im2(0,0)[0], 1 );
    TS_ASSERT_EQUALS( im2(0,1)[0], 4 );
    TS_ASSERT_EQUALS( im2(0,0)[1], 2 );
    TS_ASSERT_EQUALS( im2(0,1)[1], 5 );
    TS_ASSERT_EQUALS( im2(0,0)[2], 3 );
    TS_ASSERT_EQUALS( im2(0,1)[2], 6 );
  }

};
