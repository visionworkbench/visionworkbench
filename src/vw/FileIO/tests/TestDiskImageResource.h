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

// TestDiskImageResource.h
#define CXXTEST_ABORT_TEST_ON_FAIL
#include <cxxtest/TestSuite.h>
#include <vw/FileIO/DiskImageResource.h>

#include <vw/config.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/ImageView.h>

using namespace vw;

class TestDiskImageResource : public CxxTest::TestSuite
{
public:

  void test_read_image_rgb_png_uint8() {
#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1
    ImageView<PixelRGB<uint8> > image;
    TS_ASSERT_THROWS_NOTHING( read_image( image, "rgb2x2.png" ) );
    TS_ASSERT_EQUALS( image.cols(), 2 );
    TS_ASSERT_EQUALS( image.rows(), 2 );
    TS_ASSERT_EQUALS( image.planes(), 1 );
    TS_ASSERT_EQUALS( image(0,0).r(), 128 );
    TS_ASSERT_EQUALS( image(0,0).g(), 128 );
    TS_ASSERT_EQUALS( image(0,0).b(), 128 );
    TS_ASSERT_EQUALS( image(1,0).r(), 85 );
    TS_ASSERT_EQUALS( image(1,0).g(), 0 );
    TS_ASSERT_EQUALS( image(1,0).b(), 0 );
    TS_ASSERT_EQUALS( image(0,1).r(), 0 );
    TS_ASSERT_EQUALS( image(0,1).g(), 170 );
    TS_ASSERT_EQUALS( image(0,1).b(), 0 );
    TS_ASSERT_EQUALS( image(1,1).r(), 0 );
    TS_ASSERT_EQUALS( image(1,1).g(), 0 );
    TS_ASSERT_EQUALS( image(1,1).b(), 255 );
#endif
  }

  void test_read_image_rgb_png_float() {
#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1
    ImageView<PixelRGB<float> > image;
    TS_ASSERT_THROWS_NOTHING( read_image( image, "rgb2x2.png" ) );
    TS_ASSERT_EQUALS( image.cols(), 2 );
    TS_ASSERT_EQUALS( image.rows(), 2 );
    TS_ASSERT_EQUALS( image.planes(), 1 );
    TS_ASSERT_DELTA( image(0,0).r(), 0.50196, 1e-5 );
    TS_ASSERT_DELTA( image(0,0).g(), 0.50196, 1e-5 );
    TS_ASSERT_DELTA( image(0,0).b(), 0.50196, 1e-5 );
    TS_ASSERT_DELTA( image(1,0).r(), 0.33333, 1e-5 );
    TS_ASSERT_EQUALS( image(1,0).g(), 0 );
    TS_ASSERT_EQUALS( image(1,0).b(), 0 );
    TS_ASSERT_EQUALS( image(0,1).r(), 0 );
    TS_ASSERT_DELTA( image(0,1).g(), 0.66667, 1e-5 );
    TS_ASSERT_EQUALS( image(0,1).b(), 0 );
    TS_ASSERT_EQUALS( image(1,1).r(), 0 );
    TS_ASSERT_EQUALS( image(1,1).g(), 0 );
    TS_ASSERT_DELTA( image(1,1).b(), 1.0, 1e-5 );
#endif
  }

  void test_read_image_rgb_tif_uint8() {
#if defined(VW_HAVE_PKG_TIFF) && VW_HAVE_PKG_TIFF==1
    ImageView<PixelRGB<uint8> > image;
    TS_ASSERT_THROWS_NOTHING( read_image( image, "rgb2x2.tif" ) );
    TS_ASSERT_EQUALS( image.cols(), 2 );
    TS_ASSERT_EQUALS( image.rows(), 2 );
    TS_ASSERT_EQUALS( image.planes(), 1 );
    TS_ASSERT_EQUALS( image(0,0).r(), 128 );
    TS_ASSERT_EQUALS( image(0,0).g(), 128 );
    TS_ASSERT_EQUALS( image(0,0).b(), 128 );
    TS_ASSERT_EQUALS( image(1,0).r(), 85 );
    TS_ASSERT_EQUALS( image(1,0).g(), 0 );
    TS_ASSERT_EQUALS( image(1,0).b(), 0 );
    TS_ASSERT_EQUALS( image(0,1).r(), 0 );
    TS_ASSERT_EQUALS( image(0,1).g(), 170 );
    TS_ASSERT_EQUALS( image(0,1).b(), 0 );
    TS_ASSERT_EQUALS( image(1,1).r(), 0 );
    TS_ASSERT_EQUALS( image(1,1).g(), 0 );
    TS_ASSERT_EQUALS( image(1,1).b(), 255 );
#endif
  }

  void test_read_image_rgb_tif_float() {
#if defined(VW_HAVE_PKG_TIFF) && VW_HAVE_PKG_TIFF==1
    ImageView<PixelRGB<float> > image;
    TS_ASSERT_THROWS_NOTHING( read_image( image, "rgb2x2.tif" ) );
    TS_ASSERT_EQUALS( image.cols(), 2 );
    TS_ASSERT_EQUALS( image.rows(), 2 );
    TS_ASSERT_EQUALS( image.planes(), 1 );
    TS_ASSERT_DELTA( image(0,0).r(), 0.50196, 1e-5 );
    TS_ASSERT_DELTA( image(0,0).g(), 0.50196, 1e-5 );
    TS_ASSERT_DELTA( image(0,0).b(), 0.50196, 1e-5 );
    TS_ASSERT_DELTA( image(1,0).r(), 0.33333, 1e-5 );
    TS_ASSERT_EQUALS( image(1,0).g(), 0 );
    TS_ASSERT_EQUALS( image(1,0).b(), 0 );
    TS_ASSERT_EQUALS( image(0,1).r(), 0 );
    TS_ASSERT_DELTA( image(0,1).g(), 0.66667, 1e-5 );
    TS_ASSERT_EQUALS( image(0,1).b(), 0 );
    TS_ASSERT_EQUALS( image(1,1).r(), 0 );
    TS_ASSERT_EQUALS( image(1,1).g(), 0 );
    TS_ASSERT_DELTA( image(1,1).b(), 1.0, 1e-5 );
#endif
  }

}; // class TestDiskImageResource
