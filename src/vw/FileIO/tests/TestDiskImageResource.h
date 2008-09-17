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

// TestDiskImageResource.h
//#define CXXTEST_ABORT_TEST_ON_FAIL
#include <cxxtest/TestSuite.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageResource_internal.h>

#include <vw/config.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/ImageView.h>

using namespace vw;
using namespace vw::internal;

using std::string;
using std::set;

template <class PixelT>
static void test_extension(string fn)
{
  ImageView<PixelT> img1(1,1), img2;

  try {
    write_image(fn, img1);
    read_image(img2, fn);
  } catch (vw::NoImplErr &e) {
    // this doesn't represent a test failure, just a lack of test coverage
    TS_WARN(string("Failed to test ") + fn + " : " + e.what());
    return;
  }

  TS_ASSERT_EQUALS(img1.cols()        , img2.cols());
  TS_ASSERT_EQUALS(img1.rows()        , img2.rows());
  TS_ASSERT_EQUALS(img1.planes()      , img2.planes());
  TS_ASSERT_EQUALS(img1.channels()    , img2.channels());
  TS_ASSERT_EQUALS(img1.channel_type(), img2.channel_type());
}

class TestDiskImageResource : public CxxTest::TestSuite
{
public:

  void test_write_read_view() {
    set<string> exclude;
    const char *ex_list[] = {"img", "lbl", "pds"}; // skip the ro PDS formats
    exclude.insert(ex_list, ex_list+3);

    foreach_ext("rwtest", test_extension<PixelRGB<float>  >, exclude);
    foreach_ext("rwtest", test_extension<PixelRGB<uint8>  >, exclude);
    foreach_ext("rwtest", test_extension<PixelRGBA<uint8> >, exclude);
    foreach_ext("rwtest", test_extension<uint8>, exclude);
    foreach_ext("rwtest", test_extension<float>, exclude);
  }

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

  void test_write_image_rgb_png_uint8() {
#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1
    ImageView<PixelRGB<uint8> > image(2,2);
    image(0,0) = PixelRGB<uint8>(120,240,180);
    image(1,0) = PixelRGB<uint8>(36,89,79);
    image(0,1) = PixelRGB<uint8>(190,25,34);
    image(1,1) = PixelRGB<uint8>(23,13,189);
    TS_ASSERT_THROWS_NOTHING( write_image( "tmp.png", image ) );
    ImageView<PixelRGB<uint8> > image2;
    TS_ASSERT_THROWS_NOTHING( read_image( image2, "tmp.png" ) );
    TS_ASSERT_EQUALS( image2.cols(), image.cols() );
    TS_ASSERT_EQUALS( image2.rows(), image.rows() );
    TS_ASSERT_EQUALS( image2.planes(), image.planes() );
    TS_ASSERT_EQUALS( image2(0,0).r(), image(0,0).r() );
    TS_ASSERT_EQUALS( image2(0,0).g(), image(0,0).g() );
    TS_ASSERT_EQUALS( image2(0,0).b(), image(0,0).b() );
    TS_ASSERT_EQUALS( image2(1,0).r(), image(1,0).r() );
    TS_ASSERT_EQUALS( image2(1,0).g(), image(1,0).g() );
    TS_ASSERT_EQUALS( image2(1,0).b(), image(1,0).b() );
    TS_ASSERT_EQUALS( image2(0,1).r(), image(0,1).r() );
    TS_ASSERT_EQUALS( image2(0,1).g(), image(0,1).g() );
    TS_ASSERT_EQUALS( image2(0,1).b(), image(0,1).b() );
    TS_ASSERT_EQUALS( image2(1,1).r(), image(1,1).r() );
    TS_ASSERT_EQUALS( image2(1,1).g(), image(1,1).g() );
    TS_ASSERT_EQUALS( image2(1,1).b(), image(1,1).b() );
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

  void test_write_image_rgb_png_float() {
#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1
    ImageView<PixelRGB<float> > image(2,2);
    image(0,0) = PixelRGB<float>(0.7,0.3,0.2);
    image(1,0) = PixelRGB<float>(0.8,0.1,0.2);
    image(0,1) = PixelRGB<float>(0.2,0.7,0.4);
    image(1,1) = PixelRGB<float>(0.3,0.4,0.5);
    TS_ASSERT_THROWS_NOTHING( write_image( "tmp.png", image ) );
    ImageView<PixelRGB<float> > image2;
    TS_ASSERT_THROWS_NOTHING( read_image( image2, "tmp.png" ) );
    TS_ASSERT_EQUALS( image2.cols(), image.cols() );
    TS_ASSERT_EQUALS( image2.rows(), image.rows() );
    TS_ASSERT_EQUALS( image2.planes(), image.planes() );
    // PNG is an 8-bit format, so values may be corrupted
    TS_ASSERT_DELTA( image2(0,0).r(), image(0,0).r(), 1e-2 );
    TS_ASSERT_DELTA( image2(0,0).g(), image(0,0).g(), 1e-2 );
    TS_ASSERT_DELTA( image2(0,0).b(), image(0,0).b(), 1e-2 );
    TS_ASSERT_DELTA( image2(1,0).r(), image(1,0).r(), 1e-2 );
    TS_ASSERT_DELTA( image2(1,0).g(), image(1,0).g(), 1e-2 );
    TS_ASSERT_DELTA( image2(1,0).b(), image(1,0).b(), 1e-2 );
    TS_ASSERT_DELTA( image2(0,1).r(), image(0,1).r(), 1e-2 );
    TS_ASSERT_DELTA( image2(0,1).g(), image(0,1).g(), 1e-2 );
    TS_ASSERT_DELTA( image2(0,1).b(), image(0,1).b(), 1e-2 );
    TS_ASSERT_DELTA( image2(1,1).r(), image(1,1).r(), 1e-2 );
    TS_ASSERT_DELTA( image2(1,1).g(), image(1,1).g(), 1e-2 );
    TS_ASSERT_DELTA( image2(1,1).b(), image(1,1).b(), 1e-2 );
#endif
  }

  void test_read_image_rgb_tif_uint8() {
#if (defined(VW_HAVE_PKG_TIFF) && VW_HAVE_PKG_TIFF==1) || (defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1)
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

  void test_write_image_rgb_tif_uint8() {
#if (defined(VW_HAVE_PKG_TIFF) && VW_HAVE_PKG_TIFF==1) || (defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1)
    ImageView<PixelRGB<uint8> > image(2,2);
    image(0,0) = PixelRGB<uint8>(120,240,180);
    image(1,0) = PixelRGB<uint8>(36,89,79);
    image(0,1) = PixelRGB<uint8>(190,25,34);
    image(1,1) = PixelRGB<uint8>(23,13,189);
    TS_ASSERT_THROWS_NOTHING( write_image( "tmp.tif", image ) );
    ImageView<PixelRGB<uint8> > image2;
    TS_ASSERT_THROWS_NOTHING( read_image( image2, "tmp.tif" ) );
    TS_ASSERT_EQUALS( image2.cols(), image.cols() );
    TS_ASSERT_EQUALS( image2.rows(), image.rows() );
    TS_ASSERT_EQUALS( image2.planes(), image.planes() );
    TS_ASSERT_EQUALS( image2(0,0).r(), image(0,0).r() );
    TS_ASSERT_EQUALS( image2(0,0).g(), image(0,0).g() );
    TS_ASSERT_EQUALS( image2(0,0).b(), image(0,0).b() );
    TS_ASSERT_EQUALS( image2(1,0).r(), image(1,0).r() );
    TS_ASSERT_EQUALS( image2(1,0).g(), image(1,0).g() );
    TS_ASSERT_EQUALS( image2(1,0).b(), image(1,0).b() );
    TS_ASSERT_EQUALS( image2(0,1).r(), image(0,1).r() );
    TS_ASSERT_EQUALS( image2(0,1).g(), image(0,1).g() );
    TS_ASSERT_EQUALS( image2(0,1).b(), image(0,1).b() );
    TS_ASSERT_EQUALS( image2(1,1).r(), image(1,1).r() );
    TS_ASSERT_EQUALS( image2(1,1).g(), image(1,1).g() );
    TS_ASSERT_EQUALS( image2(1,1).b(), image(1,1).b() );
#endif
  }

  void test_read_image_rgb_tif_float() {
#if (defined(VW_HAVE_PKG_TIFF) && VW_HAVE_PKG_TIFF==1) || (defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1)
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

  void test_write_image_rgb_tif_float() {
#if (defined(VW_HAVE_PKG_TIFF) && VW_HAVE_PKG_TIFF==1) || (defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1)
    ImageView<PixelRGB<float> > image(2,2);
    image(0,0) = PixelRGB<float>(0.7,0.3,0.2);
    image(1,0) = PixelRGB<float>(0.8,0.1,0.2);
    image(0,1) = PixelRGB<float>(0.2,0.7,0.4);
    image(1,1) = PixelRGB<float>(0.3,0.4,0.5);
    TS_ASSERT_THROWS_NOTHING( write_image( "tmp.tif", image ) );
    ImageView<PixelRGB<float> > image2;
    TS_ASSERT_THROWS_NOTHING( read_image( image2, "tmp.tif" ) );
    TS_ASSERT_EQUALS( image2.cols(), image.cols() );
    TS_ASSERT_EQUALS( image2.rows(), image.rows() );
    TS_ASSERT_EQUALS( image2.planes(), image.planes() );
    // TIFF supports floating-point data, so values should
    // be reasonably precise.
    TS_ASSERT_DELTA( image2(0,0).r(), image(0,0).r(), 1e-8 );
    TS_ASSERT_DELTA( image2(0,0).g(), image(0,0).g(), 1e-8 );
    TS_ASSERT_DELTA( image2(0,0).b(), image(0,0).b(), 1e-8 );
    TS_ASSERT_DELTA( image2(1,0).r(), image(1,0).r(), 1e-8 );
    TS_ASSERT_DELTA( image2(1,0).g(), image(1,0).g(), 1e-8 );
    TS_ASSERT_DELTA( image2(1,0).b(), image(1,0).b(), 1e-8 );
    TS_ASSERT_DELTA( image2(0,1).r(), image(0,1).r(), 1e-8 );
    TS_ASSERT_DELTA( image2(0,1).g(), image(0,1).g(), 1e-8 );
    TS_ASSERT_DELTA( image2(0,1).b(), image(0,1).b(), 1e-8 );
    TS_ASSERT_DELTA( image2(1,1).r(), image(1,1).r(), 1e-8 );
    TS_ASSERT_DELTA( image2(1,1).g(), image(1,1).g(), 1e-8 );
    TS_ASSERT_DELTA( image2(1,1).b(), image(1,1).b(), 1e-8 );
#endif
  }

  void test_read_image_rgb_jpeg_uint8() {
#if defined(VW_HAVE_PKG_JPEG) && VW_HAVE_PKG_JPEG==1
    ImageView<PixelRGB<uint8> > image;
    TS_ASSERT_THROWS_NOTHING( read_image( image, "rgb2x2.jpg" ) );
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
    TS_ASSERT_DELTA( image(0,1).b(), 0, 1 );
    TS_ASSERT_EQUALS( image(1,1).r(), 0 );
    TS_ASSERT_EQUALS( image(1,1).g(), 0 );
    TS_ASSERT_DELTA( image(1,1).b(), 255, 1 );
#endif
  }

  void test_read_image_rgb_jpeg_float() {
#if defined(VW_HAVE_PKG_JPEG) && VW_HAVE_PKG_JPEG==1
    ImageView<PixelRGB<float> > image;
    TS_ASSERT_THROWS_NOTHING( read_image( image, "rgb2x2.jpg" ) );
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
    TS_ASSERT_DELTA( image(0,1).b(), 0, 0.004 );
    TS_ASSERT_EQUALS( image(1,1).r(), 0 );
    TS_ASSERT_EQUALS( image(1,1).g(), 0 );
    TS_ASSERT_DELTA( image(1,1).b(), 1.0, 0.004 );
#endif
  }

}; // class TestDiskImageResource
