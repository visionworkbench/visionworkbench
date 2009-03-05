// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestDiskImageResource.h
#define CXXTEST_ABORT_TEST_ON_FAIL
#include <cxxtest/TestSuite.h>
#include <vw/FileIO/DiskImageView.h>

#include <vw/config.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageMath.h>

using namespace vw;

class TestDiskImageView : public CxxTest::TestSuite
{
public:

  void test_disk_image_view() {
#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1
    // These contortions ensure that a copy of a DiskImageView will
    // still work even after the original has gone out of scope.
    boost::shared_ptr<DiskImageView<PixelRGB<uint8> > > v1( new DiskImageView<PixelRGB<uint8> >( "rgb2x2.png" ) );
    boost::shared_ptr<DiskImageView<PixelRGB<uint8> > > v2( new DiskImageView<PixelRGB<uint8> >( *v1 ) );
    v1.reset();
    ImageView<PixelRGB<uint8> > image = *v2;
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

  void test_disk_cache_image_view() {
#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1 &&((defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1) || (defined(VW_HAVE_PKG_TIFF) && VW_HAVE_PKG_TIFF==1)) 
    ImageView<PixelRGB<uint8> > orig_image;
    TS_ASSERT_THROWS_NOTHING( read_image( orig_image, "rgb2x2.png" ) );

    DiskCacheImageView<PixelRGB<uint8> > image = orig_image;

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

    // Test copying
    if (1) {
      DiskCacheImageView<PixelRGB<uint8> > cache2 = image;
      DiskCacheImageView<PixelRGB<uint8> > cache3 = cache2;
    }

    // Test re-assignment
    image = orig_image + PixelRGB<uint8>(10,10,10);

    // Test re-assignment after a copy
    DiskCacheImageView<PixelRGB<uint8> > cache4 = image;
    image = orig_image + PixelRGB<uint8>(20,20,20);
#endif
  }
}; // class TestDiskImageResource
