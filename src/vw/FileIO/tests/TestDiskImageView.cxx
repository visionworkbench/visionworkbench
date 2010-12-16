// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>
#include <vw/FileIO/DiskImageView.h>

#include <vw/config.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageMath.h>

using namespace vw;

#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1
TEST( DiskImageView, Construction ) {

  uint64 misses_before = vw_system_cache().misses();

  // These contortions ensure that a copy of a DiskImageView will
  // still work even after the original has gone out of scope.
  boost::shared_ptr<DiskImageView<PixelRGB<uint8> > > v1( new DiskImageView<PixelRGB<uint8> >( TEST_SRCDIR"/rgb2x2.png" ) );
  boost::shared_ptr<DiskImageView<PixelRGB<uint8> > > v2( new DiskImageView<PixelRGB<uint8> >( *v1 ) );
  v1.reset();
  ImageView<PixelRGB<uint8> > image = *v2;

  // Expect that DiskImageView actually wrote something to the cache
  EXPECT_EQ( 1, vw_system_cache().misses() - misses_before );

  ASSERT_EQ( image.cols(), 2 );
  ASSERT_EQ( image.rows(), 2 );
  ASSERT_EQ( image.planes(), 1 );
  EXPECT_EQ( image(0,0).r(), 128 );
  EXPECT_EQ( image(0,0).g(), 128 );
  EXPECT_EQ( image(0,0).b(), 128 );
  EXPECT_EQ( image(1,0).r(), 85 );
  EXPECT_EQ( image(1,0).g(), 0 );
  EXPECT_EQ( image(1,0).b(), 0 );
  EXPECT_EQ( image(0,1).r(), 0 );
  EXPECT_EQ( image(0,1).g(), 170 );
  EXPECT_EQ( image(0,1).b(), 0 );
  EXPECT_EQ( image(1,1).r(), 0 );
  EXPECT_EQ( image(1,1).g(), 0 );
  EXPECT_EQ( image(1,1).b(), 255 );
}
#endif

#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1 &&((defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1) || (defined(VW_HAVE_PKG_TIFF) && VW_HAVE_PKG_TIFF==1))
TEST( DiskCacheImageView, Construction ) {
  ImageView<PixelRGB<uint8> > orig_image;
  ASSERT_NO_THROW( read_image( orig_image, TEST_SRCDIR"/rgb2x2.png" ) );

  uint64 misses_before = vw_system_cache().misses();
  uint64 hits_before   = vw_system_cache().hits();

  DiskCacheImageView<PixelRGB<uint8> > image = orig_image;

  ASSERT_EQ( image.cols(), 2 );
  ASSERT_EQ( image.rows(), 2 );
  ASSERT_EQ( image.planes(), 1 );
  EXPECT_EQ( image(0,0).r(), 128 ); // miss
  EXPECT_EQ( image(0,0).g(), 128 ); // hit
  EXPECT_EQ( image(0,0).b(), 128 );
  EXPECT_EQ( image(1,0).r(), 85 );
  EXPECT_EQ( image(1,0).g(), 0 );
  EXPECT_EQ( image(1,0).b(), 0 );
  EXPECT_EQ( image(0,1).r(), 0 );
  EXPECT_EQ( image(0,1).g(), 170 );
  EXPECT_EQ( image(0,1).b(), 0 );
  EXPECT_EQ( image(1,1).r(), 0 );
  EXPECT_EQ( image(1,1).g(), 0 );
  EXPECT_EQ( image(1,1).b(), 255 ); // hit 11

  // Expect that DiskImageView actually wrote something to cache and
  // then used it.
  EXPECT_EQ( 1,  vw_system_cache().misses() - misses_before );
  EXPECT_EQ( 11, vw_system_cache().hits() - hits_before );

  // Test copying
  {
    DiskCacheImageView<PixelRGB<uint8> > cache2 = image;
    DiskCacheImageView<PixelRGB<uint8> > cache3 = cache2;
  }

  // Test re-assignment
  image = orig_image + PixelRGB<uint8>(10,10,10);

  // Test re-assignment after a copy
  DiskCacheImageView<PixelRGB<uint8> > cache4 = image;
  image = orig_image + PixelRGB<uint8>(20,20,20);

}
#endif
