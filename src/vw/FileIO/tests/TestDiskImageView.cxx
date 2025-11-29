// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#include <gtest/gtest_VW.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/FileIO/DiskCacheImageView.h>
#include <vw/FileIO/DiskImageResourceRaw.h>
#include <vw/FileIO/DiskImageResourceGDAL.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/ImageIO.h>

#include <vw/config.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageMath.h>

using namespace vw;

#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1
TEST( DiskImageView, Construction ) {

  // These contortions ensure that a copy of a DiskImageView will
  // still work even after the original has gone out of scope.
  boost::shared_ptr<DiskImageView<PixelRGB<uint8> > > v1( new DiskImageView<PixelRGB<uint8> >( TEST_SRCDIR"/rgb2x2.png" ) );
  boost::shared_ptr<DiskImageView<PixelRGB<uint8> > > v2( new DiskImageView<PixelRGB<uint8> >( *v1 ) );
  v1.reset();
  ImageView<PixelRGB<uint8> > image = *v2;

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

  DiskCacheImageView<PixelRGB<uint8> > image = orig_image;

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

// This DiskImageResource test is located in this file so that we can
//  use DiskImageView to help us test it.
TEST( DiskImageResource , Raw ) {

  // Create a read-only image reader for the sample image
  ImageFormat format;
  format.cols = 300;
  format.rows = 300;
  format.planes = 1;
  format.pixel_format = VW_PIXEL_GRAY;
  format.channel_type = VW_CHANNEL_UINT8; // Open read-only with a tiny block size
  DiskImageResourceRaw resource("sample.BIL", format, true, Vector2i(128, 128));
  DiskImageView<PixelGray<unsigned char> > view(resource);
  
  // Check the image values
  EXPECT_EQ(view.rows(), format.rows);
  EXPECT_EQ(64, view(82, 109));
  
  // Check that we can write the image using the block writer.
  // To ensure the image is completely written, force gdal_resource
  // to go out of scope.
  {
    DiskImageResourceGDAL gdal_resource("bil_test.tif", view.impl().format());
    ASSERT_NO_THROW(block_write_image(gdal_resource, view));
  }

  // Now read it back
  boost::shared_ptr<DiskImageResource> gdal_resource2;
  ASSERT_NO_THROW(gdal_resource2.reset(DiskImageResource::open("bil_test.tif")));

  // Test that the default loader behaves as expected. It should fail
  // because sample.DIM is missing. 
  boost::shared_ptr<DiskImageResource> generic_resource_ptr;
  EXPECT_THROW(generic_resource_ptr.reset(DiskImageResource::open("sample.BIL")), vw::ArgumentErr); 
}


