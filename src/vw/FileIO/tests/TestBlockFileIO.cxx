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


// TestBlockFileIO.h
#include <gtest/gtest_VW.h>
#include <vw/config.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Image/ImageIO.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Math/BBox.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageView.h>
#include <test/Helpers.h>

#include <string>

#include <boost/scoped_ptr.hpp>
#include <boost/filesystem/path_traits.hpp>
namespace fs = boost::filesystem;

using namespace vw;
using namespace vw::test;

static void test_read_crop(const std::string& input, const UnlinkName& output) {
  boost::scoped_ptr<DiskImageResource> dir;
  ASSERT_NO_THROW( dir.reset(DiskImageResource::open( input ) ) );

  ImageView<PixelRGB<uint8> > image;
  read_image( image, *dir, BBox2i(100,100,100,100) );
  write_image( output, image );
}

TEST( BlockFileIO, PNG_Crop ) {
  test_read_crop(TEST_SRCDIR"/mural.png", "cropped.mural.png");
}

#if defined(VW_HAVE_PKG_JPEG) && VW_HAVE_PKG_JPEG==1
TEST( BlockFileIO, JPG_Crop ) {
  test_read_crop(TEST_SRCDIR"/mural.jpg", "cropped.mural.jpg");
}
#endif

TEST( BlockFileIO, TIF_Post_Crop ) {
  UnlinkName fn1("mural.tif");
  UnlinkName fn2("cropped.mural.tif");

  boost::scoped_ptr<DiskImageResource> dir;
  ASSERT_NO_THROW( dir.reset(DiskImageResource::open( TEST_SRCDIR"/mural.png" ) ) );

  ImageView<PixelRGB<uint8> > image;
  ASSERT_NO_THROW( read_image( image, *dir ) );
  write_image( fn1, image );
  DiskImageView<PixelRGB<uint8> > div( fn1 );
  ImageView<PixelRGB<uint8> > result = crop(div,100,100,100,100);
  write_image(fn2, result );
}
