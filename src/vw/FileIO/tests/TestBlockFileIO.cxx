// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestBlockFileIO.h
#include <gtest/gtest.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageView.h>
#include <test/Helpers.h>

#include <boost/scoped_ptr.hpp>
#include <string>

using namespace vw;
using namespace vw::test;

#include <boost/filesystem/path.hpp>
namespace fs = boost::filesystem;

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
