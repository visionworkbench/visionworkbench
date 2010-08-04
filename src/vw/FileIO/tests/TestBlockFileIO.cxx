// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestBlockFileIO.h
#include <gtest/gtest.h>
#include <string>
#include <vw/FileIO/DiskImageResource.h>

#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/FileIO/DiskImageView.h>
using namespace vw;

#include <boost/filesystem/path.hpp>
namespace fs = boost::filesystem;


static void test_read_crop(const char *fn) {
  DiskImageResource *dir = 0;
  ASSERT_NO_THROW( dir = DiskImageResource::open( fn ) );

  fs::path crop_name(fn);
  crop_name = crop_name.branch_path() / (std::string("cropped.") + crop_name.leaf());

  ImageView<PixelRGB<uint8> > image;
  read_image( image, *dir, BBox2i(100,100,100,100) );
  write_image( crop_name.string(), image );
 }

TEST( BlockFileIO, PNG_Crop ) {
  test_read_crop(TEST_SRCDIR"/mural.png");
}

TEST( BlockFileIO, JPG_Crop ) {
  test_read_crop(TEST_SRCDIR"/mural.jpg");
}

TEST( BlockFileIO, TIF_Post_Crop ) {
  DiskImageResource *dir = 0;
  ASSERT_NO_THROW( dir = DiskImageResource::open( TEST_SRCDIR"/mural.png" ) );

  ImageView<PixelRGB<uint8> > image;
  ASSERT_NO_THROW( read_image( image, *dir ) );
  write_image( "mural.tif", image );
  DiskImageView<PixelRGB<uint8> > div( "mural.tif" );
  ImageView<PixelRGB<uint8> > result = crop(div,100,100,100,100);
  write_image( "cropped.mural.tif", result );
}
