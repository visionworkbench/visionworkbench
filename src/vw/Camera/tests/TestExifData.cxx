// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestExifData
#include <gtest/gtest.h>

#include <cstdio>
#include <vw/Camera/Exif.h>

using namespace vw;
using namespace vw::camera;

TEST( Exif, ReadData ) {

  ExifData exif_data;
  exif_data.import_data(TEST_SRCDIR"/test.jpg");
  {
    double exposure;
    bool found = exif_data.get_tag_value(EXIF_ExposureTime, exposure);
    EXPECT_TRUE( found );
    EXPECT_NEAR( exposure, 0.0125, 1e-8 );
  }
  {
    double fstop;
    bool found = exif_data.get_tag_value(EXIF_FNumber, fstop);
    EXPECT_TRUE( found );
    EXPECT_NEAR( fstop, 5.8, 1e-8 );
  }

  {
    int thumbnail_offset;
    int thumbnail_length;
    int exif_location;
    bool found = exif_data.get_tag_value(EXIF_ThumbnailOffset, thumbnail_offset);
    EXPECT_TRUE(found);
    found = exif_data.get_tag_value(EXIF_ThumbnailLength, thumbnail_length);
    EXPECT_TRUE(found);
    exif_location = exif_data.get_exif_location();

    printf("exif location= %d, offset=%d, length=%d\n",
           exif_location, thumbnail_offset, thumbnail_length);

  }
}
