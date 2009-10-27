// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestLinearPushbroomModel.h
#include <cxxtest/TestSuite.h>

#include <cstdio>
#include <vw/Camera/Exif.h>

using namespace vw;
using namespace vw::camera;

class TestExifData : public CxxTest::TestSuite
{
public:
  void test_exif_data()
  {
    ExifData exif_data;
    exif_data.import_data(TEST_SRCDIR"/test.jpg");
    //exif_data.print_debug();
    {
      double exposure;
      bool found = exif_data.get_tag_value(EXIF_ExposureTime, exposure);
      TS_ASSERT(found);
      TS_ASSERT_DELTA(exposure, 0.0125, 1e-8);
    }
    {
      double fstop;
      bool found = exif_data.get_tag_value(EXIF_FNumber, fstop);
      TS_ASSERT(found);
      TS_ASSERT_DELTA(fstop, 5.8, 1e-8);
    }
    {
      int thumbnail_offset;
      int thumbnail_length;
      int exif_location;
      bool found = exif_data.get_tag_value(EXIF_ThumbnailOffset, thumbnail_offset);
      TS_ASSERT(found);
      found = exif_data.get_tag_value(EXIF_ThumbnailLength, thumbnail_length);
      TS_ASSERT(found);
      exif_location = exif_data.get_exif_location();

      printf("exif location= %d, offset=%d, length=%d\n",
             exif_location, thumbnail_offset, thumbnail_length);

    }
  }
};
