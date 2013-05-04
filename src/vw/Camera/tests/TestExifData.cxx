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


// TestExifData
#include <gtest/gtest_VW.h>

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
