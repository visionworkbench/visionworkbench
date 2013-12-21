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
#include <test/Helpers.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/ImageView.h>
#include <vw/FileIO/DiskImageResource.h>

using namespace vw;
using namespace vw::test;

TEST( Endianness, SameEndPNG ) {
  // this tests to make sure VW can read images back in in the same
  // endianness it wrote them out in

  UnlinkName fn("test-png16.png");

  typedef PixelRGBA<uint16> P16;

  // numbers chosen so endianness matters
  static const P16 one   = P16(0x0102, 0x0304, 0x0506, 0xffff),
    two   = P16(0x0708, 0x0910, 0x1213, 0xffff),
    three = P16(0x1415, 0x1617, 0x1819, 0xffff),
    four  = P16(0x2021, 0x2324, 0x2526, 0xffff);

  ImageView<P16> img3(2,2), img4;
  img3(0,0) = pixel_cast<P16>(one);
  img3(0,1) = pixel_cast<P16>(two);
  img3(1,0) = pixel_cast<P16>(three);
  img3(1,1) = pixel_cast<P16>(four);

  EXPECT_EQ(img3(0,0), pixel_cast<P16>(  one));
  EXPECT_EQ(img3(0,1), pixel_cast<P16>(  two));
  EXPECT_EQ(img3(1,0), pixel_cast<P16>(three));
  EXPECT_EQ(img3(1,1), pixel_cast<P16>( four));

  write_image(fn, img3);

  read_image(img4, fn);

  EXPECT_EQ(img4(0,0), pixel_cast<P16>(  one));
  EXPECT_EQ(img4(0,1), pixel_cast<P16>(  two));
  EXPECT_EQ(img4(1,0), pixel_cast<P16>(three));
  EXPECT_EQ(img4(1,1), pixel_cast<P16>( four));
}

TEST( Endianness, DifferentEndPNG ) {
  // This tests to make sure a ground-truth png is read in with the right
  // endianness.

  typedef PixelGray<uint16> P16;
  ImageView<P16> img;

  read_image(img, TEST_SRCDIR"/png16.png");

  EXPECT_EQ(img(0,0), P16(0x0102));
  EXPECT_EQ(img(1,0), P16(0x0304));
  EXPECT_EQ(img(0,1), P16(0x0506));
  EXPECT_EQ(img(1,1), P16(0x0708));
}
