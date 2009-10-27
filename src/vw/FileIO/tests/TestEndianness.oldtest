// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <cxxtest/TestSuite.h>
#include <vw/Image.h>
#include <vw/FileIO.h>

using namespace vw;

class TestEndianness : public CxxTest::TestSuite
{public:

  void testSameEndianPNG()
  {
    typedef PixelRGBA<uint16> P16;

    // this tests to make sure VW can read images back in in the same
    // endianness it wrote them out in

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

    TS_ASSERT_EQUALS(img3(0,0), pixel_cast<P16>(  one));
    TS_ASSERT_EQUALS(img3(0,1), pixel_cast<P16>(  two));
    TS_ASSERT_EQUALS(img3(1,0), pixel_cast<P16>(three));
    TS_ASSERT_EQUALS(img3(1,1), pixel_cast<P16>( four));

    write_image(TEST_SRCDIR"/test-png16.png", img3);

    read_image(img4, TEST_SRCDIR"/test-png16.png");

    TS_ASSERT_EQUALS(img4(0,0), pixel_cast<P16>(  one));
    TS_ASSERT_EQUALS(img4(0,1), pixel_cast<P16>(  two));
    TS_ASSERT_EQUALS(img4(1,0), pixel_cast<P16>(three));
    TS_ASSERT_EQUALS(img4(1,1), pixel_cast<P16>( four));
  }

  void testDifferentEndianPNG()
  {
    // This tests to make sure a ground-truth png is read in with the right
    // endianness.

    typedef PixelGray<uint16> P16;
    ImageView<P16> img;

    read_image(img, TEST_SRCDIR"/png16.png");

    TS_ASSERT_EQUALS(img(0,0), P16(0x0102));
    TS_ASSERT_EQUALS(img(1,0), P16(0x0304));
    TS_ASSERT_EQUALS(img(0,1), P16(0x0506));
    TS_ASSERT_EQUALS(img(1,1), P16(0x0708));

  }

};
