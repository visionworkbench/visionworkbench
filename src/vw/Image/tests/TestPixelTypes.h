// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__

// TestPixelTypes.h
#include <cxxtest/TestSuite.h>

#include <vw/Image/PixelTypes.h>

using namespace std;
using namespace vw;

class TestPixelTypes : public CxxTest::TestSuite
{
public:

  void test_channel_types()
  {
    TS_ASSERT_EQUALS( ChannelRange<uint8>::max(), 255 );
    TS_ASSERT_EQUALS( ChannelRange<uint16>::max(), 65535 );
    TS_ASSERT_EQUALS( ChannelRange<float32>::max(), 1.0 );
    TS_ASSERT_EQUALS( ChannelRange<float64>::max(), 1.0 );

    TS_ASSERT_EQUALS( channel_cast<uint16>(uint8(255)), 255 );
    TS_ASSERT_EQUALS( channel_cast_rescale<uint16>(uint8(255)), 65535 );
    TS_ASSERT_EQUALS( channel_cast<uint8>(float32(17.0)), 17 );
    TS_ASSERT_EQUALS( channel_cast_rescale<uint8>(float32(0.333334)), 85 );
  }

  void test_pixel_gray()
  {
    // Test default-construction and size with all supported channel types
    { PixelGray<int8> p; TS_ASSERT( p.v()==0 ); TS_ASSERT( sizeof(p)==1 ); }
    { PixelGray<uint8> p; TS_ASSERT( p.v()==0 ); TS_ASSERT( sizeof(p)==1 ); }
    { PixelGray<int16> p; TS_ASSERT( p.v()==0 ); TS_ASSERT( sizeof(p)==2 ); }
    { PixelGray<uint16> p; TS_ASSERT( p.v()==0 ); TS_ASSERT( sizeof(p)==2 ); }
    { PixelGray<int32> p; TS_ASSERT( p.v()==0 ); TS_ASSERT( sizeof(p)==4 ); }
    { PixelGray<uint32> p; TS_ASSERT( p.v()==0 ); TS_ASSERT( sizeof(p)==4 ); }
    { PixelGray<int64> p; TS_ASSERT( p.v()==0 ); TS_ASSERT( sizeof(p)==8 ); }
    { PixelGray<uint64> p; TS_ASSERT( p.v()==0 ); TS_ASSERT( sizeof(p)==8 ); }
    { PixelGray<float32> p; TS_ASSERT( p.v()==0.0 ); TS_ASSERT( sizeof(p)==4 ); }
    { PixelGray<float64> p; TS_ASSERT( p.v()==0.0 ); TS_ASSERT( sizeof(p)==8 ); }
    // Test channel-value-construction and accessors
    { PixelGray<int8> p(1); TS_ASSERT( p.v()==1 ); TS_ASSERT( p[0]==1 ); TS_ASSERT( p(0)==1 ); }
    { PixelGray<uint8> p(1); TS_ASSERT( p.v()==1 ); TS_ASSERT( p[0]==1 ); TS_ASSERT( p(0)==1 ); }
    { PixelGray<int16> p(1); TS_ASSERT( p.v()==1 ); TS_ASSERT( p[0]==1 ); TS_ASSERT( p(0)==1 ); }
    { PixelGray<uint16> p(1); TS_ASSERT( p.v()==1 ); TS_ASSERT( p[0]==1 ); TS_ASSERT( p(0)==1 ); }
    { PixelGray<int32> p(1); TS_ASSERT( p.v()==1 ); TS_ASSERT( p[0]==1 ); TS_ASSERT( p(0)==1 ); }
    { PixelGray<uint32> p(1); TS_ASSERT( p.v()==1 ); TS_ASSERT( p[0]==1 ); TS_ASSERT( p(0)==1 ); }
    { PixelGray<int64> p(1); TS_ASSERT( p.v()==1 ); TS_ASSERT( p[0]==1 ); TS_ASSERT( p(0)==1 ); }
    { PixelGray<uint64> p(1); TS_ASSERT( p.v()==1 ); TS_ASSERT( p[0]==1 ); TS_ASSERT( p(0)==1 ); }
    { PixelGray<float32> p(1.0); TS_ASSERT( p.v()==1.0 ); TS_ASSERT( p[0]==1.0 ); TS_ASSERT( p(0)==1.0 ); }
    { PixelGray<float64> p(1.0); TS_ASSERT( p.v()==1.0 ); TS_ASSERT( p[0]==1.0 ); TS_ASSERT( p(0)==1.0 ); }
  }

  void test_pixel_graya()
  {
    // Test default-construction and size with all supported channel types
    { PixelGrayA<int8> p; TS_ASSERT( p.v()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==2 ); }
    { PixelGrayA<uint8> p; TS_ASSERT( p.v()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==2 ); }
    { PixelGrayA<int16> p; TS_ASSERT( p.v()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==4 ); }
    { PixelGrayA<uint16> p; TS_ASSERT( p.v()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==4 ); }
    { PixelGrayA<int32> p; TS_ASSERT( p.v()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==8 ); }
    { PixelGrayA<uint32> p; TS_ASSERT( p.v()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==8 ); }
    { PixelGrayA<int64> p; TS_ASSERT( p.v()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==16 ); }
    { PixelGrayA<uint64> p; TS_ASSERT( p.v()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==16 ); }
    { PixelGrayA<float32> p; TS_ASSERT( p.v()==0.0 && p.a()==0.0 ); TS_ASSERT( sizeof(p)==8 ); }
    { PixelGrayA<float64> p; TS_ASSERT( p.v()==0.0 && p.a()==0.0 ); TS_ASSERT( sizeof(p)==16 ); }
    // Test channel-value-construction and accessors
    { PixelGrayA<int8> p(1,2); TS_ASSERT( p.v()==1 && p.a()==2 ); TS_ASSERT( p[0]==1 && p[1]==2 ); TS_ASSERT( p(0)==1 && p(1)==2 ); }
    { PixelGrayA<uint8> p(1,2); TS_ASSERT( p.v()==1 && p.a()==2 ); TS_ASSERT( p[0]==1 && p[1]==2 ); TS_ASSERT( p(0)==1 && p(1)==2 ); }
    { PixelGrayA<int16> p(1,2); TS_ASSERT( p.v()==1 && p.a()==2 ); TS_ASSERT( p[0]==1 && p[1]==2 ); TS_ASSERT( p(0)==1 && p(1)==2 ); }
    { PixelGrayA<uint16> p(1,2); TS_ASSERT( p.v()==1 && p.a()==2 ); TS_ASSERT( p[0]==1 && p[1]==2 ); TS_ASSERT( p(0)==1 && p(1)==2 ); }
    { PixelGrayA<int32> p(1,2); TS_ASSERT( p.v()==1 && p.a()==2 ); TS_ASSERT( p[0]==1 && p[1]==2 ); TS_ASSERT( p(0)==1 && p(1)==2 ); }
    { PixelGrayA<uint32> p(1,2); TS_ASSERT( p.v()==1 && p.a()==2 ); TS_ASSERT( p[0]==1 && p[1]==2 ); TS_ASSERT( p(0)==1 && p(1)==2 ); }
    { PixelGrayA<int64> p(1,2); TS_ASSERT( p.v()==1 && p.a()==2 ); TS_ASSERT( p[0]==1 && p[1]==2 ); TS_ASSERT( p(0)==1 && p(1)==2 ); }
    { PixelGrayA<uint64> p(1,2); TS_ASSERT( p.v()==1 && p.a()==2 ); TS_ASSERT( p[0]==1 && p[1]==2 ); TS_ASSERT( p(0)==1 && p(1)==2 ); }
    { PixelGrayA<float32> p(1.0,2.0); TS_ASSERT( p.v()==1.0 && p.a()==2.0 ); TS_ASSERT( p[0]==1.0 && p[1]==2.0 ); TS_ASSERT( p(0)==1.0 && p(1)==2.0 ); }
    { PixelGrayA<float64> p(1.0,2.0); TS_ASSERT( p.v()==1.0 && p.a()==2.0 ); TS_ASSERT( p[0]==1.0 && p[1]==2.0 ); TS_ASSERT( p(0)==1.0 && p(1)==2.0 ); }
  }

  void test_pixel_rgb()
  {
    // Test default-construction and size with all supported channel types
    { PixelRGB<int8> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 ); TS_ASSERT( sizeof(p)==3 ); }
    { PixelRGB<uint8> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 ); TS_ASSERT( sizeof(p)==3 ); }
    { PixelRGB<int16> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 ); TS_ASSERT( sizeof(p)==6 ); }
    { PixelRGB<uint16> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 ); TS_ASSERT( sizeof(p)==6 ); }
    { PixelRGB<int32> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 ); TS_ASSERT( sizeof(p)==12 ); }
    { PixelRGB<uint32> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 ); TS_ASSERT( sizeof(p)==12 ); }
    { PixelRGB<int64> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 ); TS_ASSERT( sizeof(p)==24 ); }
    { PixelRGB<uint64> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 ); TS_ASSERT( sizeof(p)==24 ); }
    { PixelRGB<float32> p; TS_ASSERT( p.r()==0.0 && p.g()==0.0 && p.b()==0.0 ); TS_ASSERT( sizeof(p)==12 ); }
    { PixelRGB<float64> p; TS_ASSERT( p.r()==0.0 && p.g()==0.0 && p.b()==0.0 ); TS_ASSERT( sizeof(p)==24 ); }
    // Test channel-value-construction and accessors
    { PixelRGB<int8> p(1,2,3); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelRGB<uint8> p(1,2,3); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelRGB<int16> p(1,2,3); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelRGB<uint16> p(1,2,3); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelRGB<int32> p(1,2,3); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelRGB<uint32> p(1,2,3); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelRGB<int64> p(1,2,3); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelRGB<uint64> p(1,2,3); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelRGB<float32> p(1.0,2.0,3.0); TS_ASSERT( p.r()==1.0 && p.g()==2.0 && p.b()==3.0 ); TS_ASSERT( p[0]==1.0 && p[1]==2.0 && p[2]==3.0 ); TS_ASSERT( p(0)==1.0 && p(1)==2.0 && p(2)==3.0 ); }
    { PixelRGB<float64> p(1.0,2.0,3.0); TS_ASSERT( p.r()==1.0 && p.g()==2.0 && p.b()==3.0 ); TS_ASSERT( p[0]==1.0 && p[1]==2.0 && p[2]==3.0 ); TS_ASSERT( p(0)==1.0 && p(1)==2.0 && p(2)==3.0 ); }
  }

  void test_pixel_rgba()
  {
    // Test default-construction and size with all supported channel types
    { PixelRGBA<int8> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==4 ); }
    { PixelRGBA<uint8> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==4 ); }
    { PixelRGBA<int16> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==8 ); }
    { PixelRGBA<uint16> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==8 ); }
    { PixelRGBA<int32> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==16 ); }
    { PixelRGBA<uint32> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==16 ); }
    { PixelRGBA<int64> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==32 ); }
    { PixelRGBA<uint64> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==32 ); }
    { PixelRGBA<float32> p; TS_ASSERT( p.r()==0.0 && p.g()==0.0 && p.b()==0.0 && p.a()==0.0 ); TS_ASSERT( sizeof(p)==16 ); }
    { PixelRGBA<float64> p; TS_ASSERT( p.r()==0.0 && p.g()==0.0 && p.b()==0.0 && p.a()==0.0 ); TS_ASSERT( sizeof(p)==32 ); }
    // Test channel-value-construction and accessors
    { PixelRGBA<int8> p(1,2,3,4); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 && p.a()==4 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 && p[3]==4 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 && p(3)==4 ); }
    { PixelRGBA<uint8> p(1,2,3,4); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 && p.a()==4 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 && p[3]==4 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 && p(3)==4 ); }
    { PixelRGBA<int16> p(1,2,3,4); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 && p.a()==4 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 && p[3]==4 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 && p(3)==4 ); }
    { PixelRGBA<uint16> p(1,2,3,4); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 && p.a()==4 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 && p[3]==4 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 && p(3)==4 ); }
    { PixelRGBA<int32> p(1,2,3,4); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 && p.a()==4 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 && p[3]==4 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 && p(3)==4 ); }
    { PixelRGBA<uint32> p(1,2,3,4); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 && p.a()==4 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 && p[3]==4 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 && p(3)==4 ); }
    { PixelRGBA<int64> p(1,2,3,4); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 && p.a()==4 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 && p[3]==4 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 && p(3)==4 ); }
    { PixelRGBA<uint64> p(1,2,3,4); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 && p.a()==4 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 && p[3]==4 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 && p(3)==4 ); }
    { PixelRGBA<float32> p(1.0,2.0,3.0,4.0); TS_ASSERT( p.r()==1.0 && p.g()==2.0 && p.b()==3.0 && p.a()==4.0 ); TS_ASSERT( p[0]==1.0 && p[1]==2.0 && p[2]==3.0 && p[3]==4.0 ); TS_ASSERT( p(0)==1.0 && p(1)==2.0 && p(2)==3.0 && p(3)==4.0 ); }
    { PixelRGBA<float64> p(1.0,2.0,3.0,4.0); TS_ASSERT( p.r()==1.0 && p.g()==2.0 && p.b()==3.0 && p.a()==4.0 ); TS_ASSERT( p[0]==1.0 && p[1]==2.0 && p[2]==3.0 && p[3]==4.0 ); TS_ASSERT( p(0)==1.0 && p(1)==2.0 && p(2)==3.0 && p(3)==4.0 ); }
  }

  void test_rgb_to_gray()
  {
    int val0, val1, val2;
    val0 = val1 = val2 = 40;
   
    // Standard case

    PixelRGB<int8> test_rgb(val0, val1, val2);
    PixelGray<int8> test_gray(test_rgb);
    TS_ASSERT_EQUALS(test_gray.v(),
		     (test_rgb.r() + test_rgb.g() + test_rgb.b())/3 );

    // Case with type incompatibility

    PixelRGB<int16> test_rgb16(val0, val1, val2);
    PixelGray<int8> test_gray8(test_rgb);
    TS_ASSERT_EQUALS(test_gray8.v(),
		     (test_rgb16.r() + test_rgb16.g() + test_rgb16.b())/3 );
  }

  void test_weighted_rgb_to_gray()
  {
    PixelRGB<float> rgbf(0.8,0.4,0.7);
    PixelGray<float> gf = weighted_rgb_to_gray(rgbf);
    TS_ASSERT_DELTA( gf.v(), 0.5530, 1e-4 );

    PixelRGB<uint8> rgbi(180,56,212);
    PixelGray<uint8> gi = weighted_rgb_to_gray(rgbi);
    TS_ASSERT_DELTA( gi.v(), 110, 1 );

    PixelRGBA<float> rgbaf(0.8,0.4,0.7,1.0);
    PixelGrayA<float> gaf = weighted_rgb_to_gray(rgbaf);
    TS_ASSERT_DELTA( gaf.v(), 0.5530, 1e-4 );

    PixelRGBA<uint8> rgbai(180,56,212,255);
    PixelGrayA<uint8> gai = weighted_rgb_to_gray(rgbai);
    TS_ASSERT_DELTA( gai.v(), 110, 1 );
  }

  void test_pixel_hsv()
  {
    // Test default-construction and size with all supported channel types
    { PixelHSV<int8> p; TS_ASSERT( p.h()==0 && p.s()==0 && p.v()==0 ); TS_ASSERT( sizeof(p)==3 ); }
    { PixelHSV<uint8> p; TS_ASSERT( p.h()==0 && p.s()==0 && p.v()==0 ); TS_ASSERT( sizeof(p)==3 ); }
    { PixelHSV<int16> p; TS_ASSERT( p.h()==0 && p.s()==0 && p.v()==0 ); TS_ASSERT( sizeof(p)==6 ); }
    { PixelHSV<uint16> p; TS_ASSERT( p.h()==0 && p.s()==0 && p.v()==0 ); TS_ASSERT( sizeof(p)==6 ); }
    { PixelHSV<int32> p; TS_ASSERT( p.h()==0 && p.s()==0 && p.v()==0 ); TS_ASSERT( sizeof(p)==12 ); }
    { PixelHSV<uint32> p; TS_ASSERT( p.h()==0 && p.s()==0 && p.v()==0 ); TS_ASSERT( sizeof(p)==12 ); }
    { PixelHSV<int64> p; TS_ASSERT( p.h()==0 && p.s()==0 && p.v()==0 ); TS_ASSERT( sizeof(p)==24 ); }
    { PixelHSV<uint64> p; TS_ASSERT( p.h()==0 && p.s()==0 && p.v()==0 ); TS_ASSERT( sizeof(p)==24 ); }
    { PixelHSV<float32> p; TS_ASSERT( p.h()==0.0 && p.s()==0.0 && p.v()==0.0 ); TS_ASSERT( sizeof(p)==12 ); }
    { PixelHSV<float64> p; TS_ASSERT( p.h()==0.0 && p.s()==0.0 && p.v()==0.0 ); TS_ASSERT( sizeof(p)==24 ); }
    // Test channel-value-construction and accessors
    { PixelHSV<int8> p(1,2,3); TS_ASSERT( p.h()==1 && p.s()==2 && p.v()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelHSV<uint8> p(1,2,3); TS_ASSERT( p.h()==1 && p.s()==2 && p.v()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelHSV<int16> p(1,2,3); TS_ASSERT( p.h()==1 && p.s()==2 && p.v()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelHSV<uint16> p(1,2,3); TS_ASSERT( p.h()==1 && p.s()==2 && p.v()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelHSV<int32> p(1,2,3); TS_ASSERT( p.h()==1 && p.s()==2 && p.v()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelHSV<uint32> p(1,2,3); TS_ASSERT( p.h()==1 && p.s()==2 && p.v()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelHSV<int64> p(1,2,3); TS_ASSERT( p.h()==1 && p.s()==2 && p.v()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelHSV<uint64> p(1,2,3); TS_ASSERT( p.h()==1 && p.s()==2 && p.v()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelHSV<float32> p(1.0,2.0,3.0); TS_ASSERT( p.h()==1.0 && p.s()==2.0 && p.v()==3.0 ); TS_ASSERT( p[0]==1.0 && p[1]==2.0 && p[2]==3.0 ); TS_ASSERT( p(0)==1.0 && p(1)==2.0 && p(2)==3.0 ); }
    { PixelHSV<float64> p(1.0,2.0,3.0); TS_ASSERT( p.h()==1.0 && p.s()==2.0 && p.v()==3.0 ); TS_ASSERT( p[0]==1.0 && p[1]==2.0 && p[2]==3.0 ); TS_ASSERT( p(0)==1.0 && p(1)==2.0 && p(2)==3.0 ); }
  }

  void test_rgb_to_hsv()
  {
    PixelRGB<float> input_rgb(1, 1, 1);
    PixelHSV<float> test_hsv(input_rgb);
    TS_ASSERT_EQUALS(test_hsv.h(), 0.0);
    TS_ASSERT_EQUALS(test_hsv.s(), 0.0);
    TS_ASSERT_EQUALS(test_hsv.v(), 1.0);
    
    PixelRGB<uint8> input_rgb8(100, 100, 100);
    PixelHSV<uint8> test_hsv8(input_rgb8);
    TS_ASSERT_EQUALS(test_hsv8.h(), 0);
    TS_ASSERT_EQUALS(test_hsv8.s(), 0);
    TS_ASSERT_EQUALS(test_hsv8.v(), 100);

    PixelRGB<uint16> input_rgb16(100, 100, 100);
    PixelHSV<uint16> test_hsv16(input_rgb16);
    TS_ASSERT_EQUALS(test_hsv16.h(), 0);
    TS_ASSERT_EQUALS(test_hsv16.s(), 0);
    TS_ASSERT_EQUALS(test_hsv16.v(), 100);
  }

  void test_hsv_to_rgb()
  {
    PixelHSV<float> input_hsv(0, 0, 1);
    PixelRGB<float> test_rgb(input_hsv);
    TS_ASSERT_EQUALS(test_rgb.r(), 1);
    TS_ASSERT_EQUALS(test_rgb.g(), 1);
    TS_ASSERT_EQUALS(test_rgb.b(), 1);

    PixelHSV<float> input_hsv_wrap_h(1, 0, 1);
    PixelRGB<float> test_rgb_wrap(input_hsv_wrap_h);
    TS_ASSERT_EQUALS(test_rgb_wrap.r(), 1);
    TS_ASSERT_EQUALS(test_rgb_wrap.g(), 1);
    TS_ASSERT_EQUALS(test_rgb_wrap.b(), 1);

    PixelHSV<uint8> input_hsv8(0, 0, 100);
    PixelRGB<uint8> test_rgb8(input_hsv8);
    TS_ASSERT_EQUALS(test_rgb8.r(), 100);
    TS_ASSERT_EQUALS(test_rgb8.g(), 100);
    TS_ASSERT_EQUALS(test_rgb8.b(), 100);

    PixelHSV<uint16> input_hsv16(0, 0, 100);
    PixelRGB<uint16> test_rgb16(input_hsv16);
    TS_ASSERT_EQUALS(test_rgb16.r(), 100);
    TS_ASSERT_EQUALS(test_rgb16.g(), 100);
    TS_ASSERT_EQUALS(test_rgb16.b(), 100);
  }

  void test_hsv_to_rgb_to_hsv()
  {
    for( double h=0; h<1; h+=0.15 ) {
      for( double s=0.2; s<=1; s+=0.2 ) {
        for( double v=0.2; v<=1; v+=0.2 ) {
	  PixelHSV<double> tmp1(h,s,v);
	  PixelRGB<double> tmp2(tmp1);
          PixelHSV<double> hsv(tmp2);
          TS_ASSERT_DELTA( hsv.h(), h, 1e-4 );
          TS_ASSERT_DELTA( hsv.s(), s, 1e-4 );
          TS_ASSERT_DELTA( hsv.v(), v, 1e-4 );
        }
      }
    }
    // Rounding error can become significant for small numbers, so we  
    // restrict this test to sufficiently bright and saturated pixels.
    for( uint8 h=0; ; h+=5 ) {
      for( uint8 s=60; ; s+=5 ) {
        for( uint8 v=80; ; v+=5 ) {
	  PixelHSV<uint8> tmp1(h,s,v);
	  PixelRGB<uint8> tmp2(tmp1);
          PixelHSV<uint8> hsv(tmp2);
          TS_ASSERT( abs(hsv.h()-h)<=2 || abs(abs(hsv.h()-h)-256)<=2 );
          TS_ASSERT_DELTA( hsv.s(), s, 2 );
          TS_ASSERT_EQUALS( hsv.v(), v );
          if( v == 255 ) break;
        }
        if( s == 255 ) break;
      }
      if( h == 255 ) break;
    }
  }

  void test_rgb_to_hsv_to_rgb()
  {
    for( double r=0; r<=1; r+=0.1 ) {
      for( double g=0; g<=1; g+=0.1 ) {
        for( double b=0; b<=1; b+=0.1 ) {
	  PixelRGB<double> tmp1(r,g,b);
	  PixelHSV<double> tmp2(tmp1);
          PixelRGB<double> rgb(tmp2);
          TS_ASSERT_DELTA( rgb.r(), r, 1e-4 );
          TS_ASSERT_DELTA( rgb.g(), g, 1e-4 );
          TS_ASSERT_DELTA( rgb.b(), b, 1e-4 );
        }
      }
    }
    // Rounding error can become significant for small numbers, so 
    // we restrict this test to sufficiently bright pixels.
    for( uint8 r=70; ; r+=5 ) {
      for( uint8 g=70; ; g+=5 ) {
        for( uint8 b=70; ; b+=5 ) {
	  PixelRGB<uint8> tmp1(r,g,b);
	  PixelHSV<uint8> tmp2(tmp1);
          PixelRGB<uint8> rgb(tmp2);
          TS_ASSERT_DELTA( rgb.r(), r, 2 );
          TS_ASSERT_DELTA( rgb.g(), g, 2 );
          TS_ASSERT_DELTA( rgb.b(), b, 2 );
          if( b == 255 ) break;
        }
        if( g == 255 ) break;
      }
      if( r == 255 ) break;
    }
  }

  void test_pixel_xyz()
  {
    // Test default-construction and size with all supported channel types
    { PixelXYZ<int8> p; TS_ASSERT( p.x()==0 && p.y()==0 && p.z()==0 ); TS_ASSERT( sizeof(p)==3 ); }
    { PixelXYZ<uint8> p; TS_ASSERT( p.x()==0 && p.y()==0 && p.z()==0 ); TS_ASSERT( sizeof(p)==3 ); }
    { PixelXYZ<int16> p; TS_ASSERT( p.x()==0 && p.y()==0 && p.z()==0 ); TS_ASSERT( sizeof(p)==6 ); }
    { PixelXYZ<uint16> p; TS_ASSERT( p.x()==0 && p.y()==0 && p.z()==0 ); TS_ASSERT( sizeof(p)==6 ); }
    { PixelXYZ<int32> p; TS_ASSERT( p.x()==0 && p.y()==0 && p.z()==0 ); TS_ASSERT( sizeof(p)==12 ); }
    { PixelXYZ<uint32> p; TS_ASSERT( p.x()==0 && p.y()==0 && p.z()==0 ); TS_ASSERT( sizeof(p)==12 ); }
    { PixelXYZ<int64> p; TS_ASSERT( p.x()==0 && p.y()==0 && p.z()==0 ); TS_ASSERT( sizeof(p)==24 ); }
    { PixelXYZ<uint64> p; TS_ASSERT( p.x()==0 && p.y()==0 && p.z()==0 ); TS_ASSERT( sizeof(p)==24 ); }
    { PixelXYZ<float32> p; TS_ASSERT( p.x()==0.0 && p.y()==0.0 && p.z()==0.0 ); TS_ASSERT( sizeof(p)==12 ); }
    { PixelXYZ<float64> p; TS_ASSERT( p.x()==0.0 && p.y()==0.0 && p.z()==0.0 ); TS_ASSERT( sizeof(p)==24 ); }
    // Test channel-value-construction and accessors
    { PixelXYZ<int8> p(1,2,3); TS_ASSERT( p.x()==1 && p.y()==2 && p.z()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelXYZ<uint8> p(1,2,3); TS_ASSERT( p.x()==1 && p.y()==2 && p.z()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelXYZ<int16> p(1,2,3); TS_ASSERT( p.x()==1 && p.y()==2 && p.z()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelXYZ<uint16> p(1,2,3); TS_ASSERT( p.x()==1 && p.y()==2 && p.z()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelXYZ<int32> p(1,2,3); TS_ASSERT( p.x()==1 && p.y()==2 && p.z()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelXYZ<uint32> p(1,2,3); TS_ASSERT( p.x()==1 && p.y()==2 && p.z()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelXYZ<int64> p(1,2,3); TS_ASSERT( p.x()==1 && p.y()==2 && p.z()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelXYZ<uint64> p(1,2,3); TS_ASSERT( p.x()==1 && p.y()==2 && p.z()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelXYZ<float32> p(1.0,2.0,3.0); TS_ASSERT( p.x()==1.0 && p.y()==2.0 && p.z()==3.0 ); TS_ASSERT( p[0]==1.0 && p[1]==2.0 && p[2]==3.0 ); TS_ASSERT( p(0)==1.0 && p(1)==2.0 && p(2)==3.0 ); }
    { PixelXYZ<float64> p(1.0,2.0,3.0); TS_ASSERT( p.x()==1.0 && p.y()==2.0 && p.z()==3.0 ); TS_ASSERT( p[0]==1.0 && p[1]==2.0 && p[2]==3.0 ); TS_ASSERT( p(0)==1.0 && p(1)==2.0 && p(2)==3.0 ); }
  }

  void test_pixel_rgb_to_xyz_to_rgb()
  {
    for( double r=0; r<=1; r+=0.1 ) {
      for( double g=0; g<=1; g+=0.1 ) {
        for( double b=0; b<=1; b+=0.1 ) {
	  PixelRGB<double> tmp1(r,g,b);
	  PixelXYZ<double> tmp2(tmp1);
          PixelRGB<double> rgb(tmp2);
          TS_ASSERT_DELTA( rgb.r(), r, 1e-4 );
          TS_ASSERT_DELTA( rgb.g(), g, 1e-4 );
          TS_ASSERT_DELTA( rgb.b(), b, 1e-4 );
        }
      }
    }
    // Values near the top of the range can clamp in XYZ 
    // space, so we stop at 235.
    for( uint8 r=0; r<=235; r+=5 ) {
      for( uint8 g=0; g<=235; g+=5 ) {
        for( uint8 b=0; b<=235; b+=5 ) {
	  PixelRGB<uint8> tmp1(r,g,b);
	  PixelXYZ<uint8> tmp2(tmp1);
          PixelRGB<uint8> rgb(tmp2);
          TS_ASSERT_DELTA( rgb.r(), r, 3 );
          TS_ASSERT_DELTA( rgb.g(), g, 3 );
          TS_ASSERT_DELTA( rgb.b(), b, 3 );
        }
      }
    }
  }

};


