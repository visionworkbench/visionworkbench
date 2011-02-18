// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestPixelTypes.h
#include <gtest/gtest.h>

#include <vw/Image/PixelTypes.h>
#include <vw/Image/PixelMask.h>

#include <test/Helpers.h>

using namespace vw;

TEST( PixelTypes, ChannelTypes ) {
  EXPECT_EQ( ChannelRange<vw::uint8>::max(), 255 );
  EXPECT_EQ( ChannelRange<vw::uint16>::max(), 65535 );
  EXPECT_EQ( ChannelRange<vw::float32>::max(), 1.0 );
  EXPECT_EQ( ChannelRange<vw::float64>::max(), 1.0 );

  EXPECT_EQ( channel_cast<vw::uint16>(vw::uint8(255)), 255 );
  EXPECT_EQ( channel_cast_rescale<vw::uint16>(vw::uint8(255)), 65535 );
  EXPECT_EQ( channel_cast<vw::uint8>(vw::float32(17.0)), 17 );
  EXPECT_EQ( channel_cast_rescale<vw::uint8>(vw::float32(0.333334)), 85 );
}

TEST( PixelTypes, Gray ) {
  // Test default-construction and size with all supported channel types
  { PixelGray<vw::int8> p; EXPECT_EQ( p.v(),0 ); EXPECT_EQ( sizeof(p),1u ); }
  { PixelGray<vw::uint8> p; EXPECT_EQ( p.v(),0 ); EXPECT_EQ( sizeof(p),1u ); }
  { PixelGray<vw::int16> p; EXPECT_EQ( p.v(),0 ); EXPECT_EQ( sizeof(p),2u ); }
  { PixelGray<vw::uint16> p; EXPECT_EQ( p.v(),0u ); EXPECT_EQ( sizeof(p),2u ); }
  { PixelGray<vw::int32> p; EXPECT_EQ( p.v(),0 ); EXPECT_EQ( sizeof(p),4u ); }
  { PixelGray<vw::uint32> p; EXPECT_EQ( p.v(),0u ); EXPECT_EQ( sizeof(p),4u ); }
  { PixelGray<vw::int64> p; EXPECT_EQ( p.v(),0 ); EXPECT_EQ( sizeof(p),8u ); }
  { PixelGray<vw::uint64> p; EXPECT_EQ( p.v(),0u ); EXPECT_EQ( sizeof(p),8u ); }
  { PixelGray<vw::float32> p; EXPECT_EQ( p.v(),0.0 ); EXPECT_EQ( sizeof(p),4u ); }
  { PixelGray<vw::float64> p; EXPECT_EQ( p.v(),0.0 ); EXPECT_EQ( sizeof(p),8u ); }
  // Test channel-value-construction and accessors
  { PixelGray<vw::int8> p(1); EXPECT_EQ( p.v(),1 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p(0),1 ); }
  { PixelGray<vw::uint8> p(1); EXPECT_EQ( p.v(),1 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p(0),1 ); }
  { PixelGray<vw::int16> p(1); EXPECT_EQ( p.v(),1 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p(0),1 ); }
  { PixelGray<vw::uint16> p(1); EXPECT_EQ( p.v(),1u ); EXPECT_EQ( p[0],1u ); EXPECT_EQ( p(0),1u ); }
  { PixelGray<vw::int32> p(1); EXPECT_EQ( p.v(),1 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p(0),1 ); }
  { PixelGray<vw::uint32> p(1); EXPECT_EQ( p.v(),1u ); EXPECT_EQ( p[0],1u ); EXPECT_EQ( p(0),1u ); }
  { PixelGray<vw::int64> p(1); EXPECT_EQ( p.v(),1 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p(0),1 ); }
  { PixelGray<vw::uint64> p(1); EXPECT_EQ( p.v(),1u ); EXPECT_EQ( p[0],1u ); EXPECT_EQ( p(0),1u ); }
  { PixelGray<vw::float32> p(1.0); EXPECT_EQ( p.v(),1.0 ); EXPECT_EQ( p[0],1.0 ); EXPECT_EQ( p(0),1.0 ); }
  { PixelGray<vw::float64> p(1.0); EXPECT_EQ( p.v(),1.0 ); EXPECT_EQ( p[0],1.0 ); EXPECT_EQ( p(0),1.0 ); }
}

TEST( PixelTypes, GrayA ) {
  // Test default-construction and size with all supported channel types
  { PixelGrayA<vw::int8> p; EXPECT_EQ( p.v(),0 ); EXPECT_EQ( p.a(),0 ); EXPECT_EQ( sizeof(p),2u ); }
  { PixelGrayA<vw::uint8> p; EXPECT_EQ( p.v(),0 ); EXPECT_EQ( p.a(),0 ); EXPECT_EQ( sizeof(p),2u ); }
  { PixelGrayA<vw::int16> p; EXPECT_EQ( p.v(),0 ); EXPECT_EQ( p.a(),0 ); EXPECT_EQ( sizeof(p),4u ); }
  { PixelGrayA<vw::uint16> p; EXPECT_EQ( p.v(),0u ); EXPECT_EQ( p.a(),0u ); EXPECT_EQ( sizeof(p),4u ); }
  { PixelGrayA<vw::int32> p; EXPECT_EQ( p.v(),0 ); EXPECT_EQ( p.a(),0 ); EXPECT_EQ( sizeof(p),8u ); }
  { PixelGrayA<vw::uint32> p; EXPECT_EQ( p.v(),0u ); EXPECT_EQ( p.a(),0u ); EXPECT_EQ( sizeof(p),8u ); }
  { PixelGrayA<vw::int64> p; EXPECT_EQ( p.v(),0 ); EXPECT_EQ( p.a(),0 ); EXPECT_EQ( sizeof(p),16u ); }
  { PixelGrayA<vw::uint64> p; EXPECT_EQ( p.v(),0u ); EXPECT_EQ( p.a(),0u ); EXPECT_EQ( sizeof(p),16u ); }
  { PixelGrayA<vw::float32> p; EXPECT_EQ( p.v(),0.0 ); EXPECT_EQ( p.a(),0.0 ); EXPECT_EQ( sizeof(p),8u ); }
  { PixelGrayA<vw::float64> p; EXPECT_EQ( p.v(),0.0 ); EXPECT_EQ( p.a(),0.0 ); EXPECT_EQ( sizeof(p),16u ); }
  // Test channel-value-construction and accessors
  { PixelGrayA<vw::int8> p(1,2); EXPECT_EQ( p.v(),1 ); EXPECT_EQ( p.a(),2 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); }
  { PixelGrayA<vw::uint8> p(1,2); EXPECT_EQ( p.v(),1 ); EXPECT_EQ( p.a(),2 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); }
  { PixelGrayA<vw::int16> p(1,2); EXPECT_EQ( p.v(),1 ); EXPECT_EQ( p.a(),2 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); }
  { PixelGrayA<vw::uint16> p(1,2); EXPECT_EQ( p.v(),1 ); EXPECT_EQ( p.a(),2 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); }
  { PixelGrayA<vw::int32> p(1,2); EXPECT_EQ( p.v(),1 ); EXPECT_EQ( p.a(),2 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); }
  { PixelGrayA<vw::uint32> p(1,2); EXPECT_EQ( p.v(),1u ); EXPECT_EQ( p.a(),2u ); EXPECT_EQ( p[0],1u ); EXPECT_EQ( p[1],2u ); EXPECT_EQ( p(0),1u ); EXPECT_EQ( p(1),2u ); }
  { PixelGrayA<vw::int64> p(1,2); EXPECT_EQ( p.v(),1 ); EXPECT_EQ( p.a(),2 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); }
  { PixelGrayA<vw::uint64> p(1,2); EXPECT_EQ( p.v(),1u ); EXPECT_EQ( p.a(),2u ); EXPECT_EQ( p[0],1u ); EXPECT_EQ( p[1],2u ); EXPECT_EQ( p(0),1u ); EXPECT_EQ( p(1),2u ); }
  { PixelGrayA<vw::float32> p(1.0,2.0); EXPECT_EQ( p.v(),1.0 ); EXPECT_EQ( p.a(),2.0 ); EXPECT_EQ( p[0],1.0 ); EXPECT_EQ( p[1],2.0 ); EXPECT_EQ( p(0),1.0 ); EXPECT_EQ( p(1),2.0 ); }
  { PixelGrayA<vw::float64> p(1.0,2.0); EXPECT_EQ( p.v(),1.0 ); EXPECT_EQ( p.a(),2.0 ); EXPECT_EQ( p[0],1.0 ); EXPECT_EQ( p[1],2.0 ); EXPECT_EQ( p(0),1.0 ); EXPECT_EQ( p(1),2.0 ); }
}

TEST( PixelTypes, RGB ) {
  // Test default-construction and size with all supported channel types
  { PixelRGB<vw::int8> p; EXPECT_EQ( p.r(),0 ); EXPECT_EQ( p.g(),0 ); EXPECT_EQ( p.b(),0 ); EXPECT_EQ( sizeof(p),3u ); }
  { PixelRGB<vw::uint8> p; EXPECT_EQ( p.r(),0 ); EXPECT_EQ( p.g(),0 ); EXPECT_EQ( p.b(),0 ); EXPECT_EQ( sizeof(p),3u ); }
  { PixelRGB<vw::int16> p; EXPECT_EQ( p.r(),0 ); EXPECT_EQ( p.g(),0 ); EXPECT_EQ( p.b(),0 ); EXPECT_EQ( sizeof(p),6u ); }
  { PixelRGB<vw::uint16> p; EXPECT_EQ( p.r(),0u ); EXPECT_EQ( p.g(),0u ); EXPECT_EQ( p.b(),0u ); EXPECT_EQ( sizeof(p),6u ); }
  { PixelRGB<vw::int32> p; EXPECT_EQ( p.r(),0 ); EXPECT_EQ( p.g(),0 ); EXPECT_EQ( p.b(),0 ); EXPECT_EQ( sizeof(p),12u ); }
  { PixelRGB<vw::uint32> p; EXPECT_EQ( p.r(),0u ); EXPECT_EQ( p.g(),0u ); EXPECT_EQ( p.b(),0u ); EXPECT_EQ( sizeof(p),12u ); }
  { PixelRGB<vw::int64> p; EXPECT_EQ( p.r(),0 ); EXPECT_EQ( p.g(),0 ); EXPECT_EQ( p.b(),0 ); EXPECT_EQ( sizeof(p),24u ); }
  { PixelRGB<vw::uint64> p; EXPECT_EQ( p.r(),0u ); EXPECT_EQ( p.g(),0u ); EXPECT_EQ( p.b(),0u ); EXPECT_EQ( sizeof(p),24u ); }
  { PixelRGB<vw::float32> p; EXPECT_EQ( p.r(),0.0 ); EXPECT_EQ( p.g(),0.0 ); EXPECT_EQ( p.b(),0.0 ); EXPECT_EQ( sizeof(p),12u ); }
  { PixelRGB<vw::float64> p; EXPECT_EQ( p.r(),0.0 ); EXPECT_EQ( p.g(),0.0 ); EXPECT_EQ( p.b(),0.0 ); EXPECT_EQ( sizeof(p),24u ); }
  // Test channel-value-construction and accessors
  { PixelRGB<vw::int8> p(1,2,3); EXPECT_EQ( p.r(),1 ); EXPECT_EQ( p.g(),2 ); EXPECT_EQ( p.b(),3 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); }
  { PixelRGB<vw::uint8> p(1,2,3); EXPECT_EQ( p.r(),1 ); EXPECT_EQ( p.g(),2 ); EXPECT_EQ( p.b(),3 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); }
  { PixelRGB<vw::int16> p(1,2,3); EXPECT_EQ( p.r(),1 ); EXPECT_EQ( p.g(),2 ); EXPECT_EQ( p.b(),3 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); }
  { PixelRGB<vw::uint16> p(1,2,3); EXPECT_EQ( p.r(),1 ); EXPECT_EQ( p.g(),2 ); EXPECT_EQ( p.b(),3 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); }
  { PixelRGB<vw::int32> p(1,2,3); EXPECT_EQ( p.r(),1 ); EXPECT_EQ( p.g(),2 ); EXPECT_EQ( p.b(),3 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); }
  { PixelRGB<vw::uint32> p(1,2,3); EXPECT_EQ( p.r(),1u ); EXPECT_EQ( p.g(),2u ); EXPECT_EQ( p.b(),3u ); EXPECT_EQ( p[0],1u ); EXPECT_EQ( p[1],2u ); EXPECT_EQ( p[2],3u ); EXPECT_EQ( p(0),1u ); EXPECT_EQ( p(1),2u ); EXPECT_EQ( p(2),3u ); }
  { PixelRGB<vw::int64> p(1,2,3); EXPECT_EQ( p.r(),1 ); EXPECT_EQ( p.g(),2 ); EXPECT_EQ( p.b(),3 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); }
  { PixelRGB<vw::uint64> p(1,2,3); EXPECT_EQ( p.r(),1u ); EXPECT_EQ( p.g(),2u ); EXPECT_EQ( p.b(),3u ); EXPECT_EQ( p[0],1u ); EXPECT_EQ( p[1],2u ); EXPECT_EQ( p[2],3u ); EXPECT_EQ( p(0),1u ); EXPECT_EQ( p(1),2u ); EXPECT_EQ( p(2),3u ); }
  { PixelRGB<vw::float32> p(1.0,2.0,3.0); EXPECT_EQ( p.r(),1.0 ); EXPECT_EQ( p.g(),2.0 ); EXPECT_EQ( p.b(),3.0 ); EXPECT_EQ( p[0],1.0 ); EXPECT_EQ( p[1],2.0 ); EXPECT_EQ( p[2],3.0 ); EXPECT_EQ( p(0),1.0 ); EXPECT_EQ( p(1),2.0 ); EXPECT_EQ( p(2),3.0 ); }
  { PixelRGB<vw::float64> p(1.0,2.0,3.0); EXPECT_EQ( p.r(),1.0 ); EXPECT_EQ( p.g(),2.0 ); EXPECT_EQ( p.b(),3.0 ); EXPECT_EQ( p[0],1.0 ); EXPECT_EQ( p[1],2.0 ); EXPECT_EQ( p[2],3.0 ); EXPECT_EQ( p(0),1.0 ); EXPECT_EQ( p(1),2.0 ); EXPECT_EQ( p(2),3.0 ); }
}

TEST( PixelTypes, RGBA ) {
  // Test default-construction and size with all supported channel types
  { PixelRGBA<vw::int8> p; EXPECT_EQ( p.r(),0 ); EXPECT_EQ( p.g(),0 ); EXPECT_EQ( p.b(),0 ); EXPECT_EQ( p.a(),0 ); EXPECT_EQ( sizeof(p),4u ); }
  { PixelRGBA<vw::uint8> p; EXPECT_EQ( p.r(),0 ); EXPECT_EQ( p.g(),0 ); EXPECT_EQ( p.b(),0 ); EXPECT_EQ( p.a(),0 ); EXPECT_EQ( sizeof(p),4u ); }
  { PixelRGBA<vw::int16> p; EXPECT_EQ( p.r(),0 ); EXPECT_EQ( p.g(),0 ); EXPECT_EQ( p.b(),0 ); EXPECT_EQ( p.a(),0 ); EXPECT_EQ( sizeof(p),8u ); }
  { PixelRGBA<vw::uint16> p; EXPECT_EQ( p.r(),0 ); EXPECT_EQ( p.g(),0 ); EXPECT_EQ( p.b(),0 ); EXPECT_EQ( p.a(),0 ); EXPECT_EQ( sizeof(p),8u ); }
  { PixelRGBA<vw::int32> p; EXPECT_EQ( p.r(),0 ); EXPECT_EQ( p.g(),0 ); EXPECT_EQ( p.b(),0 ); EXPECT_EQ( p.a(),0 ); EXPECT_EQ( sizeof(p),16u ); }
  { PixelRGBA<vw::uint32> p; EXPECT_EQ( p.r(),0u ); EXPECT_EQ( p.g(),0u ); EXPECT_EQ( p.b(),0u ); EXPECT_EQ( p.a(),0u ); EXPECT_EQ( sizeof(p),16u ); }
  { PixelRGBA<vw::int64> p; EXPECT_EQ( p.r(),0 ); EXPECT_EQ( p.g(),0 ); EXPECT_EQ( p.b(),0 ); EXPECT_EQ( p.a(),0 ); EXPECT_EQ( sizeof(p),32u ); }
  { PixelRGBA<vw::uint64> p; EXPECT_EQ( p.r(),0u ); EXPECT_EQ( p.g(),0u ); EXPECT_EQ( p.b(),0u ); EXPECT_EQ( p.a(),0u ); EXPECT_EQ( sizeof(p),32u ); }
  { PixelRGBA<vw::float32> p; EXPECT_EQ( p.r(),0.0 ); EXPECT_EQ( p.g(),0.0 ); EXPECT_EQ( p.b(),0.0 ); EXPECT_EQ( p.a(),0.0 ); EXPECT_EQ( sizeof(p),16u ); }
  { PixelRGBA<vw::float64> p; EXPECT_EQ( p.r(),0.0 ); EXPECT_EQ( p.g(),0.0 ); EXPECT_EQ( p.b(),0.0 ); EXPECT_EQ( p.a(),0.0 ); EXPECT_EQ( sizeof(p),32u ); }
  // Test channel-value-construction and accessors
  { PixelRGBA<vw::int8> p(1,2,3,4); EXPECT_EQ( p.r(),1 ); EXPECT_EQ( p.g(),2 ); EXPECT_EQ( p.b(),3 ); EXPECT_EQ( p.a(),4 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p[3],4 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); EXPECT_EQ( p(3),4 ); }
  { PixelRGBA<vw::uint8> p(1,2,3,4); EXPECT_EQ( p.r(),1 ); EXPECT_EQ( p.g(),2 ); EXPECT_EQ( p.b(),3 ); EXPECT_EQ( p.a(),4 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p[3],4 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); EXPECT_EQ( p(3),4 ); }
  { PixelRGBA<vw::int16> p(1,2,3,4); EXPECT_EQ( p.r(),1 ); EXPECT_EQ( p.g(),2 ); EXPECT_EQ( p.b(),3 ); EXPECT_EQ( p.a(),4 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p[3],4 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); EXPECT_EQ( p(3),4 ); }
  { PixelRGBA<vw::uint16> p(1,2,3,4); EXPECT_EQ( p.r(),1 ); EXPECT_EQ( p.g(),2 ); EXPECT_EQ( p.b(),3 ); EXPECT_EQ( p.a(),4 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p[3],4 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); EXPECT_EQ( p(3),4 ); }
  { PixelRGBA<vw::int32> p(1,2,3,4); EXPECT_EQ( p.r(),1 ); EXPECT_EQ( p.g(),2 ); EXPECT_EQ( p.b(),3 ); EXPECT_EQ( p.a(),4 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p[3],4 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); EXPECT_EQ( p(3),4 ); }
  { PixelRGBA<vw::uint32> p(1,2,3,4); EXPECT_EQ( p.r(),1u ); EXPECT_EQ( p.g(),2u ); EXPECT_EQ( p.b(),3u ); EXPECT_EQ( p.a(),4u ); EXPECT_EQ( p[0],1u ); EXPECT_EQ( p[1],2u ); EXPECT_EQ( p[2],3u ); EXPECT_EQ( p[3],4u ); EXPECT_EQ( p(0),1u ); EXPECT_EQ( p(1),2u ); EXPECT_EQ( p(2),3u ); EXPECT_EQ( p(3),4u ); }
  { PixelRGBA<vw::int64> p(1,2,3,4); EXPECT_EQ( p.r(),1 ); EXPECT_EQ( p.g(),2 ); EXPECT_EQ( p.b(),3 ); EXPECT_EQ( p.a(),4 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p[3],4 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); EXPECT_EQ( p(3),4 ); }
  { PixelRGBA<vw::uint64> p(1,2,3,4); EXPECT_EQ( p.r(),1u ); EXPECT_EQ( p.g(),2u ); EXPECT_EQ( p.b(),3u ); EXPECT_EQ( p.a(),4u ); EXPECT_EQ( p[0],1u ); EXPECT_EQ( p[1],2u ); EXPECT_EQ( p[2],3u ); EXPECT_EQ( p[3],4u ); EXPECT_EQ( p(0),1u ); EXPECT_EQ( p(1),2u ); EXPECT_EQ( p(2),3u ); EXPECT_EQ( p(3),4u ); }
  { PixelRGBA<vw::float32> p(1.0,2.0,3.0,4.0); EXPECT_EQ( p.r(),1.0 ); EXPECT_EQ( p.g(),2.0 ); EXPECT_EQ( p.b(),3.0 ); EXPECT_EQ( p.a(),4.0 ); EXPECT_EQ( p[0],1.0 ); EXPECT_EQ( p[1],2.0 ); EXPECT_EQ( p[2],3.0 ); EXPECT_EQ( p[3],4.0 ); EXPECT_EQ( p(0),1.0 ); EXPECT_EQ( p(1),2.0 ); EXPECT_EQ( p(2),3.0 ); EXPECT_EQ( p(3),4.0 ); }
  { PixelRGBA<vw::float64> p(1.0,2.0,3.0,4.0); EXPECT_EQ( p.r(),1.0 ); EXPECT_EQ( p.g(),2.0 ); EXPECT_EQ( p.b(),3.0 ); EXPECT_EQ( p.a(),4.0 ); EXPECT_EQ( p[0],1.0 ); EXPECT_EQ( p[1],2.0 ); EXPECT_EQ( p[2],3.0 ); EXPECT_EQ( p[3],4.0 ); EXPECT_EQ( p(0),1.0 ); EXPECT_EQ( p(1),2.0 ); EXPECT_EQ( p(2),3.0 ); EXPECT_EQ( p(3),4.0 ); }
}

TEST( PixelTypes, RGB2Gray) {
  int val0, val1, val2;
  val0 = val1 = val2 = 40;

  // Standard case
  PixelRGB<vw::int8> test_rgb(val0, val1, val2);
  PixelGray<vw::int8> test_gray(test_rgb);
  EXPECT_EQ(test_gray.v(),
            (test_rgb.r() + test_rgb.g() + test_rgb.b())/3 );

  // Case with type incompatibility
  PixelRGB<vw::int16> test_rgb16(val0, val1, val2);
  PixelGray<vw::int8> test_gray8(test_rgb);
  EXPECT_EQ(test_gray8.v(),
            (test_rgb16.r() + test_rgb16.g() + test_rgb16.b())/3 );
}

TEST( PixelTypes, WeightedRGB2Gray) {
  PixelRGB<float> rgbf(0.8,0.4,0.7);
  PixelGray<float> gf = weighted_rgb_to_gray(rgbf);
  EXPECT_NEAR( gf.v(), 0.5530, 1e-4 );

  PixelRGB<vw::uint8> rgbi(180,56,212);
  PixelGray<vw::uint8> gi = weighted_rgb_to_gray(rgbi);
  EXPECT_NEAR( gi.v(), 110, 1 );

  PixelRGBA<float> rgbaf(0.8,0.4,0.7,1.0);
  PixelGrayA<float> gaf = weighted_rgb_to_gray(rgbaf);
  EXPECT_NEAR( gaf.v(), 0.5530, 1e-4 );

  PixelRGBA<vw::uint8> rgbai(180,56,212,255);
  PixelGrayA<vw::uint8> gai = weighted_rgb_to_gray(rgbai);
  EXPECT_NEAR( gai.v(), 110, 1 );
}

TEST( PixelTypes, HSV ) {
  // Test default-construction and size with all supported channel types
  { PixelHSV<vw::int8> p; EXPECT_EQ( p.h(),0 ); EXPECT_EQ( p.s(),0 ); EXPECT_EQ( p.v(),0 ); EXPECT_EQ( sizeof(p),3u ); }
  { PixelHSV<vw::uint8> p; EXPECT_EQ( p.h(),0 ); EXPECT_EQ( p.s(),0 ); EXPECT_EQ( p.v(),0 ); EXPECT_EQ( sizeof(p),3u ); }
  { PixelHSV<vw::int16> p; EXPECT_EQ( p.h(),0 ); EXPECT_EQ( p.s(),0 ); EXPECT_EQ( p.v(),0 ); EXPECT_EQ( sizeof(p),6u ); }
  { PixelHSV<vw::uint16> p; EXPECT_EQ( p.h(),0 ); EXPECT_EQ( p.s(),0 ); EXPECT_EQ( p.v(),0 ); EXPECT_EQ( sizeof(p),6u ); }
  { PixelHSV<vw::int32> p; EXPECT_EQ( p.h(),0 ); EXPECT_EQ( p.s(),0 ); EXPECT_EQ( p.v(),0 ); EXPECT_EQ( sizeof(p),12u ); }
  { PixelHSV<vw::uint32> p; EXPECT_EQ( p.h(),0u ); EXPECT_EQ( p.s(),0u ); EXPECT_EQ( p.v(),0u ); EXPECT_EQ( sizeof(p),12u ); }
  { PixelHSV<vw::int64> p; EXPECT_EQ( p.h(),0 ); EXPECT_EQ( p.s(),0 ); EXPECT_EQ( p.v(),0 ); EXPECT_EQ( sizeof(p),24u ); }
  { PixelHSV<vw::uint64> p; EXPECT_EQ( p.h(),0u ); EXPECT_EQ( p.s(),0u ); EXPECT_EQ( p.v(),0u ); EXPECT_EQ( sizeof(p),24u ); }
  { PixelHSV<vw::float32> p; EXPECT_EQ( p.h(),0.0 ); EXPECT_EQ( p.s(),0.0 ); EXPECT_EQ( p.v(),0.0 ); EXPECT_EQ( sizeof(p),12u ); }
  { PixelHSV<vw::float64> p; EXPECT_EQ( p.h(),0.0 ); EXPECT_EQ( p.s(),0.0 ); EXPECT_EQ( p.v(),0.0 ); EXPECT_EQ( sizeof(p),24u ); }
  // Test channel-value-construction and accessors
  { PixelHSV<vw::int8> p(1,2,3); EXPECT_EQ( p.h(),1 ); EXPECT_EQ( p.s(),2 ); EXPECT_EQ( p.v(),3 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); }
  { PixelHSV<vw::uint8> p(1,2,3); EXPECT_EQ( p.h(),1 ); EXPECT_EQ( p.s(),2 ); EXPECT_EQ( p.v(),3 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); }
  { PixelHSV<vw::int16> p(1,2,3); EXPECT_EQ( p.h(),1 ); EXPECT_EQ( p.s(),2 ); EXPECT_EQ( p.v(),3 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); }
  { PixelHSV<vw::uint16> p(1,2,3); EXPECT_EQ( p.h(),1 ); EXPECT_EQ( p.s(),2 ); EXPECT_EQ( p.v(),3 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); }
  { PixelHSV<vw::int32> p(1,2,3); EXPECT_EQ( p.h(),1 ); EXPECT_EQ( p.s(),2 ); EXPECT_EQ( p.v(),3 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); }
  { PixelHSV<vw::uint32> p(1,2,3); EXPECT_EQ( p.h(),1u ); EXPECT_EQ( p.s(),2u ); EXPECT_EQ( p.v(),3u ); EXPECT_EQ( p[0],1u ); EXPECT_EQ( p[1],2u ); EXPECT_EQ( p[2],3u ); EXPECT_EQ( p(0),1u ); EXPECT_EQ( p(1),2u ); EXPECT_EQ( p(2),3u ); }
  { PixelHSV<vw::int64> p(1,2,3); EXPECT_EQ( p.h(),1 ); EXPECT_EQ( p.s(),2 ); EXPECT_EQ( p.v(),3 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); }
  { PixelHSV<vw::uint64> p(1,2,3); EXPECT_EQ( p.h(),1u ); EXPECT_EQ( p.s(),2u ); EXPECT_EQ( p.v(),3u ); EXPECT_EQ( p[0],1u ); EXPECT_EQ( p[1],2u ); EXPECT_EQ( p[2],3u ); EXPECT_EQ( p(0),1u ); EXPECT_EQ( p(1),2u ); EXPECT_EQ( p(2),3u ); }
  { PixelHSV<vw::float32> p(1.0,2.0,3.0); EXPECT_EQ( p.h(),1.0 ); EXPECT_EQ( p.s(),2.0 ); EXPECT_EQ( p.v(),3.0 ); EXPECT_EQ( p[0],1.0 ); EXPECT_EQ( p[1],2.0 ); EXPECT_EQ( p[2],3.0 ); EXPECT_EQ( p(0),1.0 ); EXPECT_EQ( p(1),2.0 ); EXPECT_EQ( p(2),3.0 ); }
  { PixelHSV<vw::float64> p(1.0,2.0,3.0); EXPECT_EQ( p.h(),1.0 ); EXPECT_EQ( p.s(),2.0 ); EXPECT_EQ( p.v(),3.0 ); EXPECT_EQ( p[0],1.0 ); EXPECT_EQ( p[1],2.0 ); EXPECT_EQ( p[2],3.0 ); EXPECT_EQ( p(0),1.0 ); EXPECT_EQ( p(1),2.0 ); EXPECT_EQ( p(2),3.0 ); }
}

TEST( PixelTypes, RGB2HSV ) {
  PixelRGB<float> input_rgb(1, 1, 1);
  PixelHSV<float> test_hsv(input_rgb);
  EXPECT_PIXEL_EQ( test_hsv, PixelHSV<float>(0,0,1) );

  PixelRGB<vw::uint8> input_rgb8(100, 100, 100);
  PixelHSV<vw::uint8> test_hsv8(input_rgb8);
  EXPECT_PIXEL_EQ( test_hsv8, PixelHSV<uint8>(0,0,100) );

  PixelRGB<vw::uint16> input_rgb16(100, 100, 100);
  PixelHSV<vw::uint16> test_hsv16(input_rgb16);
  EXPECT_PIXEL_EQ( test_hsv16, PixelHSV<uint16>(0,0,100) );
}

TEST( PixelTypes, HSV2RGB ) {
  PixelHSV<float> input_hsv(0, 0, 1);
  PixelRGB<float> test_rgb(input_hsv);
  EXPECT_PIXEL_EQ( test_rgb, PixelRGB<float>(1,1,1) );

  PixelHSV<float> input_hsv_wrap_h(1, 0, 1);
  PixelRGB<float> test_rgb_wrap(input_hsv_wrap_h);
  EXPECT_PIXEL_EQ( test_rgb_wrap, PixelRGB<float>(1,1,1) );

  PixelHSV<vw::uint8> input_hsv8(0, 0, 100);
  PixelRGB<vw::uint8> test_rgb8(input_hsv8);
  EXPECT_PIXEL_EQ( test_rgb8, PixelRGB<vw::uint8>(100,100,100) );

  PixelHSV<vw::uint16> input_hsv16(0, 0, 100);
  PixelRGB<vw::uint16> test_rgb16(input_hsv16);
  EXPECT_PIXEL_EQ( test_rgb16, PixelRGB<vw::uint16>(100,100,100) );
}

TEST( PixelTypes, HSV2RGB2HSV ) {
  for( double h=0.05; h<1; h+=0.15 ) {
    for( double s=0.2; s<=1; s+=0.2 ) {
      for( double v=0.2; v<=1; v+=0.2 ) {
        PixelHSV<double> tmp1(h,s,v);
        PixelRGB<double> tmp2(tmp1);
        PixelHSV<double> hsv(tmp2);
        EXPECT_PIXEL_NEAR( hsv, tmp1, 1e-4 );
      }
    }
  }
  // Rounding error can become significant for small numbers, so we
  // restrict this test to sufficiently bright and saturated pixels.
  for( vw::uint8 h=0; ; h+=5 ) {
    for( vw::uint8 s=60; ; s+=5 ) {
      for( vw::uint8 v=80; ; v+=5 ) {
        PixelHSV<vw::uint8> tmp1(h,s,v);
        PixelRGB<vw::uint8> tmp2(tmp1);
        PixelHSV<vw::uint8> hsv(tmp2);
        EXPECT_TRUE( abs(hsv.h()-h)<=2 || abs(abs(hsv.h()-h)-256)<=2 );
        EXPECT_NEAR( hsv.s(), s, 2 );
        EXPECT_EQ( hsv.v(), v );
        if( v == 255 ) break;
      }
      if( s == 255 ) break;
    }
    if( h == 255 ) break;
  }
}

TEST( PixelTypes, RGB2HSV2RGB ) {
  for( double r=0; r<=1; r+=0.1 ) {
    for( double g=0; g<=1; g+=0.1 ) {
      for( double b=0; b<=1; b+=0.1 ) {
        PixelRGB<double> tmp1(r,g,b);
        PixelHSV<double> tmp2(tmp1);
        PixelRGB<double> rgb(tmp2);
        EXPECT_PIXEL_NEAR( rgb, tmp1, 1e-4 );
      }
    }
  }
  // Rounding error can become significant for small numbers, so
  // we restrict this test to sufficiently bright pixels.
  for( vw::uint8 r=70; ; r+=5 ) {
    for( vw::uint8 g=70; ; g+=5 ) {
      for( vw::uint8 b=70; ; b+=5 ) {
        PixelRGB<vw::uint8> tmp1(r,g,b);
        PixelHSV<vw::uint8> tmp2(tmp1);
        PixelRGB<vw::uint8> rgb(tmp2);
        EXPECT_PIXEL_NEAR( rgb, tmp1, 2 );
        if( b == 255 ) break;
      }
      if( g == 255 ) break;
    }
    if( r == 255 ) break;
  }
}

TEST( PixelTypes, XYZ ) {
  // Test default-construction and size with all supported channel types
  { PixelXYZ<vw::int8> p; EXPECT_EQ( p.x(),0 ); EXPECT_EQ( p.y(),0 ); EXPECT_EQ( p.z(),0 ); EXPECT_EQ( sizeof(p),3u ); }
  { PixelXYZ<vw::uint8> p; EXPECT_EQ( p.x(),0 ); EXPECT_EQ( p.y(),0 ); EXPECT_EQ( p.z(),0 ); EXPECT_EQ( sizeof(p),3u ); }
  { PixelXYZ<vw::int16> p; EXPECT_EQ( p.x(),0 ); EXPECT_EQ( p.y(),0 ); EXPECT_EQ( p.z(),0 ); EXPECT_EQ( sizeof(p),6u ); }
  { PixelXYZ<vw::uint16> p; EXPECT_EQ( p.x(),0 ); EXPECT_EQ( p.y(),0 ); EXPECT_EQ( p.z(),0 ); EXPECT_EQ( sizeof(p),6u ); }
  { PixelXYZ<vw::int32> p; EXPECT_EQ( p.x(),0 ); EXPECT_EQ( p.y(),0 ); EXPECT_EQ( p.z(),0 ); EXPECT_EQ( sizeof(p),12u ); }
  { PixelXYZ<vw::uint32> p; EXPECT_EQ( p.x(),0u ); EXPECT_EQ( p.y(),0u ); EXPECT_EQ( p.z(),0u ); EXPECT_EQ( sizeof(p),12u ); }
  { PixelXYZ<vw::int64> p; EXPECT_EQ( p.x(),0 ); EXPECT_EQ( p.y(),0 ); EXPECT_EQ( p.z(),0 ); EXPECT_EQ( sizeof(p),24u ); }
  { PixelXYZ<vw::uint64> p; EXPECT_EQ( p.x(),0u ); EXPECT_EQ( p.y(),0u ); EXPECT_EQ( p.z(),0u ); EXPECT_EQ( sizeof(p),24u ); }
  { PixelXYZ<vw::float32> p; EXPECT_EQ( p.x(),0.0 ); EXPECT_EQ( p.y(),0.0 ); EXPECT_EQ( p.z(),0.0 ); EXPECT_EQ( sizeof(p),12u ); }
  { PixelXYZ<vw::float64> p; EXPECT_EQ( p.x(),0.0 ); EXPECT_EQ( p.y(),0.0 ); EXPECT_EQ( p.z(),0.0 ); EXPECT_EQ( sizeof(p),24u ); }
  // Test channel-value-construction and accessors
  { PixelXYZ<vw::int8> p(1,2,3); EXPECT_EQ( p.x(),1 ); EXPECT_EQ( p.y(),2 ); EXPECT_EQ( p.z(),3 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); }
  { PixelXYZ<vw::uint8> p(1,2,3); EXPECT_EQ( p.x(),1 ); EXPECT_EQ( p.y(),2 ); EXPECT_EQ( p.z(),3 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); }
  { PixelXYZ<vw::int16> p(1,2,3); EXPECT_EQ( p.x(),1 ); EXPECT_EQ( p.y(),2 ); EXPECT_EQ( p.z(),3 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); }
  { PixelXYZ<vw::uint16> p(1,2,3); EXPECT_EQ( p.x(),1 ); EXPECT_EQ( p.y(),2 ); EXPECT_EQ( p.z(),3 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); }
  { PixelXYZ<vw::int32> p(1,2,3); EXPECT_EQ( p.x(),1 ); EXPECT_EQ( p.y(),2 ); EXPECT_EQ( p.z(),3 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); }
  { PixelXYZ<vw::uint32> p(1,2,3); EXPECT_EQ( p.x(),1u ); EXPECT_EQ( p.y(),2u ); EXPECT_EQ( p.z(),3u ); EXPECT_EQ( p[0],1u ); EXPECT_EQ( p[1],2u ); EXPECT_EQ( p[2],3u ); EXPECT_EQ( p(0),1u ); EXPECT_EQ( p(1),2u ); EXPECT_EQ( p(2),3u ); }
  { PixelXYZ<vw::int64> p(1,2,3); EXPECT_EQ( p.x(),1 ); EXPECT_EQ( p.y(),2 ); EXPECT_EQ( p.z(),3 ); EXPECT_EQ( p[0],1 ); EXPECT_EQ( p[1],2 ); EXPECT_EQ( p[2],3 ); EXPECT_EQ( p(0),1 ); EXPECT_EQ( p(1),2 ); EXPECT_EQ( p(2),3 ); }
  { PixelXYZ<vw::uint64> p(1,2,3); EXPECT_EQ( p.x(),1u ); EXPECT_EQ( p.y(),2u ); EXPECT_EQ( p.z(),3u ); EXPECT_EQ( p[0],1u ); EXPECT_EQ( p[1],2u ); EXPECT_EQ( p[2],3u ); EXPECT_EQ( p(0),1u ); EXPECT_EQ( p(1),2u ); EXPECT_EQ( p(2),3u ); }
  { PixelXYZ<vw::float32> p(1.0,2.0,3.0); EXPECT_EQ( p.x(),1.0 ); EXPECT_EQ( p.y(),2.0 ); EXPECT_EQ( p.z(),3.0 ); EXPECT_EQ( p[0],1.0 ); EXPECT_EQ( p[1],2.0 ); EXPECT_EQ( p[2],3.0 ); EXPECT_EQ( p(0),1.0 ); EXPECT_EQ( p(1),2.0 ); EXPECT_EQ( p(2),3.0 ); }
  { PixelXYZ<vw::float64> p(1.0,2.0,3.0); EXPECT_EQ( p.x(),1.0 ); EXPECT_EQ( p.y(),2.0 ); EXPECT_EQ( p.z(),3.0 ); EXPECT_EQ( p[0],1.0 ); EXPECT_EQ( p[1],2.0 ); EXPECT_EQ( p[2],3.0 ); EXPECT_EQ( p(0),1.0 ); EXPECT_EQ( p(1),2.0 ); EXPECT_EQ( p(2),3.0 ); }
}

TEST( PixelTypes, RGB2XYZ2RGB ) {
  for( double r=0; r<=1; r+=0.1 ) {
    for( double g=0; g<=1; g+=0.1 ) {
      for( double b=0; b<=1; b+=0.1 ) {
        PixelRGB<double> tmp1(r,g,b);
        PixelXYZ<double> tmp2(tmp1);
        PixelRGB<double> rgb(tmp2);
        EXPECT_PIXEL_NEAR( rgb, tmp1, 1e-4 );
      }
    }
  }
  // Values near the top of the range can clamp in XYZ
  // space, so we stop at 235.
  for( vw::uint8 r=0; r<=235; r+=5 ) {
    for( vw::uint8 g=0; g<=235; g+=5 ) {
      for( vw::uint8 b=0; b<=235; b+=5 ) {
        PixelRGB<vw::uint8> tmp1(r,g,b);
        PixelXYZ<vw::uint8> tmp2(tmp1);
        PixelRGB<vw::uint8> rgb(tmp2);
        EXPECT_PIXEL_NEAR( rgb, tmp1, 4 );
      }
    }
  }
}

TEST( PixelTypes, XYZ2RGB2XYZ ) {
  for( double x=0; x<=1; x+=0.1 ) {
    for( double y=0; y<=1; y+=0.1 ) {
      for( double z=0; z<=1; z+=0.1 ) {
        PixelXYZ<double> tmp1(x,y,z);
        PixelRGB<double> tmp2(tmp1);
        PixelXYZ<double> xyz(tmp2);
        EXPECT_PIXEL_NEAR( xyz, tmp1, 1e-4 );
      }
    }
  }
}

TEST( PixelTypes, XYZ2LUV2XYZ ) {
  // Omit solid black for Luv
  for( double x=0.1; x<=1; x+=0.1 ) {
    for( double y=0.1; y<=1; y+=0.1 ) {
      for( double z=0.1; z<=1; z+=0.1 ) {
        PixelXYZ<double> tmp1(x,y,z);
        PixelLuv<double> tmp2(tmp1);
        PixelXYZ<double> xyz(tmp2);
        EXPECT_PIXEL_NEAR( xyz, tmp1, 1e-4 );
      }
    }
  }
}

TEST( PixelTypes, LUV2XYZ2LUV ) {
  // Omit solid black for Luv
  for( double l=0.1; l<=1; l+=0.1 ) {
    for( double u=0.0; u<=1; u+=0.1 ) {
      for( double v=0.0; v<=1; v+=0.1 ) {
        PixelLuv<double> tmp1(l,u,v);
        PixelXYZ<double> tmp2(tmp1);
        PixelLuv<double> luv(tmp2);
        EXPECT_PIXEL_NEAR( luv, tmp1, 1e-4 );
      }
    }
  }
}

TEST( PixelTypes, RGB2LUV2RGB ) {
  // Omit solid black for Luv
  for( double r=0.1; r<=1; r+=0.1 ) {
    for( double g=0.1; g<=1; g+=0.1 ) {
      for( double b=0.1; b<=1; b+=0.1 ) {
        PixelRGB<double> tmp1(r,g,b);
        PixelLuv<double> tmp2(tmp1);
        PixelRGB<double> rgb(tmp2);
        EXPECT_PIXEL_NEAR( rgb, tmp1, 1e-4 );
      }
    }
  }
}

TEST( PixelTypes, LUV2RGB2LUV ) {
  // Omit solid black for Luv
  for( double l=0.1; l<=1; l+=0.1 ) {
    for( double u=0.0; u<=1; u+=0.1 ) {
      for( double v=0.0; v<=1; v+=0.1 ) {
        PixelLuv<double> tmp1(l,u,v);
        PixelRGB<double> tmp2(tmp1);
        PixelLuv<double> luv(tmp2);
        EXPECT_PIXEL_NEAR( luv, tmp1, 1e-4 );
      }
    }
  }
}

TEST( PixelTypes, XYZ2LAB2XYZ ) {
  // Omit solid black for Lab
  for( double x=0.1; x<=1; x+=0.1 ) {
    for( double y=0.1; y<=1; y+=0.1 ) {
      for( double z=0.1; z<=1; z+=0.1 ) {
        PixelXYZ<double> tmp1(x,y,z);
        PixelLab<double> tmp2(tmp1);
        PixelXYZ<double> xyz(tmp2);
        EXPECT_PIXEL_NEAR( xyz, tmp1, 1e-4 );
      }
    }
  }
}

TEST( PixelTypes, LAB2XYZ2LAB ) {
  // Omit solid black for Lab
  for( double l=0.1; l<=1; l+=0.1 ) {
    for( double a=0.0; a<=1; a+=0.1 ) {
      for( double b=0.0; b<=1; b+=0.1 ) {
        PixelLab<double> tmp1(l,a,b);
        PixelXYZ<double> tmp2(tmp1);
        PixelLab<double> lab(tmp2);
        EXPECT_PIXEL_NEAR( lab, tmp1, 1e-4 );
      }
    }
  }
}

TEST( PixelTypes, RGB2LAB2RGB ) {
  // Omit solid black for Lab
  for( double r=0.1; r<=1; r+=0.1 ) {
    for( double g=0.1; g<=1; g+=0.1 ) {
      for( double b=0.1; b<=1; b+=0.1 ) {
        PixelRGB<double> tmp1(r,g,b);
        PixelLab<double> tmp2(tmp1);
        PixelRGB<double> rgb(tmp2);
        EXPECT_PIXEL_NEAR( rgb, tmp1, 1e-4 );
      }
    }
  }
}

TEST( PixelTypes, LAB2RGB2LAB ) {
  // Omit solid black for Lab
  for( double l=0.1; l<=1; l+=0.1 ) {
    for( double a=0.0; a<=1; a+=0.1 ) {
      for( double b=0.0; b<=1; b+=0.1 ) {
        PixelLab<double> tmp1(l,a,b);
        PixelRGB<double> tmp2(tmp1);
        PixelLab<double> lab(tmp2);
        EXPECT_PIXEL_NEAR( lab, tmp1, 1e-4 );
      }
    }
  }
}

TEST( PixelTypes, PixelMask ) {
  // Default construction
  {
    PixelMask<PixelGray<vw::uint8> > test;
    EXPECT_EQ( test.valid() , 0 );
  }

  // Implicit construction from scalar
  {
    PixelMask<PixelGray<vw::uint8> > test = 5;
    EXPECT_EQ( test[0] , 5 );
    EXPECT_EQ( test[1] , 255 );
  }

  // Construction from child type
  {
    PixelGray<vw::uint8> g = 5;
    PixelMask<PixelGray<vw::uint8> > test = g;
    EXPECT_EQ( test[0] , 5 );
    EXPECT_EQ( test[1] , 255 );
  }

  // Construction from another PixelMask<> w/ same channel type
  {
    PixelMask<PixelGray<vw::uint8> > gv = 5;
    PixelMask<PixelGray<vw::uint8> > test = gv;
    EXPECT_EQ( test[0] , 5 );
    EXPECT_EQ( test[1] , 255 );
  }

  // Construction from another PixelMask<> w/ different channel type
  {
    PixelMask<PixelGray<vw::uint8> > gv = 5;
    PixelMask<PixelGray<float> > test = channel_cast<float>(gv);
    EXPECT_EQ( test[0] , 5 );
    EXPECT_EQ( test[1] , 1 );
  }

  // Construction from another PixelMask<> w/ an implicit conversion
  {
    PixelGray<vw::uint8> foo = 5;
    PixelRGB<vw::uint8> bar(foo);
    EXPECT_EQ( foo[0] , 5 );
    EXPECT_EQ( bar[0] , 5 );

    PixelMask<PixelGray<vw::uint8> > gv = 5;
    PixelMask<PixelRGB<vw::uint8> > test(gv);
    EXPECT_EQ( gv[0] , 5 );
    EXPECT_EQ( test[0] , 5 );
    EXPECT_EQ( test[3] , 255 );

    gv.invalidate();
    PixelMask<PixelRGB<vw::uint8> > test2(gv);
    EXPECT_EQ( gv[0] , 5 );
    EXPECT_EQ( test2[0] , 5 );
    EXPECT_EQ( test2[3] , 0 );
  }

  // Construction from scalar types
  {
    vw::uint8 foo = 5;
    PixelMask<vw::uint8> gv = foo;
    PixelMask<vw::uint8> test(gv);
    EXPECT_EQ( gv[0] , 5 );
    EXPECT_EQ( test[0] , 5 );
    EXPECT_EQ( test[1] , 255 );

    {
      PixelMask<float> foo(5);
      PixelMask<float> bar(foo);
      EXPECT_EQ( foo[0], 5 );
      EXPECT_EQ( bar[0], 5 );
    }

    {
      PixelMask<PixelRGB<float> > foo(0);
      EXPECT_TRUE( is_valid(foo) );
      EXPECT_EQ( foo[0], 0 );
      EXPECT_EQ( foo[1], 0 );
      EXPECT_EQ( foo[2], 0 );
    }

    // Downcast back to vw::uint8
    vw::uint8 bar = test;
    EXPECT_EQ( bar , 5 );

    test.invalidate();
    bar = test;
    EXPECT_EQ( bar , 5 );

    // The following lines should fail to compile ( throwing a boost
    // static assert error... ) because you should not be able to
    // downcast from a PixelRGB<> to a uint8.
    //PixelMask<PixelRGB<uint8> > downcast_test(4,2,6);
    //bar = downcast_test;
  }

  // Check for pixel transparency
  {
    PixelMask<float> v1(1.0);
    PixelMask<PixelGray<vw::uint8> > v2;

    EXPECT_FALSE( is_transparent(v1) );
    EXPECT_TRUE( is_transparent(v2) );

    v2.validate();
    v1.invalidate();

    EXPECT_TRUE( is_transparent(v1) );
    EXPECT_FALSE( is_transparent(v2) );
  }

  // Check to ensure that operations still occur even when the
  // pixels are masked.
  {
    PixelMask<PixelGray<vw::uint8> > v1(23);
    PixelMask<PixelGray<vw::uint8> > v2(6);
    v1.invalidate();

    EXPECT_TRUE( is_transparent(v1) );
    EXPECT_FALSE( is_transparent(v2) );
    EXPECT_EQ( v1[0] , 23 );
    EXPECT_EQ( v2[0] , 6 );

    PixelMask<PixelGray<vw::uint8> > test = v1 + v2;
    EXPECT_TRUE( is_transparent(test) );
    EXPECT_EQ( test[0] , 29 );
  }

  // Test type traits. Oddly CompoundNumChannels can not be in
  // the EXPECT_EQ with out causing linking errors.
  {
    int test = CompoundNumChannels<PixelMask<PixelGray<vw::uint8> > >::value;
    EXPECT_EQ( test, 2 );
    test =  CompoundNumChannels<PixelMask<PixelRGB<vw::uint8> > >::value;
    EXPECT_EQ( test, 4 );
    test = CompoundNumChannels<PixelMask<PixelRGBA<vw::uint8> > >::value;
    EXPECT_EQ( test, 5 );
    test = CompoundNumChannels<PixelMask<Vector3> >::value;
    EXPECT_EQ( test, 4 );
  }
}

TEST( PixelTypes, ChannelName ) {
#define test1(x) EXPECT_EQ(channel_name_to_enum(channel_type_name(x)), x)
#define test2(x) EXPECT_EQ(std::string(channel_type_name(channel_name_to_enum(x))), std::string(x))

  test1(VW_CHANNEL_BOOL);
  test1(VW_CHANNEL_CHAR);
  test1(VW_CHANNEL_INT8);
  test1(VW_CHANNEL_UINT8);
  test1(VW_CHANNEL_INT16);
  test1(VW_CHANNEL_UINT16);
  test1(VW_CHANNEL_INT32);
  test1(VW_CHANNEL_UINT32);
  test1(VW_CHANNEL_FLOAT16);
  test1(VW_CHANNEL_FLOAT32);
  test1(VW_CHANNEL_INT64);
  test1(VW_CHANNEL_UINT64);
  test1(VW_CHANNEL_FLOAT64);
  test1(VW_CHANNEL_GENERIC_1_BYTE);
  test1(VW_CHANNEL_GENERIC_2_BYTE);
  test1(VW_CHANNEL_GENERIC_4_BYTE);
  test1(VW_CHANNEL_GENERIC_8_BYTE);

  test2("BOOL");
  test2("CHAR");
  test2("INT8");
  test2("UINT8");
  test2("INT16");
  test2("UINT16");
  test2("INT32");
  test2("UINT32");
  test2("FLOAT16");
  test2("FLOAT32");
  test2("INT64");
  test2("UINT64");
  test2("FLOAT64");
  test2("GENERIC_1_BYTE");
  test2("GENERIC_2_BYTE");
  test2("GENERIC_4_BYTE");
  test2("GENERIC_8_BYTE");

#undef test2
#undef test1

}
