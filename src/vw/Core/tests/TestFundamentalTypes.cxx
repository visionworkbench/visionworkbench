// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>

#include <vw/Core/FundamentalTypes.h>
#include <sstream>

template <typename T>
T smallest() {return vw::ScalarTypeLimits<T>::smallest();}

template <typename T>
T lowest() {return vw::ScalarTypeLimits<T>::lowest();}

template <typename T>
T highest() {return vw::ScalarTypeLimits<T>::highest();}

TEST(FundamentalTypes, Basic) {
  vw::FundamentalTypeClass<vw::int32> t, t2(42);
  EXPECT_EQ( 0, t );
  EXPECT_EQ( 42, t2 );
}

TEST(FundamentalTypes, Ranges) {
  // TODO: commented out some checks because there's no reliable way to declare
  // a 64-bit literal. this should probably be extended to check whether
  // -ffast-math (... broken math) is enabled.
  EXPECT_EQ( 1, smallest<vw::int8>() );
  EXPECT_EQ( 1, smallest<vw::int16>() );
  EXPECT_EQ( 1, smallest<vw::int32>() );
  EXPECT_EQ( 1, smallest<vw::int64>() );
  EXPECT_EQ( 1u, smallest<vw::uint8>() );
  EXPECT_EQ( 1u, smallest<vw::uint16>() );
  EXPECT_EQ( 1u, smallest<vw::uint32>() );
  EXPECT_EQ( 1u, smallest<vw::uint64>() );
  //EXPECT_FLOAT_EQ( , smallest<vw::float32>() );
  //EXPECT_DOUBLE_EQ( , smallest<vw::float64>() );

  EXPECT_EQ( -128L        , lowest<vw::int8>() );
  EXPECT_EQ( -32768L      , lowest<vw::int16>() );
  EXPECT_EQ( -2147483648L , lowest<vw::int32>() );
  //EXPECT_EQ(            , lowest<vw::int64>() );
  EXPECT_EQ( 0UL          , lowest<vw::uint8>() );
  EXPECT_EQ( 0UL          , lowest<vw::uint16>() );
  EXPECT_EQ( 0UL          , lowest<vw::uint32>() );
  EXPECT_EQ( 0UL          , lowest<vw::uint64>() );
  //EXPECT_FLOAT_EQ( 0    , lowest<vw::float32>() );
  //EXPECT_DOUBLE_EQ( 0   , lowest<vw::float64>() );

  EXPECT_EQ( 127L         , highest<vw::int8>() );
  EXPECT_EQ( 32767L       , highest<vw::int16>() );
  EXPECT_EQ( 2147483647L  , highest<vw::int32>() );
  //EXPECT_EQ(            , highest<vw::int64>() );
  EXPECT_EQ( 255UL        , highest<vw::uint8>() );
  EXPECT_EQ( 65535UL      , highest<vw::uint16>() );
  EXPECT_EQ( 4294967295UL , highest<vw::uint32>() );
  //EXPECT_EQ(            , highest<vw::uint64>() );
  //EXPECT_FLOAT_EQ(      , highest<vw::float32>() );
  //EXPECT_DOUBLE_EQ(     , highest<vw::float64>() );
}

TEST(FundamentalTypes, Numeric) {
  std::stringstream s;
  s << vw::_numeric(vw::uint8(42)) << vw::_numeric(vw::int8(-42));
  EXPECT_STREQ( "42-42", s.str().c_str() );
  EXPECT_FLOAT_EQ( 4.2f, vw::_numeric(4.2f) );
}
