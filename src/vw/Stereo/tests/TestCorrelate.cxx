// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestDisparity.h
#include <gtest/gtest.h>

#include <vw/Stereo/Correlate.h>
#include <vw/Image.h>

using namespace vw;
using namespace vw::stereo;

typedef PixelMask<Vector2f> PixelDisp;

TEST( Correlate, CrossCorrConsistency ) {

  ImageView<PixelDisp> r2l(3,3), l2r(3,3);
  fill( r2l, PixelDisp(Vector2f(0,0)) );
  fill( l2r, PixelDisp(Vector2f(0,0)) );

  fill( crop(l2r,2,0,1,3), PixelDisp(Vector2f(2,2)) );
  l2r(0,0) = PixelDisp(Vector2f(1,1));
  r2l(1,1) = PixelDisp(Vector2f(-1,-1));
  l2r(1,0) = PixelDisp(Vector2f(1,1));

  ImageView<PixelDisp> l2r_copy = copy(l2r);
  cross_corr_consistency_check( l2r_copy, r2l, 0 );
  EXPECT_FALSE( is_valid(l2r_copy(2,0)) );
  EXPECT_FALSE( is_valid(l2r_copy(2,1)) );
  EXPECT_FALSE( is_valid(l2r_copy(2,2)) );
  EXPECT_TRUE( is_valid( l2r_copy(0,0)) );
  EXPECT_FALSE( is_valid(l2r_copy(1,0)) );

  l2r_copy = l2r;
  cross_corr_consistency_check( l2r_copy, r2l, 2 );
  EXPECT_FALSE( is_valid(l2r_copy(2,0)) );
  EXPECT_FALSE( is_valid(l2r_copy(2,1)) );
  EXPECT_FALSE( is_valid(l2r_copy(2,2)) );
  EXPECT_TRUE( is_valid( l2r_copy(0,0)) );
  EXPECT_TRUE( is_valid( l2r_copy(1,0)) );
}
