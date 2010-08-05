// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

// TestBoxFilter.h
#include <gtest/gtest.h>

#include <vw/InterestPoint/IntegralImage.h>
#include <vw/InterestPoint/BoxFilter.h>
#include <vw/Image/ImageView.h>

using namespace vw;
using namespace vw::ip;

TEST( BoxFilter, ApplyBoxFilter ) {
  ImageView<float> image(3,3);

  float count = 0;
  for ( uint i = 0; i < 3; i++ ) {
    for ( uint j = 0; j < 3; j++ ) {
      image(i,j) = count;
      count++;
    }
  }

  ImageView<float> integral = IntegralImage( image );
  EXPECT_EQ( 4, integral.cols() );
  EXPECT_EQ( 4, integral.rows() );

  BoxFilter filter;
  filter.resize(2);
  filter[0].start = Vector2i(0,0);
  filter[0].size = Vector2i(1,1);
  filter[0].weight = -9;
  filter[1].start = Vector2i(-1,-1);
  filter[1].size = Vector2i(3,3);
  filter[1].weight = 1;

  float response = apply_box_filter_at_point( integral.origin().advance(1,1),
                                              filter );
  EXPECT_NEAR( 0, response, 1e-5 );
}

TEST( BoxFilter, FilterView ) {
  ImageView<float> image(4,4);
  float count = 0;
  for ( uint8 i = 0; i < 4; i++ ) {
    for ( uint8 j = 0; j < 4; j++ ) {
      image(i,j) = count;
      count++;
    }
  }
  image(3,3) = 20;
  ImageView<float> integral = IntegralImage( image );
  EXPECT_EQ( 5, integral.cols() );
  EXPECT_EQ( 5, integral.rows() );

  BoxFilter filter;
  filter.resize(2);
  filter[0].start = Vector2i(0,0);
  filter[0].size = Vector2i(1,1);
  filter[0].weight = -9;
  filter[1].start = Vector2i(-1,-1);
  filter[1].size = Vector2i(3,3);
  filter[1].weight = 1;

  ImageView<float> applied = box_filter( integral, filter );
  EXPECT_EQ( 4, applied.cols() );
  EXPECT_EQ( 4, applied.rows() );

  EXPECT_NEAR( 0, applied(1,1), 1e-5 );
  EXPECT_NEAR( 5, applied(2,2), 1e-5 );
  EXPECT_NEAR( 0, applied(1,2), 1e-5 );
  EXPECT_NEAR( 0, applied(2,1), 1e-5 );
  EXPECT_NEAR( 0, applied(0,0), 1e-5 );
  EXPECT_NEAR( 0, applied(3,3), 1e-5 );
}
