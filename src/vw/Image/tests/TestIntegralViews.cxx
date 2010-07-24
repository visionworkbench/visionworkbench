// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestIntegralViews.h
#include <gtest/gtest.h>

#include <vw/Image/IntegralView.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/Algorithms.h>
#include <iostream>
#include <stdio.h>

#include <test/Helpers.h>

using namespace std;
using namespace vw;

typedef ImageView<double> image_t;

// | 1 2 |
// | 1 4 |
// | 1 1 |
image_t test_image() {
  image_t im(2,3);
  im(0,0)=1;
  im(1,0)=1;
  im(0,1)=2;
  im(1,1)=4;
  im(0,2)=1;
  im(1,2)=1;
  return im;
}

TEST( IntegralViews, BlockSumViewBasic ){
  image_t im(test_image());
  BlockSumView<image_t> cv(im);

  ASSERT_EQ( im.cols(),   cv.cols() );
  ASSERT_EQ( im.rows(),   cv.rows() );
  ASSERT_EQ( im.planes(), cv.planes() );
}

TEST( IntegralViews, BlockSumViewAccess ) {

  image_t im1 = test_image();
  BlockSumView<image_t> cv(im1);

  //test full rasterization & pixel access
  image_t im2 = block_sum(im1,3);
  ASSERT_EQ( im2.cols(),   cv.cols() );
  ASSERT_EQ( im2.rows(),   cv.rows() );
  ASSERT_EQ( im2.planes(), cv.planes() );

  for ( int r=0; r < im2.rows(); ++r)
    for ( int c=0; c < im2.cols(); ++c)
      EXPECT_EQ( im2(c,r), cv(c,r) ) << " at (c,r) = (" << c << "," << r << ")";

  //test partial rasterization
  image_t im3(im1.cols()-1,im1.rows()-1);
  ASSERT_NO_THROW( cv.rasterize(im3, BBox2i(1,1,cv.cols()-1,cv.rows()-1) ) );
  for ( int r=0; r < im3.rows(); ++r)
    for ( int c=0; c < im3.cols(); ++c)
      EXPECT_EQ( im3(c,r), cv(c+1,r+1) ) << " at (c,r) = (" << c << "," << r << ")";
}

// | 1 1 1 1 |     | 1 2  3  4 |
// | 1 1 1 1 |  I  | 2 4  6  8 |
// | 1 1 1 1 | ->  | 3 6  9 12 | out(i,j) = (i+1) * (j+1)
// | 1 1 1 1 |     | 4 8 12 16 |
TEST( IntegralViews, IntegralImage ) {
  image_t im(4,4);
  vw::fill(im, 1);

  image_t result = integral_image(im);
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      EXPECT_EQ( (i+1) * (j+1), result(i,j)) << " at (i,j) = (" << i << "," << j << ")";
}

TEST( IntegralViews, IntegralEdgeExtension ) {
  image_t im(4,4);
  vw::fill(im, 1);
  typedef EdgeExtensionView<image_t,IntegralEdgeExtension> result_type;
  result_type result =
    edge_extend(integral_image(im), IntegralEdgeExtension());

  EXPECT_EQ( 0, result(-1,-1) );
  EXPECT_EQ( 0, result(0,-1) );
  EXPECT_EQ( 0, result(-1,0) );
  EXPECT_EQ( 0, result(-2,-4) );
  EXPECT_EQ( 1, result(0,0) );
  EXPECT_EQ( 4, result(1,1) );
  EXPECT_EQ( 9, result(2,2) );
  EXPECT_EQ( 16, result(3,3) );
  EXPECT_EQ( 16, result(4,4) );
  EXPECT_EQ( 16, result(5,5) );
  EXPECT_EQ( 12, result(2,5) );
  EXPECT_EQ( 12, result(5,2) );
  EXPECT_EQ( 4,  result(8,0) );

  image_t prerastered = edge_extend(integral_image(im),
                                    BBox2i(-1,-1,6,6),
                                    IntegralEdgeExtension());
  EXPECT_EQ( 6, prerastered.cols() );
  EXPECT_EQ( 6, prerastered.rows() );
  for ( unsigned i = 0; i < 6; i++ ) {
    EXPECT_EQ( 0, prerastered(0,i) ) << " at (0," << i << ")";
    EXPECT_EQ( 0, prerastered(i,0) ) << " at (" << i << ",0)";
  }
}

// | 1 1 1 1 |     | 1 2  3  4 |     | 4 6 6 4 |
// | 1 1 1 1 |  I  | 2 4  6  8 |  B  | 6 9 9 6 |
// | 1 1 1 1 | ->  | 3 6  9 12 | ->  | 6 9 9 6 |
// | 1 1 1 1 |     | 4 8 12 16 |     | 4 6 6 4 |
TEST( IntegralViews, BlockSum ) {
  image_t im(10,10);
  vw::fill(im, 1);

  image_t expected(10,10);
  vw::fill(expected, 9);
  vw::fill(crop(expected,0,1,1,8), 6);
  vw::fill(crop(expected,1,0,8,1), 6);
  vw::fill(crop(expected,9,1,1,8), 6);
  vw::fill(crop(expected,1,9,8,1), 6);
  expected(0,0) = expected(9,9) = expected(9,0) = expected(0,9) = 4;

  image_t result = block_sum(im,3);
  for ( unsigned i = 0; i < 10; i++ ) {
    for ( unsigned j = 0; j < 10; j++ ) {
      EXPECT_EQ( expected(i,j), result(i,j) ) << " at (i,j) = (" << i << "," << j << ")";
    }
  }
}

