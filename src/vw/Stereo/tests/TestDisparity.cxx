
// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestDisparity.h
#include <gtest/gtest.h>

#include <vw/Stereo/DisparityMap.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/Transform.h>
#include <test/Helpers.h>

using namespace vw;
using namespace vw::stereo;

typedef PixelMask<Vector2f> PixelDisp;

TEST( DisparityMap, Transform1 ) {

  vw::Matrix<double,3,3> align_matrix;
  align_matrix.set_identity();
  // Adding just a translation.
  align_matrix(0,2) = 45;
  align_matrix(1,2) = -30;

  // Building disparity map
  ImageView<PixelDisp > map(5,5);
  for (size_t i = 0; i < 5; i++)
    for (size_t j = 0; j < 5; j++)
      map(i,j) = PixelDisp(i*5+j,j*7+i);

  // Applying the inverse of the align matrix
  ImageViewRef<PixelDisp > result;
  result = transform_disparities(map, HomographyTransform(align_matrix));

  // Comparing results
  for (size_t i = 0; i < 5; i++)
    for (size_t j = 0; j < 5; j++) {
      Vector3 t_disparity(result(i,j)[0], result(i,j)[1], 1);
      Vector3 location(i,j,0);
      Vector3 check = align_matrix*(t_disparity + location) - location;
      EXPECT_VECTOR_NEAR(map(i,j).child(),
                         subvector(check,0,2), 1e-1);
    }
}

TEST( DisparityMap, Transform2 ) {

  vw::Matrix<double,3,3> align_matrix;
  align_matrix.set_identity();
  align_matrix(0,0) = 1.00679;
  align_matrix(0,1) = -0.0125401;
  align_matrix(0,2) = 116.812;
  align_matrix(1,0) = 0.00788373;
  align_matrix(1,1) = .996033;
  align_matrix(1,2) = -1.93039; //Homography affine

  // Building disparity map
  ImageView<PixelDisp > map(5,5);
  for (size_t i = 0; i < 5; i++)
    for (size_t j = 0; j < 5; j++)
      map(i,j) = PixelDisp(i*5+j,j*7+i);

  // Applying the inverse of the align matrix
  ImageViewRef<PixelDisp > result;
  result = transform_disparities(map, HomographyTransform(align_matrix));

  // Comparing results
  for (size_t i = 0; i < 5; i++)
    for (size_t j = 0; j < 5; j++) {
      Vector3 t_disparity(result(i,j)[0], result(i,j)[1], 1);
      Vector3 location(i,j,0);
      Vector3 check = align_matrix*(location + t_disparity) - location;
      EXPECT_VECTOR_NEAR(map(i,j).child(),
                         subvector(check,0,2), 1e-1);
    }
}

TEST( DisparityMap, DisparitySubsample ) {
  ImageView<PixelDisp > map(4,4);
  map(0,0) = PixelDisp(Vector2f(3,1));
  map(2,0) = PixelDisp(Vector2f(4,2));
  map(2,1) = PixelDisp(Vector2f(2,2));
  ASSERT_TRUE( is_valid(map(0,0)) );
  ASSERT_FALSE( is_valid(map(1,1)) );

  ImageView<PixelDisp > submap =
    disparity_subsample( map );
  ASSERT_EQ( submap.cols(), 2 );
  ASSERT_EQ( submap.rows(), 2 );
  EXPECT_TRUE( is_valid(submap(0,0)) );
  EXPECT_TRUE( is_valid(submap(1,0)) );
  EXPECT_FALSE( is_valid(submap(0,1)) );
  EXPECT_TRUE( is_valid(submap(1,1)) );
  EXPECT_VECTOR_NEAR( submap(0,0).child(),
                      Vector2f(1.5,0.5), 1e-3 );
  EXPECT_VECTOR_NEAR( submap(1,0).child(),
                      Vector2f(1.75,1.0), 1e-3 );
  EXPECT_VECTOR_NEAR( submap(1,1).child(),
                      Vector2f(1,1), 1e-3 );

  ImageViewRef<PixelDisp > submapref =
    disparity_subsample( ImageViewRef<PixelDisp >(map) );
  ASSERT_EQ( submapref.cols(), 2 );
  ASSERT_EQ( submapref.rows(), 2 );
  EXPECT_TRUE( is_valid(submapref(0,0)) );
  EXPECT_TRUE( is_valid(submapref(1,0)) );
  EXPECT_FALSE( is_valid(submapref(0,1)) );
  EXPECT_TRUE( is_valid(submapref(1,1)) );
  EXPECT_VECTOR_NEAR( submapref(0,0).child(),
                      Vector2f(1.5,0.5), 1e-3 );
  EXPECT_VECTOR_NEAR( submapref(1,0).child(),
                      Vector2f(1.75,1.0), 1e-3 );
  EXPECT_VECTOR_NEAR( submapref(1,1).child(),
                      Vector2f(1,1), 1e-3 );

  ImageView<PixelMask<Vector2i> > imap(3,1);
  imap(0,0) = PixelMask<Vector2i>(Vector2i(4,2));
  imap(1,0) = PixelMask<Vector2i>(Vector2i(10,-8));
  ImageViewRef<PixelMask<Vector2i> > simap =
    disparity_subsample( imap );
  ASSERT_EQ( simap.cols(), 2 );
  ASSERT_EQ( simap.rows(), 1 );
  EXPECT_TRUE( is_valid(simap(0,0)) );
  EXPECT_TRUE( is_valid(simap(1,0)) );
  EXPECT_VECTOR_NEAR( simap(0,0).child(),
                      Vector2i(2, 0), 1e-3 );
  EXPECT_VECTOR_NEAR( simap(1,0).child(),
                      Vector2i(5,-4), 1e-3 );
}

TEST( DisparityMap, DisparityUpsample ) {
  ImageView<PixelDisp > map(2,2);
  map(0,0) = PixelDisp(Vector2f(3,1));
  map(1,1) = PixelDisp(Vector2f(5,5));

  ImageView<PixelDisp > upmap =
    disparity_upsample( map );
  ASSERT_EQ( upmap.cols(), 4 );
  ASSERT_EQ( upmap.rows(), 4 );
  for (size_t i = 0; i < 4; i++ )
    EXPECT_TRUE( is_valid(upmap(i,i)) );
  for (size_t i = 0; i < 4; i++ )
    EXPECT_FALSE( is_valid(upmap(3-i,i)) );
  EXPECT_VECTOR_NEAR( upmap(0,0).child(),
                      Vector2f(6,2), 1e-3 );
  EXPECT_VECTOR_NEAR( upmap(1,1).child(),
                      Vector2f(6,2), 1e-3 );
  EXPECT_VECTOR_NEAR( upmap(2,2).child(),
                      Vector2f(10,10), 1e-3 );
  EXPECT_VECTOR_NEAR( upmap(3,2).child(),
                      Vector2f(10,10), 1e-3 );

  ImageViewRef<PixelDisp > upmapref =
    disparity_upsample( ImageViewRef<PixelDisp >(map) );
  ASSERT_EQ( upmapref.cols(), 4 );
  ASSERT_EQ( upmapref.rows(), 4 );
  for (size_t i = 0; i < 4; i++ )
    EXPECT_TRUE( is_valid(upmapref(i,i)) );
  for (size_t i = 0; i < 4; i++ )
    EXPECT_FALSE( is_valid(upmapref(3-i,i)) );
  EXPECT_VECTOR_NEAR( upmapref(0,0).child(),
                      Vector2f(6,2), 1e-3 );
  EXPECT_VECTOR_NEAR( upmapref(1,1).child(),
                      Vector2f(6,2), 1e-3 );
  EXPECT_VECTOR_NEAR( upmapref(2,2).child(),
                      Vector2f(10,10), 1e-3 );
  EXPECT_VECTOR_NEAR( upmapref(3,2).child(),
                      Vector2f(10,10), 1e-3 );
}

TEST( DisparityMap, DisparityTransform ) {
  // Disparity map will be delta.x = 2 + (left.x - 50) * 0.1
  ImageView<PixelDisp > disparity(100,1);
  for ( int32 i = 0; i < 100; i++ ) {
    float delta = 2 + (i - 50)*0.1;
    disparity(i,0) = PixelDisp(Vector2f(delta,0));
  }

  //Test Disparity Transform
  DisparityTransform trans( disparity );
  for ( int32 i = 0; i < 100; i++ ) {
    float delta = 2 + (i - 50)*0.1;
    EXPECT_NEAR( delta+float(i), trans.reverse(Vector2(i,0))[0], 1e-5 );
  }
}

TEST( DisparityMap, GetDisparityRange ) {
  ImageView<PixelDisp> disparity(4,1);
  disparity(0,0) = PixelDisp(Vector2f(2,2));
  disparity(1,0) = PixelDisp(Vector2f(3,5));
  disparity(2,0) = PixelDisp(Vector2f(-4,-1));
  disparity(2,0).invalidate();

  BBox2f range = get_disparity_range(disparity);
  EXPECT_VECTOR_EQ( Vector2f(2,2), range.min() );
  EXPECT_VECTOR_EQ( Vector2f(3,5), range.max() );

  disparity(0,0).invalidate(); disparity(1,0).invalidate();
  range = get_disparity_range(disparity);
  EXPECT_VECTOR_EQ( Vector2f(), range.min() );
  EXPECT_VECTOR_EQ( Vector2f(), range.max() );
}
