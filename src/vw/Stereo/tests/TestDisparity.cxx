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

TEST( DisparityMap, Transform1 ) {

  vw::Matrix<double,3,3> align_matrix;
  align_matrix.set_identity();
  // Adding just a translation.
  align_matrix(0,2) = 45;
  align_matrix(1,2) = -30;

  // Building disparity map
  ImageView<PixelMask<Vector2f> > map(5,5);
  for (unsigned i = 0; i < 5; i++)
    for (unsigned j = 0; j < 5; j++)
      map(i,j) = PixelMask<Vector2f>(i*5+j,j*7+i);

  // Applying the inverse of the align matrix
  ImageViewRef<PixelMask<Vector2f> > result;
  result = transform_disparities(map, HomographyTransform(align_matrix));

  // Comparing results
  for (unsigned i = 0; i < 5; i++)
    for (unsigned j = 0; j < 5; j++) {
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
  ImageView<PixelMask<Vector2f> > map(5,5);
  for (unsigned i = 0; i < 5; i++)
    for (unsigned j = 0; j < 5; j++)
      map(i,j) = PixelMask<Vector2f>(i*5+j,j*7+i);

  // Applying the inverse of the align matrix
  ImageViewRef<PixelMask<Vector2f> > result;
  result = transform_disparities(map, HomographyTransform(align_matrix));

  // Comparing results
  for (unsigned i = 0; i < 5; i++)
    for (unsigned j = 0; j < 5; j++) {
      Vector3 t_disparity(result(i,j)[0], result(i,j)[1], 1);
      Vector3 location(i,j,0);
      Vector3 check = align_matrix*(location + t_disparity) - location;
      EXPECT_VECTOR_NEAR(map(i,j).child(),
                         subvector(check,0,2), 1e-1);
    }
}

TEST( DisparityMap, DisparitySubsample ) {
  ImageView<PixelMask<Vector2f> > map(4,4);
  map(0,0) = PixelMask<Vector2f>(Vector2f(3,1));
  map(2,0) = PixelMask<Vector2f>(Vector2f(4,2));
  map(2,1) = PixelMask<Vector2f>(Vector2f(2,2));

  ImageView<PixelMask<Vector2f> > submap =
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
}

TEST( DisparityMap, DisparityUpsample ) {
  ImageView<PixelMask<Vector2f> > map(2,2);
  map(0,0) = PixelMask<Vector2f>(Vector2f(3,1));
  map(1,1) = PixelMask<Vector2f>(Vector2f(5,5));

  ImageView<PixelMask<Vector2f> > upmap =
    disparity_upsample( map );
  ASSERT_EQ( upmap.cols(), 4 );
  ASSERT_EQ( upmap.rows(), 4 );
  for (unsigned i = 0; i < 4; i++ )
    EXPECT_TRUE( is_valid(upmap(i,i)) );
  for (unsigned i = 0; i < 4; i++ )
    EXPECT_FALSE( is_valid(upmap(3-i,i)) );
  EXPECT_VECTOR_NEAR( upmap(0,0).child(),
                      Vector2f(6,2), 1e-3 );
  EXPECT_VECTOR_NEAR( upmap(1,1).child(),
                      Vector2f(6,2), 1e-3 );
  EXPECT_VECTOR_NEAR( upmap(2,2).child(),
                      Vector2f(10,10), 1e-3 );
  EXPECT_VECTOR_NEAR( upmap(3,2).child(),
                      Vector2f(10,10), 1e-3 );
}

TEST( DisparityMap, DisparityTransform ) {
  // Disparity map will be delta.x = 2 + (left.x - 50) * 0.1
  ImageView<PixelMask<Vector2f> > disparity(100,1);
  for ( int32 i = 0; i < 100; i++ ) {
    float delta = 2 + (i - 50)*0.1;
    disparity(i,0) = PixelMask<Vector2f>(Vector2f(delta,0));
  }

  //Test Disparity Transform
  DisparityTransform trans( disparity );
  for ( int32 i = 0; i < 100; i++ ) {
    float delta = 2 + (i - 50)*0.1;
    EXPECT_NEAR( delta+float(i), trans.reverse(Vector2(i,0))[0], 1e-5 );
  }
}
