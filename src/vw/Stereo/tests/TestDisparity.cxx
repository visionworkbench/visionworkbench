// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestDisparity.h
#include <gtest/gtest.h>

#include <vw/Math.h>
#include <vw/Stereo.h>
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


