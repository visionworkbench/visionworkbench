// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestUtilityViews.h
#include <gtest/gtest.h>

#include <vw/Image/PixelTypes.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/UtilityViews.h>

// For rand48 (we use rand48 because it's fastest)
#include <boost/random/linear_congruential.hpp>

#include <test/Helpers.h>

using namespace vw;

TEST( UtilityViews, ConstantView ) {
  PixelGray<float32> px(3.5);
  ImageViewRef<PixelGray<float32> > view = constant_view<PixelGray<float32> >(px, 50, 50);

  EXPECT_PIXEL_EQ( view(10, 20), px );
  view = constant_view(px, view);
  EXPECT_PIXEL_EQ( view(10, 20), px );
}

TEST( UtilityViews, VectorIndexView ) {
  ImageViewRef<Vector2> view = pixel_index_view(50, 50);

  EXPECT_VECTOR_DOUBLE_EQ( view(10, 20), Vector2(10, 20) );
  view = pixel_index_view(copy(view));
  EXPECT_VECTOR_DOUBLE_EQ( view(10, 20), Vector2(10, 20) );
}

TEST( UtilityViews, Vector3IndexView ) {
  ImageViewRef<Vector3> view = pixel_index3_view(50, 50, 3);

  EXPECT_VECTOR_DOUBLE_EQ( view(10, 20, 2), Vector3(10, 20, 2) );
  view = pixel_index3_view(view);
  EXPECT_VECTOR_DOUBLE_EQ( view(10, 20, 2), Vector3(10, 20, 2) );
}

TEST( UtilityViews, UniformNoiseView ) {
  boost::rand48 gen(10);

  ImageViewRef<double> view = uniform_noise_view(gen, 10, 10);
  ImageViewRef<double> view2 = uniform_noise_view(gen, view);

  // Insure gen isn't copy constructed anywhere
  EXPECT_NE(view(5, 5), view2(5, 5));

  // Insure each access gives a different value
  EXPECT_NE(view(5, 5), view(5, 5));
}

TEST( UtilityViews, GaussianNoiseView ) {
  boost::rand48 gen(10);

  ImageViewRef<double> view = gaussian_noise_view(gen, 1, 0.5, 10, 10);
  ImageViewRef<double> view2 = gaussian_noise_view(gen, 1, 0.5, view);

  // Insure gen isn't copy constructed anywhere
  EXPECT_NE(view(5, 5), view2(5, 5));

  // Insure each access gives a different value
  EXPECT_NE(view(5, 5), view(5, 5));
}
