// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestInterpolation.h
#include <gtest/gtest.h>

#include <vw/Image/ImageView.h>
#include <vw/Image/Interpolation.h>

using namespace vw;

template <template<class> class TraitT, class T>
static bool bool_trait( T const& /*arg*/ ) {
  return TraitT<T>::value;
}

TEST( Interpolation, Bilinear ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  InterpolationView<EdgeExtensionView<ImageView<double>, ConstantEdgeExtension>, BilinearInterpolation> im2 = interpolate(im, BilinearInterpolation());
  ASSERT_EQ( im2.cols(), 2 );
  ASSERT_EQ( im2.rows(), 3 );
  EXPECT_EQ( im2(1,1), 4 );
  EXPECT_EQ( im2(1,1.5), 5 );
  EXPECT_EQ( im2(0.5,1), 3.5 );
  EXPECT_EQ( im2(0.5,0.5), 2.5 );

  // Teste accessor
  EXPECT_EQ( *(im2.origin().advance(1,1.5)), 5 );
  EXPECT_EQ( *(im2.origin().advance(1,1)), 4 );

  // Test the traits
  ASSERT_TRUE( bool_trait<IsFloatingPointIndexable>( im2 ) );
  ASSERT_TRUE( bool_trait<IsFloatingPointIndexable>( interpolate(im, BilinearInterpolation()) ) );
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( interpolate(im, BilinearInterpolation()) ) );
}

TEST( Interpolation, Bicubic ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  InterpolationView<EdgeExtensionView<ImageView<double>, ZeroEdgeExtension>, BicubicInterpolation> im2 = interpolate(im, BicubicInterpolation(), ZeroEdgeExtension());
  ASSERT_EQ( im2.cols(), 2 );
  ASSERT_EQ( im2.rows(), 3 );
  EXPECT_EQ( im2(1,1), 4 );
  EXPECT_EQ( im2(1,1.5), 5.5 );
  EXPECT_EQ( im2(0.5,1), 3.9375);
  EXPECT_NEAR( im2(0.5,0.5), 2.7773, 0.001);

  // Teste accessor
  EXPECT_EQ( *(im2.origin().advance(1,1.5)), 5.5 );
  EXPECT_EQ( *(im2.origin().advance(1,1)), 4 );

  // Test the traits
  ASSERT_TRUE( bool_trait<IsFloatingPointIndexable>( im2 ) );
  ASSERT_TRUE( bool_trait<IsFloatingPointIndexable>( interpolate(im, BicubicInterpolation()) ) );
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( interpolate(im, BicubicInterpolation()) ) );
}

TEST( Interpolation, Nearest ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  InterpolationView<EdgeExtensionView<ImageView<double>, ConstantEdgeExtension>, NearestPixelInterpolation> im2 = interpolate(im, NearestPixelInterpolation());
  ASSERT_EQ( im2.cols(), 2 );
  ASSERT_EQ( im2.rows(), 3 );
  EXPECT_EQ( im2(1,1), 4 );
  EXPECT_EQ( im2(1,1.5), 6 );
  EXPECT_EQ( im2(1,1.2), 4 );
  EXPECT_EQ( im2(0.7,1), 4 );
  EXPECT_EQ( im2(0.4,1.8), 5 );

  // Teste accessor
  EXPECT_EQ( *(im2.origin().advance(1,1.5)), 6 );
  EXPECT_EQ( *(im2.origin().advance(1,1)), 4 );

  // Test the traits
  ASSERT_TRUE( bool_trait<IsFloatingPointIndexable>( im2 ) );
  ASSERT_TRUE( bool_trait<IsFloatingPointIndexable>( interpolate(im, NearestPixelInterpolation()) ) );
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( interpolate(im, NearestPixelInterpolation()) ) );
}

