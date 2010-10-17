// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestEdgeExtension.h
#include <gtest/gtest.h>

#include <vw/Image/ImageView.h>
#include <vw/Image/EdgeExtension.h>

#include <boost/utility/result_of.hpp>
#include <boost/type_traits.hpp>

using namespace vw;

#define EXPECT_BBOX(B,X,Y,W,H)    \
  { BBox2i b = B;                    \
    EXPECT_EQ(b.min().x(),X); \
    EXPECT_EQ(b.min().y(),Y); \
    EXPECT_EQ(b.width(),W);   \
    EXPECT_EQ(b.height(),H);  \
  }

template <template<class> class TraitT, class T>
static bool bool_trait( T const& /*arg*/ ) {
  return TraitT<T>::value;
}

class SomeType {};

TEST( EdgeExtension, No ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  EdgeExtensionView<ImageView<double>, NoEdgeExtension> im2 = edge_extend(im, NoEdgeExtension() );
  EXPECT_EQ( im2.cols(), 2 );
  EXPECT_EQ( im2.rows(), 3 );
  EXPECT_EQ( im2(1,1), 4 );

  // Test the accessor
  EXPECT_EQ( *(im2.origin().advance(1,1)), 4 );

  // Test the traits
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( edge_extend(im, ZeroEdgeExtension() ) ) );
  ASSERT_TRUE( (boost::is_same<boost::result_of<NoEdgeExtension(ImageView<SomeType>,int,int,int)>::type,SomeType>::value) );

  // Test the prerasterization bbox
  NoEdgeExtension ee;
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(0,0,2,1)), 0,0,2,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,1,1,2)), 1,1,1,2 );
}

TEST( EdgeExtension, Zero ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  EdgeExtensionView<ImageView<double>, ZeroEdgeExtension> im2 = edge_extend(im, ZeroEdgeExtension() );
  EXPECT_EQ( im2.cols(), 2 );
  EXPECT_EQ( im2.rows(), 3 );

  EXPECT_EQ( im2(-1,0), 0 );
  EXPECT_EQ( im2(0,-1), 0 );
  EXPECT_EQ( im2(2,0), 0 );
  EXPECT_EQ( im2(0,3), 0 );
  EXPECT_EQ( im2(-1,-1), 0 );
  EXPECT_EQ( im2(2,3), 0 );

  EXPECT_EQ( im2(5,0), 0 );
  EXPECT_EQ( im2(0,6), 0 );
  EXPECT_EQ( im2(-4,0), 0 );
  EXPECT_EQ( im2(1,-4), 0 );

  // Test the accessor
  EXPECT_EQ( *(im2.origin().advance(-1,-1)), 0 );
  EXPECT_EQ( *(im2.origin().advance(1,1)), 4 );

  // Test the traits
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( edge_extend(im, ZeroEdgeExtension() ) ) );
  ASSERT_TRUE( (boost::is_same<boost::result_of<ZeroEdgeExtension(ImageView<SomeType>,int,int,int)>::type,SomeType>::value) );

  // Test the prerasterization bbox
  ZeroEdgeExtension ee;
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(0,0,2,1)), 0,0,2,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,1,1,2)), 1,1,1,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,2,2)), 0,0,1,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,-1,2,2)), 1,0,1,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-1,2,2,2)), 0,2,1,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,2,2,2)), 1,2,1,1 );
  ASSERT_TRUE( ee.source_bbox(im,BBox2i(-2,-2,2,2)).empty() );
  ASSERT_TRUE( ee.source_bbox(im,BBox2i(2,3,2,2)).empty() );
}

TEST( EdgeExtension, Constant ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  EdgeExtensionView<ImageView<double>, ConstantEdgeExtension> im2 = edge_extend(im, ConstantEdgeExtension() );
  EXPECT_EQ( im2.cols(), 2 );
  EXPECT_EQ( im2.rows(), 3 );

  EXPECT_EQ( im2(-1,0), 1 );
  EXPECT_EQ( im2(1,-1), 2 );
  EXPECT_EQ( im2(2,0), 2 );
  EXPECT_EQ( im2(0,3), 5 );
  EXPECT_EQ( im2(-1,-1), 1 );
  EXPECT_EQ( im2(2,3), 6 );

  EXPECT_EQ( im2(5,0), 2 );
  EXPECT_EQ( im2(0,6), 5 );
  EXPECT_EQ( im2(-4,0), 1 );
  EXPECT_EQ( im2(1,-4), 2 );

  // Test the accessor
  EXPECT_EQ( *(im2.origin().advance(-1,-1)), 1 );
  EXPECT_EQ( *(im2.origin().advance(1,1)), 4 );

  // Test the traits
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( edge_extend(im, ZeroEdgeExtension() ) ) );
  ASSERT_TRUE( (boost::is_same<boost::result_of<ConstantEdgeExtension(ImageView<SomeType>,int,int,int)>::type,SomeType>::value) );

  // Test the prerasterization bbox
  ConstantEdgeExtension ee;
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(0,0,2,1)), 0,0,2,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,1,1,2)), 1,1,1,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,2,2)), 0,0,1,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,-1,2,2)), 1,0,1,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-1,2,2,2)), 0,2,1,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,2,2,2)), 1,2,1,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-2,-2,2,2)), 0,0,1,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(2,3,2,2)), 1,2,1,1 );
}

TEST( EdgeExtension, Periodic ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  EdgeExtensionView<ImageView<double>, PeriodicEdgeExtension> im2 = edge_extend(im, PeriodicEdgeExtension() );
  EXPECT_EQ( im2.cols(), 2 );
  EXPECT_EQ( im2.rows(), 3 );

  EXPECT_EQ( im2(-1,0), 2 );
  EXPECT_EQ( im2(1,-1), 6 );
  EXPECT_EQ( im2(2,0), 1 );
  EXPECT_EQ( im2(0,3), 1 );
  EXPECT_EQ( im2(-1,-1), 6 );
  EXPECT_EQ( im2(2,3), 1 );

  EXPECT_EQ( im2(3,0), 2 );
  EXPECT_EQ( im2(4,0), 1 );
  EXPECT_EQ( im2(5,0), 2 );

  EXPECT_EQ( im2(0,4), 3 );
  EXPECT_EQ( im2(0,5), 5 );
  EXPECT_EQ( im2(0,6), 1 );

  EXPECT_EQ( im2(-2,0), 1 );
  EXPECT_EQ( im2(-3,0), 2 );
  EXPECT_EQ( im2(-4,0), 1 );

  EXPECT_EQ( im2(1,-2), 4 );
  EXPECT_EQ( im2(1,-3), 2 );
  EXPECT_EQ( im2(1,-4), 6 );

  // Test the accessor
  EXPECT_EQ( *(im2.origin().advance(-1,-1)), 6 );
  EXPECT_EQ( *(im2.origin().advance(1,1)), 4 );

  // Test the traits
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( edge_extend(im, PeriodicEdgeExtension() ) ) );
  ASSERT_TRUE( (boost::is_same<boost::result_of<PeriodicEdgeExtension(ImageView<SomeType>,int,int,int)>::type,SomeType>::value) );

  // Test the prerasterization bbox
  PeriodicEdgeExtension ee;
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(0,0,2,1)), 0,0,2,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,1,1,2)), 1,1,1,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,2,2)), 0,0,2,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,-1,2,2)), 0,0,2,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-1,2,2,2)), 0,0,2,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,2,2,2)), 0,0,2,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-2,-2,2,2)), 0,1,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(2,3,2,2)), 0,0,2,2 );
  im.set_size(4,4);
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-5,-5,2,2)), 0,0,4,4 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-4,-4,2,2)), 0,0,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-3,-3,2,2)), 1,1,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-2,-2,2,2)), 2,2,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,2,2)), 0,0,4,4 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(0,0,2,2)), 0,0,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,1,2,2)), 1,1,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(2,2,2,2)), 2,2,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(3,3,2,2)), 0,0,4,4 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(4,4,2,2)), 0,0,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(5,5,2,2)), 1,1,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-5,-5,3,3)), 0,0,4,4 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-4,-4,3,3)), 0,0,3,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-3,-3,3,3)), 1,1,3,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-2,-2,3,3)), 0,0,4,4 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,3,3)), 0,0,4,4 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(0,0,3,3)), 0,0,3,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,1,3,3)), 1,1,3,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(2,2,3,3)), 0,0,4,4 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(3,3,3,3)), 0,0,4,4 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(4,4,3,3)), 0,0,3,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(5,5,3,3)), 1,1,3,3 );
}

TEST( EdgeExtension, Cylindrical ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  EdgeExtensionView<ImageView<double>, CylindricalEdgeExtension> im2 = edge_extend(im, CylindricalEdgeExtension() );
  EXPECT_EQ( im2.cols(), 2 );
  EXPECT_EQ( im2.rows(), 3 );

  EXPECT_EQ( im2(-1,0), 2 );
  EXPECT_EQ( im2(1,-1), 2 );
  EXPECT_EQ( im2(2,0), 1 );
  EXPECT_EQ( im2(0,3), 5 );
  EXPECT_EQ( im2(-1,-1), 2 );
  EXPECT_EQ( im2(2,3), 5 );

  EXPECT_EQ( im2(3,0), 2 );
  EXPECT_EQ( im2(4,0), 1 );
  EXPECT_EQ( im2(5,0), 2 );

  EXPECT_EQ( im2(0,4), 5 );
  EXPECT_EQ( im2(0,5), 5 );
  EXPECT_EQ( im2(0,6), 5 );

  EXPECT_EQ( im2(-2,0), 1 );
  EXPECT_EQ( im2(-3,0), 2 );
  EXPECT_EQ( im2(-4,0), 1 );

  EXPECT_EQ( im2(1,-2), 2 );
  EXPECT_EQ( im2(1,-3), 2 );
  EXPECT_EQ( im2(1,-4), 2 );

  // Test the accessor
  EXPECT_EQ( *(im2.origin().advance(-1,-1)), 2 );
  EXPECT_EQ( *(im2.origin().advance(1,1)), 4 );

  // Test the traits
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( edge_extend(im, CylindricalEdgeExtension() ) ) );
  ASSERT_TRUE( (boost::is_same<boost::result_of<CylindricalEdgeExtension(ImageView<SomeType>,int,int,int)>::type,SomeType>::value) );

  // Test the prerasterization bbox
  CylindricalEdgeExtension ee;
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(0,0,2,1)), 0,0,2,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,1,1,2)), 1,1,1,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,2,2)), 0,0,2,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,-1,2,2)), 0,0,2,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-1,2,2,2)), 0,2,2,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,2,2,2)), 0,2,2,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-2,-2,2,2)), 0,0,2,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(2,3,2,2)), 0,2,2,1 );
  im.set_size(4,4);
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-5,-5,2,2)), 0,0,4,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-4,-4,2,2)), 0,0,2,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-3,-3,2,2)), 1,0,2,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-2,-2,2,2)), 2,0,2,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,2,2)), 0,0,4,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(0,0,2,2)), 0,0,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,1,2,2)), 1,1,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(2,2,2,2)), 2,2,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(3,3,2,2)), 0,3,4,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(4,4,2,2)), 0,3,2,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(5,5,2,2)), 1,3,2,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-5,-5,3,3)), 0,0,4,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-4,-4,3,3)), 0,0,3,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-3,-3,3,3)), 1,0,3,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-2,-2,3,3)), 0,0,4,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,3,3)), 0,0,4,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(0,0,3,3)), 0,0,3,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,1,3,3)), 1,1,3,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(2,2,3,3)), 0,2,4,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(3,3,3,3)), 0,3,4,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(4,4,3,3)), 0,3,3,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(5,5,3,3)), 1,3,3,1 );
}

TEST( EdgeExtension, Reflect ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  EdgeExtensionView<ImageView<double>, ReflectEdgeExtension> im2 = edge_extend(im, ReflectEdgeExtension() );
  EXPECT_EQ( im2.cols(), 2 );
  EXPECT_EQ( im2.rows(), 3 );

  EXPECT_EQ( im2(-1,0), 2 );
  EXPECT_EQ( im2(1,-1), 4 );
  EXPECT_EQ( im2(2,0), 1 );
  EXPECT_EQ( im2(0,3), 3 );
  EXPECT_EQ( im2(-1,-1), 4 );
  EXPECT_EQ( im2(2,3), 3 );

  EXPECT_EQ( im2(3,0), 2 );
  EXPECT_EQ( im2(4,0), 1 );
  EXPECT_EQ( im2(5,0), 2 );

  EXPECT_EQ( im2(0,4), 1 );
  EXPECT_EQ( im2(0,5), 3 );
  EXPECT_EQ( im2(0,6), 5 );

  EXPECT_EQ( im2(-2,0), 1 );
  EXPECT_EQ( im2(-3,0), 2 );
  EXPECT_EQ( im2(-4,0), 1 );

  EXPECT_EQ( im2(1,-2), 6 );
  EXPECT_EQ( im2(1,-3), 4 );
  EXPECT_EQ( im2(1,-4), 2 );

  // Test the accessor
  EXPECT_EQ( *(im2.origin().advance(-1,-1)), 4 );
  EXPECT_EQ( *(im2.origin().advance(1,1)), 4 );

  // Test the traits
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( edge_extend(im, ReflectEdgeExtension() ) ) );
  ASSERT_TRUE( (boost::is_same<boost::result_of<ReflectEdgeExtension(ImageView<SomeType>,int,int,int)>::type,SomeType>::value) );

  // Test the prerasterization bbox
  ReflectEdgeExtension ee;
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(0,0,2,1)), 0,0,2,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,1,1,2)), 1,1,1,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,2,2)), 0,0,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,-1,2,2)), 0,0,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-1,2,2,2)), 0,1,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,2,2,2)), 0,1,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-2,-2,2,2)), 0,1,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(2,3,2,2)), 0,0,2,2 );
  im.set_size(4,4);
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-5,-5,2,2)), 1,1,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-4,-4,2,2)), 2,2,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-3,-3,2,2)), 2,2,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-2,-2,2,2)), 1,1,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,2,2)), 0,0,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(0,0,2,2)), 0,0,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,1,2,2)), 1,1,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(2,2,2,2)), 2,2,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(3,3,2,2)), 2,2,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(4,4,2,2)), 1,1,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(5,5,2,2)), 0,0,2,2 );

  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-5,-5,3,3)), 1,1,3,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-4,-4,3,3)), 2,2,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-3,-3,3,3)), 1,1,3,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-2,-2,3,3)), 0,0,3,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,3,3)), 0,0,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(0,0,3,3)), 0,0,3,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,1,3,3)), 1,1,3,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(2,2,3,3)), 2,2,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(3,3,3,3)), 1,1,3,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(4,4,3,3)), 0,0,3,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(5,5,3,3)), 0,0,2,2 );

  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-5,-5,4,4)), 1,1,3,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-4,-4,4,4)), 1,1,3,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-3,-3,4,4)), 0,0,4,4 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-2,-2,4,4)), 0,0,3,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,4,4)), 0,0,3,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(0,0,4,4)), 0,0,4,4 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,1,4,4)), 1,1,3,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(2,2,4,4)), 1,1,3,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(3,3,4,4)), 0,0,4,4 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(4,4,4,4)), 0,0,3,3 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(5,5,4,4)), 0,0,3,3 );

  im.set_size(5,5);
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,7,7)), 0,0,5,5 );
}

TEST( EdgeExtension, Linear ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  EdgeExtensionView<ImageView<double>, LinearEdgeExtension> im2 = edge_extend(im, LinearEdgeExtension() );
  EXPECT_EQ( im2.cols(), 2 );
  EXPECT_EQ( im2.rows(), 3 );

  EXPECT_EQ( im2(-1,0), 0 );
  EXPECT_EQ( im2(1,-1), 0 );
  EXPECT_EQ( im2(2,0), 3 );
  EXPECT_EQ( im2(0,3), 7 );
  EXPECT_EQ( im2(-1,-1), -2 );
  EXPECT_EQ( im2(2,3), 9 );

  EXPECT_EQ( im2(3,0), 4 );
  EXPECT_EQ( im2(4,0), 5 );
  EXPECT_EQ( im2(5,0), 6 );

  EXPECT_EQ( im2(0,4), 9 );
  EXPECT_EQ( im2(0,5), 11 );
  EXPECT_EQ( im2(0,6), 13 );

  EXPECT_EQ( im2(-2,0), -1 );
  EXPECT_EQ( im2(-3,0), -2 );
  EXPECT_EQ( im2(-4,0), -3 );

  EXPECT_EQ( im2(1,-2), -2 );
  EXPECT_EQ( im2(1,-3), -4 );
  EXPECT_EQ( im2(1,-4), -6 );

  // Test the accessor
  EXPECT_EQ( *(im2.origin().advance(-1,-1)), -2 );
  EXPECT_EQ( *(im2.origin().advance(1,1)), 4 );

  // Test the traits
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( edge_extend(im, LinearEdgeExtension() ) ) );
  ASSERT_TRUE( (boost::is_same<boost::result_of<LinearEdgeExtension(ImageView<SomeType>,int,int,int)>::type,SomeType>::value) );

  // Test the prerasterization bbox
  LinearEdgeExtension ee;
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(0,0,2,1)), 0,0,2,1 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,1,1,2)), 1,1,1,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,2,2)), 0,0,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,-1,2,2)), 0,0,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-1,2,2,2)), 0,1,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(1,2,2,2)), 0,1,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(-2,-2,2,2)), 0,0,2,2 );
  EXPECT_BBOX( ee.source_bbox(im,BBox2i(2,3,2,2)), 0,1,2,2 );
}

template <class PixelT>
class FloatingView : public ImageViewBase<FloatingView<PixelT> > {
  int32 m_cols, m_rows, m_planes;
public:
  typedef PixelT pixel_type;
  typedef PixelT result_type;
  typedef ProceduralPixelAccessor<FloatingView> pixel_accessor;

  FloatingView( int32 cols, int32 rows, int32 planes = 1 )
    : m_cols(cols), m_rows(rows), m_planes(planes) {}

  inline int32 cols() const { return m_cols; }
  inline int32 rows() const { return m_rows; }
  inline int32 planes() const { return m_planes; }

  inline pixel_accessor origin() const { return pixel_accessor( *this ); }

  inline result_type operator()( double col, double row, double /*plane*/=0 ) const { return col*row; }

  typedef FloatingView prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i /*bbox*/ ) const { return *this; }
  template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const {
    vw::rasterize( prerasterize(bbox), dest, bbox );
  }
};

namespace vw {
template <class PixelT>
struct IsMultiplyAccessible<FloatingView<PixelT> > : public true_type {};

template <class PixelT>
struct IsFloatingPointIndexable<FloatingView<PixelT> > : public true_type {};
}

template <class PixelT>
FloatingView<PixelT> floating_view( PixelT const& /*value*/, int32 cols, int32 rows, int32 planes=1 ) {
  return FloatingView<PixelT>( cols, rows, planes );
}
