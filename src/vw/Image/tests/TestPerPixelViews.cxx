// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


// TestPerPixelView.h
#include <gtest/gtest_VW.h>

#include <vw/Image/PerPixelViews.h>
#include <vw/Image/ImageView.h>
#include <vw/Core/Functors.h>

using namespace vw;

template <template<class> class TraitT, class T>
static bool bool_trait( T const& /*arg*/ ) {
  return TraitT<T>::value;
}

template <class T1, class T2>
static bool has_pixel_type( T2 ) {
  return boost::is_same<T1,typename T2::pixel_type>::value;
}

template <class T>
static T square(T arg) {
  return arg*arg;
}

template <class T>
static T multiply(T arg1, T arg2) {
  return arg1*arg2;
}

TEST( PerPixelView, Unary ) {
  ImageView<float> im(2,2); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4;
  UnaryPerPixelView<ImageView<float>, float(*)(float)> ppv( im, square );
  EXPECT_EQ( ppv.cols(), 2 );
  EXPECT_EQ( ppv.rows(), 2 );
  EXPECT_EQ( ppv.planes(), 1 );

  // Test individual pixel access
  EXPECT_EQ( ppv(0,0), 1 );
  EXPECT_EQ( ppv(1,0), 4 );
  EXPECT_EQ( ppv(0,1), 9 );
  EXPECT_EQ( ppv(1,1), 16 );

  // Test full rasterizaion
  ImageView<float> im2 = ppv;
  EXPECT_EQ( im2.cols(), 2 );
  EXPECT_EQ( im2.rows(), 2 );
  EXPECT_EQ( im2.planes(), 1 );
  EXPECT_EQ( im2(0,0), ppv(0,0) );
  EXPECT_EQ( im2(1,0), ppv(1,0) );
  EXPECT_EQ( im2(0,1), ppv(0,1) );
  EXPECT_EQ( im2(1,1), ppv(1,1) );

  // Test partial rasterization
  ImageView<float> im3(1,2);
  ASSERT_NO_THROW( ppv.rasterize( im3, BBox2i(1,0,1,2) ) );
  EXPECT_EQ( im3(0,0), ppv(1,0) );
  EXPECT_EQ( im3(0,1), ppv(1,1) );
  ImageView<float> im4(2,1);
  ASSERT_NO_THROW( ppv.rasterize( im4, BBox2i(0,1,2,1) ) );
  EXPECT_EQ( im4(0,0), ppv(0,1) );
  EXPECT_EQ( im4(1,0), ppv(1,1) );

  // Test the accessor / generic rasterization
  ImageView<float> im5(2,2);
  vw::rasterize( ppv, im5, BBox2i(0,0,2,2) );
  EXPECT_EQ( im5(0,0), ppv(0,0) );
  EXPECT_EQ( im5(1,0), ppv(1,0) );
  EXPECT_EQ( im5(0,1), ppv(0,1) );
  EXPECT_EQ( im5(1,1), ppv(1,1) );

  // Test the iterator
  ImageView<float>::iterator im2i = im2.begin();
  UnaryPerPixelView<ImageView<float>, float(*)(float)>::iterator ppvi = ppv.begin();
  for( int i=0; i<4; ++i ) {
    EXPECT_EQ( *im2i, *ppvi );
    ASSERT_NO_THROW( ++ppvi );
    ++im2i;
  }

  // Test the types
  ASSERT_TRUE( has_pixel_type<float>( ppv ) );
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>(ppv) );
  ASSERT_TRUE( bool_trait<IsImageView>(ppv) );
}

TEST( PerPixelView, Unary2 ) {
  ImageView<float> im(2,2); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4;
  UnaryPerPixelView<ImageView<float>, ArgNegationFunctor> ppv( im );
  ASSERT_EQ( ppv.cols(), 2 );
  ASSERT_EQ( ppv.rows(), 2 );
  ASSERT_EQ( ppv.planes(), 1 );

  // Test individual pixel access
  EXPECT_EQ( ppv(0,0), -1 );
  EXPECT_EQ( ppv(1,0), -2 );
  EXPECT_EQ( ppv(0,1), -3 );
  EXPECT_EQ( ppv(1,1), -4 );

  // Test optimized rasterizaion
  ImageView<float> im2 = ppv;
  EXPECT_EQ( im2.cols(), 2 );
  EXPECT_EQ( im2.rows(), 2 );
  EXPECT_EQ( im2.planes(), 1 );
  EXPECT_EQ( im2(0,0), ppv(0,0) );
  EXPECT_EQ( im2(1,0), ppv(1,0) );
  EXPECT_EQ( im2(0,1), ppv(0,1) );
  EXPECT_EQ( im2(1,1), ppv(1,1) );

  // Test partial rasterization
  ImageView<float> im3(1,2);
  ASSERT_NO_THROW( ppv.rasterize( im3, BBox2i(1,0,1,2) ) );
  EXPECT_EQ( im3(0,0), ppv(1,0) );
  EXPECT_EQ( im3(0,1), ppv(1,1) );
  ImageView<float> im4(2,1);
  ASSERT_NO_THROW( ppv.rasterize( im4, BBox2i(0,1,2,1) ) );
  EXPECT_EQ( im4(0,0), ppv(0,1) );
  EXPECT_EQ( im4(1,0), ppv(1,1) );

  // Test the accessor / generic rasterization
  ImageView<float> im5(2,2);
  vw::rasterize( ppv, im5, BBox2i(0,0,2,2) );
  EXPECT_EQ( im5(0,0), ppv(0,0) );
  EXPECT_EQ( im5(1,0), ppv(1,0) );
  EXPECT_EQ( im5(0,1), ppv(0,1) );
  EXPECT_EQ( im5(1,1), ppv(1,1) );

  // Test the iterator
  ImageView<float>::iterator im2i = im2.begin();
  UnaryPerPixelView<ImageView<float>, ArgNegationFunctor>::iterator ppvi = ppv.begin();
  for( int i=0; i<4; ++i ) {
    EXPECT_EQ( *im2i, *ppvi );
    ASSERT_NO_THROW( ++ppvi );
    ++im2i;
  }

  // Test the types
  ASSERT_TRUE( has_pixel_type<float>( ppv ) );
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>(ppv) );
  ASSERT_TRUE( bool_trait<IsImageView>(ppv) );
}

TEST( PerPixelView, Binary ) {
  ImageView<float> ima(2,2); ima(0,0)=1; ima(1,0)=2; ima(0,1)=3; ima(1,1)=4;
  ImageView<float> imb(2,2); imb(0,0)=1; imb(1,0)=2; imb(0,1)=-1; imb(1,1)=-2;
  BinaryPerPixelView<ImageView<float>, ImageView<float>, float(*)(float,float)> ppv( ima, imb, multiply );
  ASSERT_EQ( ppv.cols(), 2 );
  ASSERT_EQ( ppv.rows(), 2 );
  ASSERT_EQ( ppv.planes(), 1 );

  // Test individual pixel access
  EXPECT_EQ( ppv(0,0), 1 );
  EXPECT_EQ( ppv(1,0), 4 );
  EXPECT_EQ( ppv(0,1), -3 );
  EXPECT_EQ( ppv(1,1), -8 );

  // Test full rasterizaion
  ImageView<float> im2 = ppv;
  EXPECT_EQ( im2.cols(), 2 );
  EXPECT_EQ( im2.rows(), 2 );
  EXPECT_EQ( im2.planes(), 1 );
  EXPECT_EQ( im2(0,0), ppv(0,0) );
  EXPECT_EQ( im2(1,0), ppv(1,0) );
  EXPECT_EQ( im2(0,1), ppv(0,1) );
  EXPECT_EQ( im2(1,1), ppv(1,1) );

  // Test partial rasterization
  ImageView<float> im3(1,2);
  ASSERT_NO_THROW( ppv.rasterize( im3, BBox2i(1,0,1,2) ) );
  EXPECT_EQ( im3(0,0), ppv(1,0) );
  EXPECT_EQ( im3(0,1), ppv(1,1) );
  ImageView<float> im4(2,1);
  ASSERT_NO_THROW( ppv.rasterize( im4, BBox2i(0,1,2,1) ) );
  EXPECT_EQ( im4(0,0), ppv(0,1) );
  EXPECT_EQ( im4(1,0), ppv(1,1) );

  // Test the accessor / generic rasterization
  ImageView<float> im5(2,2);
  vw::rasterize( ppv, im5, BBox2i(0,0,2,2) );
  EXPECT_EQ( im5(0,0), ppv(0,0) );
  EXPECT_EQ( im5(1,0), ppv(1,0) );
  EXPECT_EQ( im5(0,1), ppv(0,1) );
  EXPECT_EQ( im5(1,1), ppv(1,1) );

  // Test the iterator
  ImageView<float>::iterator im2i = im2.begin();
  BinaryPerPixelView<ImageView<float>, ImageView<float>, float(*)(float,float)>::iterator ppvi = ppv.begin();
  for( int i=0; i<4; ++i ) {
    EXPECT_EQ( *im2i, *ppvi );
    ASSERT_NO_THROW( ++ppvi );
    ++im2i;
  }

  // Test the types
  ASSERT_TRUE( has_pixel_type<float>( ppv ) );
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>(ppv) );
  ASSERT_TRUE( bool_trait<IsImageView>(ppv) );
}

TEST( PerPixelView, Binary2 ) {
  ImageView<float> ima(2,2); ima(0,0)=1; ima(1,0)=2; ima(0,1)=3; ima(1,1)=4;
  ImageView<float> imb(2,2); imb(0,0)=1; imb(1,0)=2; imb(0,1)=-1; imb(1,1)=-2;
  BinaryPerPixelView<ImageView<float>, ImageView<float>, ArgArgSumFunctor> ppv( ima, imb );
  ASSERT_EQ( ppv.cols(), 2 );
  ASSERT_EQ( ppv.rows(), 2 );
  ASSERT_EQ( ppv.planes(), 1 );

  // Test individual pixel access
  EXPECT_EQ( ppv(0,0), 2 );
  EXPECT_EQ( ppv(1,0), 4 );
  EXPECT_EQ( ppv(0,1), 2 );
  EXPECT_EQ( ppv(1,1), 2 );

  // Test optimized rasterizaion
  ImageView<float> im2 = ppv;
  ASSERT_EQ( im2.cols(), 2 );
  ASSERT_EQ( im2.rows(), 2 );
  ASSERT_EQ( im2.planes(), 1 );
  EXPECT_EQ( im2(0,0), ppv(0,0) );
  EXPECT_EQ( im2(1,0), ppv(1,0) );
  EXPECT_EQ( im2(0,1), ppv(0,1) );
  EXPECT_EQ( im2(1,1), ppv(1,1) );

  // Test partial rasterization
  ImageView<float> im3(1,2);
  ASSERT_NO_THROW( ppv.rasterize( im3, BBox2i(1,0,1,2) ) );
  EXPECT_EQ( im3(0,0), ppv(1,0) );
  EXPECT_EQ( im3(0,1), ppv(1,1) );
  ImageView<float> im4(2,1);
  ASSERT_NO_THROW( ppv.rasterize( im4, BBox2i(0,1,2,1) ) );
  EXPECT_EQ( im4(0,0), ppv(0,1) );
  EXPECT_EQ( im4(1,0), ppv(1,1) );

  // Test the accessor / generic rasterization
  ImageView<float> im5(2,2);
  vw::rasterize( ppv, im5, BBox2i(0,0,2,2) );
  EXPECT_EQ( im5(0,0), ppv(0,0) );
  EXPECT_EQ( im5(1,0), ppv(1,0) );
  EXPECT_EQ( im5(0,1), ppv(0,1) );
  EXPECT_EQ( im5(1,1), ppv(1,1) );

  // Test the iterator
  ImageView<float>::iterator im2i = im2.begin();
  BinaryPerPixelView<ImageView<float>, ImageView<float>, ArgArgSumFunctor>::iterator ppvi = ppv.begin();
  for( int i=0; i<4; ++i ) {
    EXPECT_EQ( *im2i, *ppvi );
    ASSERT_NO_THROW( ++ppvi );
    ++im2i;
  }

  // Test the types
  ASSERT_TRUE( has_pixel_type<float>( ppv ) );
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>(ppv) );
  ASSERT_TRUE( bool_trait<IsImageView>(ppv) );
}

/// Simple test function that incorporates the pixel indices
template <class T>
static T addConstant(T arg, int32 i, int32 j, int32 p=0) {
  return arg + (i + 2*j + 4*p);
}


TEST( PerPixelView, UnaryIndex ) {
  ImageView<float> im(2,2); im(0,0)=1; im(0,1)=1; 
                            im(1,0)=1; im(1,1)=1;
  UnaryPerPixelIndexView<ImageView<float>, float(*)(float, int32, int32, int32)> ppv( im, addConstant );

  // Test individual pixel access
  EXPECT_EQ( ppv(0,0), 1 );
  EXPECT_EQ( ppv(1,0), 2 );
  EXPECT_EQ( ppv(0,1), 3 );
  EXPECT_EQ( ppv(1,1), 4 );

  // Test full rasterizaion
  ImageView<float> im2 = ppv;
  EXPECT_EQ( im2.cols(),   2 );
  EXPECT_EQ( im2.rows(),   2 );
  EXPECT_EQ( im2.planes(), 1 );
  EXPECT_EQ( im2(0,0), ppv(0,0) );
  EXPECT_EQ( im2(1,0), ppv(1,0) );
  EXPECT_EQ( im2(0,1), ppv(0,1) );
  EXPECT_EQ( im2(1,1), ppv(1,1) );

  // Test partial rasterization
  ImageView<float> im3(1,2);
  ASSERT_NO_THROW( ppv.rasterize( im3, BBox2i(1,0,1,2) ) );
  EXPECT_EQ( im3(0,0), ppv(1,0) );
  EXPECT_EQ( im3(0,1), ppv(1,1) );
  ImageView<float> im4(2,1);
  ASSERT_NO_THROW( ppv.rasterize( im4, BBox2i(0,1,2,1) ) );
  EXPECT_EQ( im4(0,0), ppv(0,1) );
  EXPECT_EQ( im4(1,0), ppv(1,1) );

  // Test the accessor / generic rasterization
  ImageView<float> im5(2,2);
  vw::rasterize( ppv, im5, BBox2i(0,0,2,2) );
  EXPECT_EQ( im5(0,0), ppv(0,0) );
  EXPECT_EQ( im5(1,0), ppv(1,0) );
  EXPECT_EQ( im5(0,1), ppv(0,1) );
  EXPECT_EQ( im5(1,1), ppv(1,1) );

  // Test the types
  ASSERT_TRUE( has_pixel_type<float>( ppv ) );
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>(ppv) );
  ASSERT_TRUE( bool_trait<IsImageView>(ppv) );
}

