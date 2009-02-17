// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestPerPixelView.h
#include <cxxtest/TestSuite.h>

#include <vw/Image/PerPixelViews.h>
#include <vw/Image/ImageView.h>
#include <vw/Core/Functors.h>

using namespace vw;

class TestPerPixelViews : public CxxTest::TestSuite
{
public:

  template <template<class> class TraitT, class T>
  static bool bool_trait( T const& arg ) {
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

  void test_unary_per_pixel_view_function() {
    ImageView<float> im(2,2); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4;
    UnaryPerPixelView<ImageView<float>, float(*)(float)> ppv( im, square );
    TS_ASSERT_EQUALS( ppv.cols(), 2 );
    TS_ASSERT_EQUALS( ppv.rows(), 2 );
    TS_ASSERT_EQUALS( ppv.planes(), 1 );

    // Test individual pixel access
    TS_ASSERT_EQUALS( ppv(0,0), 1 );
    TS_ASSERT_EQUALS( ppv(1,0), 4 );
    TS_ASSERT_EQUALS( ppv(0,1), 9 );
    TS_ASSERT_EQUALS( ppv(1,1), 16 );

    // Test full rasterizaion
    ImageView<float> im2 = ppv;
    TS_ASSERT_EQUALS( im2.cols(), 2 );
    TS_ASSERT_EQUALS( im2.rows(), 2 );
    TS_ASSERT_EQUALS( im2.planes(), 1 );
    TS_ASSERT_EQUALS( im2(0,0), ppv(0,0) );
    TS_ASSERT_EQUALS( im2(1,0), ppv(1,0) );
    TS_ASSERT_EQUALS( im2(0,1), ppv(0,1) );
    TS_ASSERT_EQUALS( im2(1,1), ppv(1,1) );

    // Test partial rasterization
    ImageView<float> im3(1,2);
    TS_ASSERT_THROWS_NOTHING( ppv.rasterize( im3, BBox2i(1,0,1,2) ) );
    TS_ASSERT_EQUALS( im3(0,0), ppv(1,0) );
    TS_ASSERT_EQUALS( im3(0,1), ppv(1,1) );
    ImageView<float> im4(2,1);
    TS_ASSERT_THROWS_NOTHING( ppv.rasterize( im4, BBox2i(0,1,2,1) ) );
    TS_ASSERT_EQUALS( im4(0,0), ppv(0,1) );
    TS_ASSERT_EQUALS( im4(1,0), ppv(1,1) );

    // Test the accessor / generic rasterization
    ImageView<float> im5(2,2);
    vw::rasterize( ppv, im5, BBox2i(0,0,2,2) );
    TS_ASSERT_EQUALS( im5(0,0), ppv(0,0) );
    TS_ASSERT_EQUALS( im5(1,0), ppv(1,0) );
    TS_ASSERT_EQUALS( im5(0,1), ppv(0,1) );
    TS_ASSERT_EQUALS( im5(1,1), ppv(1,1) );

    // Test the iterator
    ImageView<float>::iterator im2i = im2.begin();
    UnaryPerPixelView<ImageView<float>, float(*)(float)>::iterator ppvi = ppv.begin();
    for( int i=0; i<4; ++i ) {
      TS_ASSERT_DIFFERS( ppvi, ppv.end() );
      TS_ASSERT_EQUALS( *im2i, *ppvi );
      TS_ASSERT_THROWS_NOTHING( ++ppvi );
      ++im2i;
    }
    TS_ASSERT_EQUALS( ppvi, ppv.end() );

    // Test the types
    TS_ASSERT( has_pixel_type<float>( ppv ) );
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>(ppv) );
    TS_ASSERT( bool_trait<IsImageView>(ppv) );
  }

  void test_unary_per_pixel_view_functor() {
    ImageView<float> im(2,2); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4;
    UnaryPerPixelView<ImageView<float>, ArgNegationFunctor> ppv( im );
    TS_ASSERT_EQUALS( ppv.cols(), 2 );
    TS_ASSERT_EQUALS( ppv.rows(), 2 );
    TS_ASSERT_EQUALS( ppv.planes(), 1 );

    // Test individual pixel access
    TS_ASSERT_EQUALS( ppv(0,0), -1 );
    TS_ASSERT_EQUALS( ppv(1,0), -2 );
    TS_ASSERT_EQUALS( ppv(0,1), -3 );
    TS_ASSERT_EQUALS( ppv(1,1), -4 );

    // Test optimized rasterizaion
    ImageView<float> im2 = ppv;
    TS_ASSERT_EQUALS( im2.cols(), 2 );
    TS_ASSERT_EQUALS( im2.rows(), 2 );
    TS_ASSERT_EQUALS( im2.planes(), 1 );
    TS_ASSERT_EQUALS( im2(0,0), ppv(0,0) );
    TS_ASSERT_EQUALS( im2(1,0), ppv(1,0) );
    TS_ASSERT_EQUALS( im2(0,1), ppv(0,1) );
    TS_ASSERT_EQUALS( im2(1,1), ppv(1,1) );

    // Test partial rasterization
    ImageView<float> im3(1,2);
    TS_ASSERT_THROWS_NOTHING( ppv.rasterize( im3, BBox2i(1,0,1,2) ) );
    TS_ASSERT_EQUALS( im3(0,0), ppv(1,0) );
    TS_ASSERT_EQUALS( im3(0,1), ppv(1,1) );
    ImageView<float> im4(2,1);
    TS_ASSERT_THROWS_NOTHING( ppv.rasterize( im4, BBox2i(0,1,2,1) ) );
    TS_ASSERT_EQUALS( im4(0,0), ppv(0,1) );
    TS_ASSERT_EQUALS( im4(1,0), ppv(1,1) );

    // Test the accessor / generic rasterization
    ImageView<float> im5(2,2);
    vw::rasterize( ppv, im5, BBox2i(0,0,2,2) );
    TS_ASSERT_EQUALS( im5(0,0), ppv(0,0) );
    TS_ASSERT_EQUALS( im5(1,0), ppv(1,0) );
    TS_ASSERT_EQUALS( im5(0,1), ppv(0,1) );
    TS_ASSERT_EQUALS( im5(1,1), ppv(1,1) );

    // Test the iterator
    ImageView<float>::iterator im2i = im2.begin();
    UnaryPerPixelView<ImageView<float>, ArgNegationFunctor>::iterator ppvi = ppv.begin();
    for( int i=0; i<4; ++i ) {
      TS_ASSERT_DIFFERS( ppvi, ppv.end() );
      TS_ASSERT_EQUALS( *im2i, *ppvi );
      TS_ASSERT_THROWS_NOTHING( ++ppvi );
      ++im2i;
    }
    TS_ASSERT_EQUALS( ppvi, ppv.end() );

    // Test the types
    TS_ASSERT( has_pixel_type<float>( ppv ) );
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>(ppv) );
    TS_ASSERT( bool_trait<IsImageView>(ppv) );
  }

  void test_binary_per_pixel_view_function() {
    ImageView<float> ima(2,2); ima(0,0)=1; ima(1,0)=2; ima(0,1)=3; ima(1,1)=4;
    ImageView<float> imb(2,2); imb(0,0)=1; imb(1,0)=2; imb(0,1)=-1; imb(1,1)=-2;
    BinaryPerPixelView<ImageView<float>, ImageView<float>, float(*)(float,float)> ppv( ima, imb, multiply );
    TS_ASSERT_EQUALS( ppv.cols(), 2 );
    TS_ASSERT_EQUALS( ppv.rows(), 2 );
    TS_ASSERT_EQUALS( ppv.planes(), 1 );

    // Test individual pixel access
    TS_ASSERT_EQUALS( ppv(0,0), 1 );
    TS_ASSERT_EQUALS( ppv(1,0), 4 );
    TS_ASSERT_EQUALS( ppv(0,1), -3 );
    TS_ASSERT_EQUALS( ppv(1,1), -8 );

    // Test full rasterizaion
    ImageView<float> im2 = ppv;
    TS_ASSERT_EQUALS( im2.cols(), 2 );
    TS_ASSERT_EQUALS( im2.rows(), 2 );
    TS_ASSERT_EQUALS( im2.planes(), 1 );
    TS_ASSERT_EQUALS( im2(0,0), ppv(0,0) );
    TS_ASSERT_EQUALS( im2(1,0), ppv(1,0) );
    TS_ASSERT_EQUALS( im2(0,1), ppv(0,1) );
    TS_ASSERT_EQUALS( im2(1,1), ppv(1,1) );

    // Test partial rasterization
    ImageView<float> im3(1,2);
    TS_ASSERT_THROWS_NOTHING( ppv.rasterize( im3, BBox2i(1,0,1,2) ) );
    TS_ASSERT_EQUALS( im3(0,0), ppv(1,0) );
    TS_ASSERT_EQUALS( im3(0,1), ppv(1,1) );
    ImageView<float> im4(2,1);
    TS_ASSERT_THROWS_NOTHING( ppv.rasterize( im4, BBox2i(0,1,2,1) ) );
    TS_ASSERT_EQUALS( im4(0,0), ppv(0,1) );
    TS_ASSERT_EQUALS( im4(1,0), ppv(1,1) );

    // Test the accessor / generic rasterization
    ImageView<float> im5(2,2);
    vw::rasterize( ppv, im5, BBox2i(0,0,2,2) );
    TS_ASSERT_EQUALS( im5(0,0), ppv(0,0) );
    TS_ASSERT_EQUALS( im5(1,0), ppv(1,0) );
    TS_ASSERT_EQUALS( im5(0,1), ppv(0,1) );
    TS_ASSERT_EQUALS( im5(1,1), ppv(1,1) );

    // Test the iterator
    ImageView<float>::iterator im2i = im2.begin();
    BinaryPerPixelView<ImageView<float>, ImageView<float>, float(*)(float,float)>::iterator ppvi = ppv.begin();
    for( int i=0; i<4; ++i ) {
      TS_ASSERT_DIFFERS( ppvi, ppv.end() );
      TS_ASSERT_EQUALS( *im2i, *ppvi );
      TS_ASSERT_THROWS_NOTHING( ++ppvi );
      ++im2i;
    }
    TS_ASSERT_EQUALS( ppvi, ppv.end() );

    // Test the types
    TS_ASSERT( has_pixel_type<float>( ppv ) );
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>(ppv) );
    TS_ASSERT( bool_trait<IsImageView>(ppv) );
  }

  void test_binary_per_pixel_view_functor() {
    ImageView<float> ima(2,2); ima(0,0)=1; ima(1,0)=2; ima(0,1)=3; ima(1,1)=4;
    ImageView<float> imb(2,2); imb(0,0)=1; imb(1,0)=2; imb(0,1)=-1; imb(1,1)=-2;
    BinaryPerPixelView<ImageView<float>, ImageView<float>, ArgArgSumFunctor> ppv( ima, imb );
    TS_ASSERT_EQUALS( ppv.cols(), 2 );
    TS_ASSERT_EQUALS( ppv.rows(), 2 );
    TS_ASSERT_EQUALS( ppv.planes(), 1 );

    // Test individual pixel access
    TS_ASSERT_EQUALS( ppv(0,0), 2 );
    TS_ASSERT_EQUALS( ppv(1,0), 4 );
    TS_ASSERT_EQUALS( ppv(0,1), 2 );
    TS_ASSERT_EQUALS( ppv(1,1), 2 );

    // Test optimized rasterizaion
    ImageView<float> im2 = ppv;
    TS_ASSERT_EQUALS( im2.cols(), 2 );
    TS_ASSERT_EQUALS( im2.rows(), 2 );
    TS_ASSERT_EQUALS( im2.planes(), 1 );
    TS_ASSERT_EQUALS( im2(0,0), ppv(0,0) );
    TS_ASSERT_EQUALS( im2(1,0), ppv(1,0) );
    TS_ASSERT_EQUALS( im2(0,1), ppv(0,1) );
    TS_ASSERT_EQUALS( im2(1,1), ppv(1,1) );

    // Test partial rasterization
    ImageView<float> im3(1,2);
    TS_ASSERT_THROWS_NOTHING( ppv.rasterize( im3, BBox2i(1,0,1,2) ) );
    TS_ASSERT_EQUALS( im3(0,0), ppv(1,0) );
    TS_ASSERT_EQUALS( im3(0,1), ppv(1,1) );
    ImageView<float> im4(2,1);
    TS_ASSERT_THROWS_NOTHING( ppv.rasterize( im4, BBox2i(0,1,2,1) ) );
    TS_ASSERT_EQUALS( im4(0,0), ppv(0,1) );
    TS_ASSERT_EQUALS( im4(1,0), ppv(1,1) );

    // Test the accessor / generic rasterization
    ImageView<float> im5(2,2);
    vw::rasterize( ppv, im5, BBox2i(0,0,2,2) );
    TS_ASSERT_EQUALS( im5(0,0), ppv(0,0) );
    TS_ASSERT_EQUALS( im5(1,0), ppv(1,0) );
    TS_ASSERT_EQUALS( im5(0,1), ppv(0,1) );
    TS_ASSERT_EQUALS( im5(1,1), ppv(1,1) );

    // Test the iterator
    ImageView<float>::iterator im2i = im2.begin();
    BinaryPerPixelView<ImageView<float>, ImageView<float>, ArgArgSumFunctor>::iterator ppvi = ppv.begin();
    for( int i=0; i<4; ++i ) {
      TS_ASSERT_DIFFERS( ppvi, ppv.end() );
      TS_ASSERT_EQUALS( *im2i, *ppvi );
      TS_ASSERT_THROWS_NOTHING( ++ppvi );
      ++im2i;
    }
    TS_ASSERT_EQUALS( ppvi, ppv.end() );

    // Test the types
    TS_ASSERT( has_pixel_type<float>( ppv ) );
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>(ppv) );
    TS_ASSERT( bool_trait<IsImageView>(ppv) );
  }

};
