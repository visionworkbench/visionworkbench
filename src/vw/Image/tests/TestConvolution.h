// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
//
// __END_LICENSE__

// TestConvolution.h
#include <cxxtest/TestSuite.h>

#include <vw/Image/Convolution.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>

#include <vector>

using namespace std;
using namespace vw;

class TestConvolution : public CxxTest::TestSuite
{
public:

  template <template<class> class TraitT, class T>
  static bool bool_trait( T const& arg ) {
    return TraitT<T>::value;
  }

  template <class T1, class T2>
  static bool is_of_type( T2 ) {
    return boost::is_same<T1,T2>::value;
  }

  template <class T1, class T2>
  static bool has_pixel_type( T2 ) {
    return boost::is_same<T1,typename T2::pixel_type>::value;
  }

  void test_correlate_1d_at_point() {
    ImageView<double> src(3,1); src(0,0)=1; src(1,0)=2; src(2,0)=3;
    double kernel[] = {2,3,-1};
    TS_ASSERT_EQUALS( correlate_1d_at_point(src.origin(),(double*)kernel,3), 5 );
    TS_ASSERT( is_of_type<double>( correlate_1d_at_point( ImageView<double>().origin(), (double*)0, 0 ) ) );
    TS_ASSERT( is_of_type<PixelRGB<float> >( correlate_1d_at_point( ImageView<PixelRGB<uint8> >().origin(), (float*)0, 0 ) ) );
  }

  void test_correlate_2d_at_point() {
    ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
    ImageView<double> krn(2,2); krn(0,0)=2; krn(1,0)=-1; krn(0,1)=0; krn(1,1)=3;
    TS_ASSERT_EQUALS( correlate_2d_at_point(src.origin(),krn.origin(),2,2), 12 );
    TS_ASSERT( is_of_type<double>( correlate_2d_at_point( ImageView<double>().origin(), ImageView<double>().origin(), 0, 0 ) ) );
    TS_ASSERT( is_of_type<PixelRGB<float> >( correlate_2d_at_point( ImageView<PixelRGB<uint8> >().origin(), ImageView<float>().origin(), 0, 0 ) ) );
  }

  void test_ConvolutionView() {
    ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
    ImageView<double> krn(2,2); krn(0,0)=2; krn(1,0)=-1; krn(0,1)=0; krn(1,1)=3;
    ConvolutionView<ImageView<double>,ImageView<double>,ZeroEdgeExtend> cnv( src, krn );

    // Test individual pixel access
    TS_ASSERT_EQUALS( cnv.cols(), 2 );
    TS_ASSERT_EQUALS( cnv.rows(), 2 );
    TS_ASSERT_EQUALS( cnv.planes(), 1 );
    TS_ASSERT_EQUALS( cnv(0,0), 2 );
    TS_ASSERT_EQUALS( cnv(1,0), 3 );
    TS_ASSERT_EQUALS( cnv(0,1), 6 );
    TS_ASSERT_EQUALS( cnv(1,1), 8 );

    // Test rasterization
    ImageView<double> dst = cnv;
    TS_ASSERT_EQUALS( dst.cols(), 2 );
    TS_ASSERT_EQUALS( dst.rows(), 2 );
    TS_ASSERT_EQUALS( dst(0,0), 2 );
    TS_ASSERT_EQUALS( dst(1,0), 3 );
    TS_ASSERT_EQUALS( dst(0,1), 6 );
    TS_ASSERT_EQUALS( dst(1,1), 8 );

    // Test the traits
    TS_ASSERT( is_of_type<double>( cnv(0,0) ) );
    TS_ASSERT( is_of_type<double>( *(cnv.origin()) ) );
    TS_ASSERT( !bool_trait<IsResizable>( cnv ) );
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( cnv ) );
    TS_ASSERT( !bool_trait<IsFloatingPointIndexable>( cnv ) );
    TS_ASSERT( bool_trait<IsImageView>( cnv ) );
  }

  void test_SeparableConvolutionView() {
    ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
    vector<double> krn; krn.push_back(1); krn.push_back(-1);
    SeparableConvolutionView<ImageView<double>,double,ZeroEdgeExtend> cnv( src, krn, krn );

    // Test individual pixel access
    TS_ASSERT_EQUALS( cnv.cols(), 2 );
    TS_ASSERT_EQUALS( cnv.rows(), 2 );
    TS_ASSERT_EQUALS( cnv.planes(), 1 );
    TS_ASSERT_EQUALS( cnv(0,0), 1 );
    TS_ASSERT_EQUALS( cnv(0,1), 2 );
    TS_ASSERT_EQUALS( cnv(1,0), 1 );
    TS_ASSERT_EQUALS( cnv(1,1), 0 );

    // Test rasterization
    ImageView<double> dst = cnv;
    TS_ASSERT_EQUALS( dst.cols(), 2 );
    TS_ASSERT_EQUALS( dst.rows(), 2 );
    TS_ASSERT_EQUALS( dst.planes(), 1 );
    TS_ASSERT_EQUALS( dst(0,0), 1 );
    TS_ASSERT_EQUALS( dst(0,1), 2 );
    TS_ASSERT_EQUALS( dst(1,0), 1 );
    TS_ASSERT_EQUALS( dst(1,1), 0 );

    // Test the traits
    TS_ASSERT( is_of_type<double>( cnv(0,0) ) );
    TS_ASSERT( is_of_type<double>( *(cnv.origin()) ) );
    TS_ASSERT( !bool_trait<IsResizable>( cnv ) );
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( cnv ) );
    TS_ASSERT( !bool_trait<IsFloatingPointIndexable>( cnv ) );
    TS_ASSERT( bool_trait<IsImageView>( cnv ) );
  }

  void test_SeparableConvolutionView_0x2() {
    ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=4; src(1,1)=6;
    vector<double> krnx;
    vector<double> krny; krny.push_back(1); krny.push_back(-1);
    SeparableConvolutionView<ImageView<double>,double,ZeroEdgeExtend> cnv( src, krnx, krny );

    // Test individual pixel access
    TS_ASSERT_EQUALS( cnv.cols(), 2 );
    TS_ASSERT_EQUALS( cnv.rows(), 2 );
    TS_ASSERT_EQUALS( cnv.planes(), 1 );
    TS_ASSERT_EQUALS( cnv(0,0), 1 );
    TS_ASSERT_EQUALS( cnv(0,1), 3 );
    TS_ASSERT_EQUALS( cnv(1,0), 2 );
    TS_ASSERT_EQUALS( cnv(1,1), 4 );

    // Test rasterization
    ImageView<double> dst = cnv;
    TS_ASSERT_EQUALS( dst.cols(), 2 );
    TS_ASSERT_EQUALS( dst.rows(), 2 );
    TS_ASSERT_EQUALS( dst.planes(), 1 );
    TS_ASSERT_EQUALS( dst(0,0), 1 );
    TS_ASSERT_EQUALS( dst(0,1), 3 );
    TS_ASSERT_EQUALS( dst(1,0), 2 );
    TS_ASSERT_EQUALS( dst(1,1), 4 );

    // Test the traits
    TS_ASSERT( is_of_type<double>( cnv(0,0) ) );
    TS_ASSERT( is_of_type<double>( *(cnv.origin()) ) );
    TS_ASSERT( !bool_trait<IsResizable>( cnv ) );
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( cnv ) );
    TS_ASSERT( !bool_trait<IsFloatingPointIndexable>( cnv ) );
    TS_ASSERT( bool_trait<IsImageView>( cnv ) );
  }

  void test_SeparableConvolutionView_2x0() {
    ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
    vector<double> krnx; krnx.push_back(1); krnx.push_back(-1);
    vector<double> krny;
    SeparableConvolutionView<ImageView<double>,double,ZeroEdgeExtend> cnv( src, krnx, krny );

    // Test individual pixel access
    TS_ASSERT_EQUALS( cnv.cols(), 2 );
    TS_ASSERT_EQUALS( cnv.rows(), 2 );
    TS_ASSERT_EQUALS( cnv.planes(), 1 );
    TS_ASSERT_EQUALS( cnv(0,0), 1 );
    TS_ASSERT_EQUALS( cnv(1,0), 1 );
    TS_ASSERT_EQUALS( cnv(0,1), 3 );
    TS_ASSERT_EQUALS( cnv(1,1), 1 );

    // Test rasterization
    ImageView<double> dst = cnv;
    TS_ASSERT_EQUALS( dst.cols(), 2 );
    TS_ASSERT_EQUALS( dst.rows(), 2 );
    TS_ASSERT_EQUALS( dst.planes(), 1 );
    TS_ASSERT_EQUALS( dst(0,0), 1 );
    TS_ASSERT_EQUALS( dst(1,0), 1 );
    TS_ASSERT_EQUALS( dst(0,1), 3 );
    TS_ASSERT_EQUALS( dst(1,1), 1 );

    // Test the traits
    TS_ASSERT( is_of_type<double>( cnv(0,0) ) );
    TS_ASSERT( is_of_type<double>( *(cnv.origin()) ) );
    TS_ASSERT( !bool_trait<IsResizable>( cnv ) );
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( cnv ) );
    TS_ASSERT( !bool_trait<IsFloatingPointIndexable>( cnv ) );
    TS_ASSERT( bool_trait<IsImageView>( cnv ) );
  }

}; // class TestConvolution
