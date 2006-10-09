// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
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

// TestFilter.h
#include <cxxtest/TestSuite.h>

#include <vw/Image/Filter.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>

#include <vector>

using namespace std;
using namespace vw;

class TestFilter : public CxxTest::TestSuite
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

  void test_generate_gaussian_kernel() {
    std::vector<double> kernel;
    generate_gaussian_kernel( kernel, 1.0, 5 );
    TS_ASSERT_EQUALS( kernel.size(), 5 );
    TS_ASSERT_DELTA( kernel[0], 0.06135958087, 1e-8 );
    TS_ASSERT_DELTA( kernel[1], 0.2447702197, 1e-7 );
    TS_ASSERT_DELTA( kernel[2], 0.3877403988, 1e-7 );
    TS_ASSERT_DELTA( kernel[3], 0.2447702197, 1e-7 );
    TS_ASSERT_DELTA( kernel[4], 0.06135958087, 1e-8 );
    generate_gaussian_kernel( kernel, 1.0, 4 );
    TS_ASSERT_EQUALS( kernel.size(), 4 );
    TS_ASSERT_DELTA( kernel[0], 0.1423836140, 1e-7 );
    TS_ASSERT_DELTA( kernel[1], 0.3576163860, 1e-7 );
    TS_ASSERT_DELTA( kernel[2], 0.3576163860, 1e-7 );
    TS_ASSERT_DELTA( kernel[3], 0.1423836140, 1e-7 );
    generate_gaussian_kernel( kernel, 1.5 );
    TS_ASSERT_EQUALS( kernel.size(), 9 );
    TS_ASSERT_DELTA( kernel[0], 0.008488347404, 1e-9 );
    TS_ASSERT_DELTA( kernel[1], 0.03807782601, 1e-8 );
    TS_ASSERT_DELTA( kernel[2], 0.1111650246, 1e-7 );
    TS_ASSERT_DELTA( kernel[3], 0.2113567063, 1e-7 );
    TS_ASSERT_DELTA( kernel[4], 0.2618241916, 1e-7 );
    TS_ASSERT_DELTA( kernel[5], 0.2113567063, 1e-7 );
    TS_ASSERT_DELTA( kernel[6], 0.1111650246, 1e-7 );
    TS_ASSERT_DELTA( kernel[7], 0.03807782601, 1e-8 );
    TS_ASSERT_DELTA( kernel[8], 0.008488347404, 1e-9 );
    generate_gaussian_kernel( kernel, 0 );
    TS_ASSERT_EQUALS( kernel.size(), 0 );
  }

  void test_generate_derivative_kernel() {
    std::vector<double> kernel;
    generate_derivative_kernel( kernel, 1, 5 );
    TS_ASSERT_EQUALS( kernel.size(), 5 );
    TS_ASSERT_DELTA( kernel[0], -0.08333333333, 1e-8 );
    TS_ASSERT_DELTA( kernel[1], 0.6666666667, 1e-7 );
    TS_ASSERT_DELTA( kernel[2], 0, 1e-7 );
    TS_ASSERT_DELTA( kernel[3], -0.6666666667, 1e-7 );
    TS_ASSERT_DELTA( kernel[4], 0.08333333333, 1e-8 );
    generate_derivative_kernel( kernel, 4, 5 );
    TS_ASSERT_EQUALS( kernel.size(), 5 );
    TS_ASSERT_DELTA( kernel[0], 1.0, 1e-7 );
    TS_ASSERT_DELTA( kernel[1], -4.0, 1e-7 );
    TS_ASSERT_DELTA( kernel[2], 6.0, 1e-7 );
    TS_ASSERT_DELTA( kernel[3], -4.0, 1e-7 );
    TS_ASSERT_DELTA( kernel[4], 1.0, 1e-7 );
    generate_derivative_kernel( kernel, 5 );
    TS_ASSERT_EQUALS( kernel.size(), 7 );
    TS_ASSERT_DELTA( kernel[0], 0.5, 1e-8 );
    TS_ASSERT_DELTA( kernel[1], -2.0, 1e-7 );
    TS_ASSERT_DELTA( kernel[2], 2.5, 1e-7 );
    TS_ASSERT_DELTA( kernel[3], 0, 1e-8 );
    TS_ASSERT_DELTA( kernel[4], -2.5, 1e-7 );
    TS_ASSERT_DELTA( kernel[5], 2.0, 1e-7 );
    TS_ASSERT_DELTA( kernel[6], -0.5, 1e-8 );
    generate_derivative_kernel( kernel, 1 );
    TS_ASSERT_EQUALS( kernel.size(), 3 );
    TS_ASSERT_EQUALS( kernel[0], 0.5 );
    TS_ASSERT_EQUALS( kernel[1], 0 );
    TS_ASSERT_EQUALS( kernel[2], -0.5 );
    generate_derivative_kernel( kernel, 2 );
    TS_ASSERT_EQUALS( kernel.size(), 3 );
    TS_ASSERT_EQUALS( kernel[0], 1 );
    TS_ASSERT_EQUALS( kernel[1], -2 );
    TS_ASSERT_EQUALS( kernel[2], 1 );
    generate_gaussian_kernel( kernel, 0 );
    TS_ASSERT_EQUALS( kernel.size(), 0 );
  }

  void test_derivative_filter() {
    ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
    ImageView<double> dst = derivative_filter( src, 1, 0 );
    TS_ASSERT_EQUALS( dst(0,0), 0.5 );
    TS_ASSERT_EQUALS( dst(1,0), 0.5 );
    TS_ASSERT_EQUALS( dst(0,1), 0.5 );
    TS_ASSERT_EQUALS( dst(1,1), 0.5 );
  }

  void test_gaussian_filter() {
    ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
    ImageView<double> dst = gaussian_filter( src, 1.0, 0, 5, 0, ZeroEdgeExtend() );
    TS_ASSERT_DELTA( dst(0,0), 0.3877403988*1+0.2447702197*2, 1e-7 );
    TS_ASSERT_DELTA( dst(1,0), 0.3877403988*2+0.2447702197*1, 1e-7 );
    TS_ASSERT_DELTA( dst(0,1), 0.3877403988*3+0.2447702197*4, 1e-7 );
    TS_ASSERT_DELTA( dst(1,1), 0.3877403988*4+0.2447702197*3, 1e-7 );
  }

  void test_laplacian_filter() {
    ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
    ImageView<double> dst = laplacian_filter( src, ZeroEdgeExtend() );
    TS_ASSERT_EQUALS( dst(0,0), 1 );
    TS_ASSERT_EQUALS( dst(1,0), -3 );
    TS_ASSERT_EQUALS( dst(0,1), -7 );
    TS_ASSERT_EQUALS( dst(1,1), -11 );
  }

  void test_per_pixel_filter() {
    ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
    ImageView<double> dst = per_pixel_filter( src, ArgValProductFunctor<double>(0.5) );
    TS_ASSERT_EQUALS( dst(0,0), 0.5 );
    TS_ASSERT_EQUALS( dst(1,0), 1.0 );
    TS_ASSERT_EQUALS( dst(0,1), 1.5 );
    TS_ASSERT_EQUALS( dst(1,1), 2.0 );
  }

  void test_per_pixel_channel_filter() {
    ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
    ImageView<double> dst = per_pixel_channel_filter( src, ArgValProductFunctor<double>(2.0) );
    TS_ASSERT_EQUALS( dst(0,0), 2 );
    TS_ASSERT_EQUALS( dst(1,0), 4 );
    TS_ASSERT_EQUALS( dst(0,1), 6 );
    TS_ASSERT_EQUALS( dst(1,1), 8 );
  }

  void test_convolution_filter() {
    ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
    ImageView<double> krn(2,2); krn(0,0)=2; krn(1,0)=-1; krn(0,1)=0; krn(1,1)=3;

    // Test rasterization
    ImageView<double> dst = convolution_filter( src, krn, ZeroEdgeExtend() );
    TS_ASSERT_EQUALS( dst.cols(), 2 );
    TS_ASSERT_EQUALS( dst.rows(), 2 );
    TS_ASSERT_EQUALS( dst(0,0), 2 );
    TS_ASSERT_EQUALS( dst(1,0), 3 );
    TS_ASSERT_EQUALS( dst(0,1), 6 );
    TS_ASSERT_EQUALS( dst(1,1), 8 );

    // Test the traits
    TS_ASSERT( is_of_type<double>( convolution_filter( src, krn )(0,0) ) );
    TS_ASSERT( is_of_type<double>( *(convolution_filter( src, krn ).origin()) ) );
  }

  void test_separable_convolution_filter() {
    ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
    vector<double> krn; krn.push_back(1); krn.push_back(-1);
    SeparableConvolutionView<ImageView<double>,double,ZeroEdgeExtend> cnv( src, krn, krn );

    // Test rasterization
    ImageView<double> dst = separable_convolution_filter( src, krn, krn, ZeroEdgeExtend() );
    TS_ASSERT_EQUALS( dst.cols(), 2 );
    TS_ASSERT_EQUALS( dst.rows(), 2 );
    TS_ASSERT_EQUALS( dst.planes(), 1 );
    TS_ASSERT_EQUALS( dst(0,0), 1 );
    TS_ASSERT_EQUALS( dst(0,1), 2 );
    TS_ASSERT_EQUALS( dst(1,0), 1 );
    TS_ASSERT_EQUALS( dst(1,1), 0 );

    // Test the traits
    TS_ASSERT( is_of_type<double>( separable_convolution_filter( src, krn, krn, ZeroEdgeExtend() )(0,0) ) );
    TS_ASSERT( is_of_type<double>( *(separable_convolution_filter( src, krn, krn, ZeroEdgeExtend() ).origin()) ) );
  }

  void test_2x2_conv_box_filter() {
    ImageView<double> kernel(2,2); kernel(0,0)=0.25; kernel(1,0)=0.25; kernel(0,1)=0.25; kernel(1,1)=0.25;
    {
    ImageView<double> src(2,2); src(0,0)=1; src(1,0)=0; src(0,1)=0; src(1,1)=0;
    ConvolutionView<ImageView<double>,ImageView<double>,ZeroEdgeExtend> cnv( src, kernel, 1, 1 );
    TS_ASSERT_EQUALS( cnv.cols(), 2 );
    TS_ASSERT_EQUALS( cnv.rows(), 2 );
    TS_ASSERT_EQUALS( cnv.planes(), 1 );
    TS_ASSERT_EQUALS( cnv(0,0), 0.25 );
    TS_ASSERT_EQUALS( cnv(0,1), 0 );
    TS_ASSERT_EQUALS( cnv(1,0), 0 );
    TS_ASSERT_EQUALS( cnv(1,1), 0 );
    ImageView<double> dst = cnv;
    TS_ASSERT_EQUALS( dst.cols(), 2 );
    TS_ASSERT_EQUALS( dst.rows(), 2 );
    TS_ASSERT_EQUALS( dst.planes(), 1 );
    TS_ASSERT_EQUALS( dst(0,0), 0.25 );
    TS_ASSERT_EQUALS( dst(0,1), 0 );
    TS_ASSERT_EQUALS( dst(1,0), 0 );
    TS_ASSERT_EQUALS( dst(1,1), 0 );
    }
    {
    ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
    ImageView<double> dst = convolution_filter( src, kernel, 1, 1, ZeroEdgeExtend() );
    TS_ASSERT_EQUALS( dst.cols(), 2 );
    TS_ASSERT_EQUALS( dst.rows(), 2 );
    TS_ASSERT_EQUALS( dst.planes(), 1 );
    TS_ASSERT_EQUALS( dst(0,0), 2.5 );
    TS_ASSERT_EQUALS( dst(0,1), 1.75 );
    TS_ASSERT_EQUALS( dst(1,0), 1.5 );
    TS_ASSERT_EQUALS( dst(1,1), 1 );
    }
  }

  void test_2x2_sepconv_box_filter() {
    std::vector<float> kernel(2); kernel[0]=0.5; kernel[1]=0.5;
    {
    ImageView<double> src(2,2); src(0,0)=1; src(1,0)=0; src(0,1)=0; src(1,1)=0;
    ImageView<double> dst = separable_convolution_filter( src, kernel, kernel, 1, 1, ZeroEdgeExtend() );
    TS_ASSERT_EQUALS( dst.cols(), 2 );
    TS_ASSERT_EQUALS( dst.rows(), 2 );
    TS_ASSERT_EQUALS( dst.planes(), 1 );
    TS_ASSERT_EQUALS( dst(0,0), 0.25 );
    TS_ASSERT_EQUALS( dst(0,1), 0 );
    TS_ASSERT_EQUALS( dst(1,0), 0 );
    TS_ASSERT_EQUALS( dst(1,1), 0 );
    }
    {
    ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
    ImageView<double> dst = separable_convolution_filter( src, kernel, kernel, 1, 1, ZeroEdgeExtend() );
    TS_ASSERT_EQUALS( dst.cols(), 2 );
    TS_ASSERT_EQUALS( dst.rows(), 2 );
    TS_ASSERT_EQUALS( dst.planes(), 1 );
    TS_ASSERT_EQUALS( dst(0,0), 2.5 );
    TS_ASSERT_EQUALS( dst(0,1), 1.75 );
    TS_ASSERT_EQUALS( dst(1,0), 1.5 );
    TS_ASSERT_EQUALS( dst(1,1), 1 );
    }
  }

}; // class TestFilter
