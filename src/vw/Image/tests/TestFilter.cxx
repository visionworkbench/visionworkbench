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


// TestFilter.h
#include <gtest/gtest_VW.h>

#include <vw/Image/Filter.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>

#include <vector>

using namespace vw;

template <template<class> class TraitT, class T>
static bool bool_trait( T const& /*arg*/ ) {
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

TEST( Filter, GaussianKernel ) {
  std::vector<double> kernel;
  generate_gaussian_kernel( kernel, 1.0, 5 );
  EXPECT_EQ( int(kernel.size()), 5 );
  EXPECT_NEAR( kernel[0], 0.06135958087, 1e-8 );
  EXPECT_NEAR( kernel[1], 0.2447702197, 1e-7 );
  EXPECT_NEAR( kernel[2], 0.3877403988, 1e-7 );
  EXPECT_NEAR( kernel[3], 0.2447702197, 1e-7 );
  EXPECT_NEAR( kernel[4], 0.06135958087, 1e-8 );
  generate_gaussian_kernel( kernel, 1.0, 4 );
  EXPECT_EQ( int(kernel.size()), 4 );
  EXPECT_NEAR( kernel[0], 0.1423836140, 1e-7 );
  EXPECT_NEAR( kernel[1], 0.3576163860, 1e-7 );
  EXPECT_NEAR( kernel[2], 0.3576163860, 1e-7 );
  EXPECT_NEAR( kernel[3], 0.1423836140, 1e-7 );
  generate_gaussian_kernel( kernel, 1.5 );
  EXPECT_EQ( int(kernel.size()), 9 );
  EXPECT_NEAR( kernel[0], 0.008488347404, 1e-9 );
  EXPECT_NEAR( kernel[1], 0.03807782601, 1e-8 );
  EXPECT_NEAR( kernel[2], 0.1111650246, 1e-7 );
  EXPECT_NEAR( kernel[3], 0.2113567063, 1e-7 );
  EXPECT_NEAR( kernel[4], 0.2618241916, 1e-7 );
  EXPECT_NEAR( kernel[5], 0.2113567063, 1e-7 );
  EXPECT_NEAR( kernel[6], 0.1111650246, 1e-7 );
  EXPECT_NEAR( kernel[7], 0.03807782601, 1e-8 );
  EXPECT_NEAR( kernel[8], 0.008488347404, 1e-9 );
  generate_gaussian_kernel( kernel, 0 );
  EXPECT_EQ( int(kernel.size()), 0 );
}

TEST( Filter, DerivativeKernel ) {
  std::vector<double> kernel;
  generate_derivative_kernel( kernel, 1, 5 );
  EXPECT_EQ( int(kernel.size()), 5 );
  EXPECT_NEAR( kernel[0], -0.08333333333, 1e-8 );
  EXPECT_NEAR( kernel[1], 0.6666666667, 1e-7 );
  EXPECT_NEAR( kernel[2], 0, 1e-7 );
  EXPECT_NEAR( kernel[3], -0.6666666667, 1e-7 );
  EXPECT_NEAR( kernel[4], 0.08333333333, 1e-8 );
  generate_derivative_kernel( kernel, 4, 5 );
  EXPECT_EQ( int(kernel.size()), 5 );
  EXPECT_NEAR( kernel[0], 1.0, 1e-7 );
  EXPECT_NEAR( kernel[1], -4.0, 1e-7 );
  EXPECT_NEAR( kernel[2], 6.0, 1e-7 );
  EXPECT_NEAR( kernel[3], -4.0, 1e-7 );
  EXPECT_NEAR( kernel[4], 1.0, 1e-7 );
  generate_derivative_kernel( kernel, 5 );
  EXPECT_EQ( int(kernel.size()), 7 );
  EXPECT_NEAR( kernel[0], 0.5, 1e-8 );
  EXPECT_NEAR( kernel[1], -2.0, 1e-7 );
  EXPECT_NEAR( kernel[2], 2.5, 1e-7 );
  EXPECT_NEAR( kernel[3], 0, 1e-8 );
  EXPECT_NEAR( kernel[4], -2.5, 1e-7 );
  EXPECT_NEAR( kernel[5], 2.0, 1e-7 );
  EXPECT_NEAR( kernel[6], -0.5, 1e-8 );
  generate_derivative_kernel( kernel, 1 );
  EXPECT_EQ( int(kernel.size()), 3 );
  EXPECT_EQ( kernel[0], 0.5 );
  EXPECT_EQ( kernel[1], 0 );
  EXPECT_EQ( kernel[2], -0.5 );
  generate_derivative_kernel( kernel, 2 );
  EXPECT_EQ( int(kernel.size()), 3 );
  EXPECT_EQ( kernel[0], 1 );
  EXPECT_EQ( kernel[1], -2 );
  EXPECT_EQ( kernel[2], 1 );
  generate_gaussian_kernel( kernel, 0 );
  EXPECT_EQ( int(kernel.size()), 0 );
}

TEST( Filter, Derivative ) {
  ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
  ImageView<double> dst = derivative_filter( src, 1, 0 );
  EXPECT_EQ( dst(0,0), 0.5 );
  EXPECT_EQ( dst(1,0), 0.5 );
  EXPECT_EQ( dst(0,1), 0.5 );
  EXPECT_EQ( dst(1,1), 0.5 );
}

TEST( Filter, Gaussian ) {
  ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
  ImageView<double> dst = gaussian_filter( src, 1.0, 0, 5, 0, ZeroEdgeExtension() );
  EXPECT_NEAR( dst(0,0), 0.3877403988*1+0.2447702197*2, 1e-7 );
  EXPECT_NEAR( dst(1,0), 0.3877403988*2+0.2447702197*1, 1e-7 );
  EXPECT_NEAR( dst(0,1), 0.3877403988*3+0.2447702197*4, 1e-7 );
  EXPECT_NEAR( dst(1,1), 0.3877403988*4+0.2447702197*3, 1e-7 );
}

TEST( Filter, Laplacian ) {
  ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
  ImageView<double> dst = laplacian_filter( src, ZeroEdgeExtension() );
  EXPECT_EQ( dst(0,0), 1 );
  EXPECT_EQ( dst(1,0), -3 );
  EXPECT_EQ( dst(0,1), -7 );
  EXPECT_EQ( dst(1,1), -11 );
  ImageView<PixelGrayA<double> > src2(2,2);
  ImageView<PixelGrayA<double> > dst2a = laplacian_filter( src, ZeroEdgeExtension() );
  ImageView<PixelGrayA<double> > dst2b = laplacian_filter( src2, ZeroEdgeExtension() );
}

TEST( Filter, Sobel ) {

  int n = 5;
  ImageView<double> fun(n, n), dx(n, n), dy(n, n);
  for (int x = 0; x < n; x++){
    for (int y = 0; y < n; y++){
      fun(x, y) = x*x;
    }
  }

  bool do_x_deriv = true;
  ImageView<double> sobel_x = sobel_filter( fun, do_x_deriv );

  do_x_deriv = false;
  ImageView<double> sobel_y = sobel_filter( fun, do_x_deriv );

  // The Sobel derivative in x is negative (we'd expect d(x^2) = 2*x > 0) due to
  // how the Sobel coefficients are chosen in the literature.
  EXPECT_EQ(sobel_x(n/2, n/2), -32);

  // The Sobel derivative in y is 0, since the input function is constant in y.
  EXPECT_EQ(sobel_y(n/2, n/2), 0);
}

TEST( Filter, PerPixel ) {
  ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
  ImageView<double> dst = per_pixel_filter( src, ArgValProductFunctor<double>(0.5) );
  EXPECT_EQ( dst(0,0), 0.5 );
  EXPECT_EQ( dst(1,0), 1.0 );
  EXPECT_EQ( dst(0,1), 1.5 );
  EXPECT_EQ( dst(1,1), 2.0 );
}

TEST( Filter, PerChannel ) {
  ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
  ImageView<double> dst = per_pixel_channel_filter( src, ArgValProductFunctor<double>(2.0) );
  EXPECT_EQ( dst(0,0), 2 );
  EXPECT_EQ( dst(1,0), 4 );
  EXPECT_EQ( dst(0,1), 6 );
  EXPECT_EQ( dst(1,1), 8 );
}

TEST( Filter, Convolution ) {
  ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
  ImageView<double> krn(2,2); krn(0,0)=2; krn(1,0)=-1; krn(0,1)=0; krn(1,1)=3;

  // Test rasterization
  ImageView<double> dst = convolution_filter( src, krn, ZeroEdgeExtension() );
  EXPECT_EQ( dst.cols(), 2 );
  EXPECT_EQ( dst.rows(), 2 );
  EXPECT_EQ( dst(0,0), 2 );
  EXPECT_EQ( dst(1,0), 3 );
  EXPECT_EQ( dst(0,1), 6 );
  EXPECT_EQ( dst(1,1), 8 );

  // Test the traits
  ASSERT_TRUE( is_of_type<double>( convolution_filter( src, krn )(0,0) ) );
  ASSERT_TRUE( is_of_type<double>( *(convolution_filter( src, krn ).origin()) ) );
}

TEST( Filter, SepConvolution ) {
  ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
  std::vector<double> krn; krn.push_back(1); krn.push_back(-1);
  SeparableConvolutionView<ImageView<double>,double,ZeroEdgeExtension> cnv( src, krn, krn );

  // Test rasterization
  ImageView<double> dst = separable_convolution_filter( src, krn, krn, ZeroEdgeExtension() );
  EXPECT_EQ( dst.cols(), 2 );
  EXPECT_EQ( dst.rows(), 2 );
  EXPECT_EQ( dst.planes(), 1 );
  EXPECT_EQ( dst(0,0), 1 );
  EXPECT_EQ( dst(0,1), 2 );
  EXPECT_EQ( dst(1,0), 1 );
  EXPECT_EQ( dst(1,1), 0 );

  // Test the traits
  ASSERT_TRUE( is_of_type<double>( separable_convolution_filter( src, krn, krn, ZeroEdgeExtension() )(0,0) ) );
  ASSERT_TRUE( is_of_type<double>( *(separable_convolution_filter( src, krn, krn, ZeroEdgeExtension() ).origin()) ) );
}

TEST( Filter, Conv22Box ) {
  ImageView<double> kernel(2,2); kernel(0,0)=0.25; kernel(1,0)=0.25; kernel(0,1)=0.25; kernel(1,1)=0.25;
  {
    ImageView<double> src(2,2); src(0,0)=1; src(1,0)=0; src(0,1)=0; src(1,1)=0;
    ConvolutionView<ImageView<double>,ImageView<double>,ZeroEdgeExtension> cnv( src, kernel, 1, 1 );
    EXPECT_EQ( cnv.cols(), 2 );
    EXPECT_EQ( cnv.rows(), 2 );
    EXPECT_EQ( cnv.planes(), 1 );
    EXPECT_EQ( cnv(0,0), 0.25 );
    EXPECT_EQ( cnv(0,1), 0 );
    EXPECT_EQ( cnv(1,0), 0 );
    EXPECT_EQ( cnv(1,1), 0 );
    ImageView<double> dst = cnv;
    EXPECT_EQ( dst.cols(), 2 );
    EXPECT_EQ( dst.rows(), 2 );
    EXPECT_EQ( dst.planes(), 1 );
    EXPECT_EQ( dst(0,0), 0.25 );
    EXPECT_EQ( dst(0,1), 0 );
    EXPECT_EQ( dst(1,0), 0 );
    EXPECT_EQ( dst(1,1), 0 );
  }
  {
    ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
    ImageView<double> dst = convolution_filter( src, kernel, 1, 1, ZeroEdgeExtension() );
    EXPECT_EQ( dst.cols(), 2 );
    EXPECT_EQ( dst.rows(), 2 );
    EXPECT_EQ( dst.planes(), 1 );
    EXPECT_EQ( dst(0,0), 2.5 );
    EXPECT_EQ( dst(0,1), 1.75 );
    EXPECT_EQ( dst(1,0), 1.5 );
    EXPECT_EQ( dst(1,1), 1 );
  }
}

TEST( Filter, SepConv22Box ) {
  std::vector<float> kernel(2); kernel[0]=0.5; kernel[1]=0.5;
  {
    ImageView<double> src(2,2); src(0,0)=1; src(1,0)=0; src(0,1)=0; src(1,1)=0;
    ImageView<double> dst = separable_convolution_filter( src, kernel, kernel, 1, 1, ZeroEdgeExtension() );
    EXPECT_EQ( dst.cols(), 2 );
    EXPECT_EQ( dst.rows(), 2 );
    EXPECT_EQ( dst.planes(), 1 );
    EXPECT_EQ( dst(0,0), 0.25 );
    EXPECT_EQ( dst(0,1), 0 );
    EXPECT_EQ( dst(1,0), 0 );
    EXPECT_EQ( dst(1,1), 0 );
  }
  {
    ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
    ImageView<double> dst = separable_convolution_filter( src, kernel, kernel, 1, 1, ZeroEdgeExtension() );
    EXPECT_EQ( dst.cols(), 2 );
    EXPECT_EQ( dst.rows(), 2 );
    EXPECT_EQ( dst.planes(), 1 );
    EXPECT_EQ( dst(0,0), 2.5 );
    EXPECT_EQ( dst(0,1), 1.75 );
    EXPECT_EQ( dst(1,0), 1.5 );
    EXPECT_EQ( dst(1,1), 1 );
  }
}
