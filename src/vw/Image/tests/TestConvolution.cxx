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


// TestConvolution.h
#include <gtest/gtest_VW.h>

#include <vw/Core/TypeDeduction.h>
#include <vw/Image/AlgorithmFunctions.h>
#include <vw/Image/Convolution.h>
#include <vw/Image/Filter.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/PixelTypes.h>

#include <test/Helpers.h>

using namespace vw;

// Simple Comparison
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

TEST( Convolution, Point1D ) {
  ImageView<double> src(3,1); src(0,0)=1; src(1,0)=2; src(2,0)=3;
  double kernel[] = {2,3,-1};
  EXPECT_EQ( correlate_1d_at_point(src.origin(), (double*)kernel,3), 5 );
  ASSERT_TRUE( is_of_type<double>( correlate_1d_at_point( ImageView<double>().origin(), (double*)0, 0 ) ) );
  ASSERT_TRUE( is_of_type<PixelRGB<float> >( correlate_1d_at_point( ImageView<PixelRGB<uint8> >().origin(), (float*)0, 0 ) ) );
}

TEST( Convolution, Point2D ) {
  ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
  ImageView<double> krn(2,2); krn(0,0)=2; krn(1,0)=-1; krn(0,1)=0; krn(1,1)=3;
  EXPECT_EQ( correlate_2d_at_point(src.origin(),krn.origin(),2,2), 12 );
  ASSERT_TRUE( is_of_type<double>( correlate_2d_at_point( ImageView<double>().origin(), ImageView<double>().origin(), 0, 0 ) ) );
  ASSERT_TRUE( is_of_type<PixelRGB<float> >( correlate_2d_at_point( ImageView<PixelRGB<uint8> >().origin(), ImageView<float>().origin(), 0, 0 ) ) );
}

TEST( Convolution, View ) {
  ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
  ImageView<double> krn(2,2); krn(0,0)=2; krn(1,0)=-1; krn(0,1)=0; krn(1,1)=3;
  ConvolutionView<ImageView<double>,ImageView<double>,ZeroEdgeExtension> cnv( src, krn );

  // Test individual pixel access
  EXPECT_EQ( cnv.cols(), 2 );
  EXPECT_EQ( cnv.rows(), 2 );
  EXPECT_EQ( cnv.planes(), 1 );
  EXPECT_EQ( cnv(0,0), 2 );
  EXPECT_EQ( cnv(1,0), 3 );
  EXPECT_EQ( cnv(0,1), 6 );
  EXPECT_EQ( cnv(1,1), 8 );

  // Test rasterization
  ImageView<double> dst = cnv;
  EXPECT_EQ( dst.cols(), 2 );
  EXPECT_EQ( dst.rows(), 2 );
  EXPECT_EQ( dst(0,0), 2 );
  EXPECT_EQ( dst(1,0), 3 );
  EXPECT_EQ( dst(0,1), 6 );
  EXPECT_EQ( dst(1,1), 8 );

  // Test the traits
  ASSERT_TRUE( is_of_type<double>( cnv(0,0) ) );
  ASSERT_TRUE( is_of_type<double>( *(cnv.origin()) ) );
  ASSERT_FALSE( bool_trait<IsResizable>( cnv ) );
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( cnv ) );
  ASSERT_FALSE( bool_trait<IsFloatingPointIndexable>( cnv ) );
  ASSERT_TRUE( bool_trait<IsImageView>( cnv ) );
}

TEST( Convolution, SeparableView ) {
  ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
  std::vector<double> krn; krn.push_back(1); krn.push_back(-1);
  SeparableConvolutionView<ImageView<double>,double,ZeroEdgeExtension> cnv( src, krn, krn );

  // Test individual pixel access
  EXPECT_EQ( cnv.cols(), 2 );
  EXPECT_EQ( cnv.rows(), 2 );
  EXPECT_EQ( cnv.planes(), 1 );
  EXPECT_EQ( cnv(0,0), 1 );
  EXPECT_EQ( cnv(0,1), 2 );
  EXPECT_EQ( cnv(1,0), 1 );
  EXPECT_EQ( cnv(1,1), 0 );

  // Test rasterization
  ImageView<double> dst = cnv;
  EXPECT_EQ( dst.cols(), 2 );
  EXPECT_EQ( dst.rows(), 2 );
  EXPECT_EQ( dst.planes(), 1 );
  EXPECT_EQ( dst(0,0), 1 );
  EXPECT_EQ( dst(0,1), 2 );
  EXPECT_EQ( dst(1,0), 1 );
  EXPECT_EQ( dst(1,1), 0 );

  // Test the traits
  ASSERT_TRUE( is_of_type<double>( cnv(0,0) ) );
  ASSERT_TRUE( is_of_type<double>( *(cnv.origin()) ) );
  ASSERT_FALSE( bool_trait<IsResizable>( cnv ) );
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( cnv ) );
  ASSERT_FALSE( bool_trait<IsFloatingPointIndexable>( cnv ) );
  ASSERT_TRUE( bool_trait<IsImageView>( cnv ) );
}

TEST( Convolution, SeparableView_0x2 ) {
  ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=4; src(1,1)=6;
  std::vector<double> krnx;
  std::vector<double> krny; krny.push_back(1); krny.push_back(-1);
  SeparableConvolutionView<ImageView<double>,double,ZeroEdgeExtension> cnv( src, krnx, krny );

  // Test individual pixel access
  EXPECT_EQ( cnv.cols(), 2 );
  EXPECT_EQ( cnv.rows(), 2 );
  EXPECT_EQ( cnv.planes(), 1 );
  EXPECT_EQ( cnv(0,0), 1 );
  EXPECT_EQ( cnv(0,1), 3 );
  EXPECT_EQ( cnv(1,0), 2 );
  EXPECT_EQ( cnv(1,1), 4 );

  // Test rasterization
  ImageView<double> dst = cnv;
  EXPECT_EQ( dst.cols(), 2 );
  EXPECT_EQ( dst.rows(), 2 );
  EXPECT_EQ( dst.planes(), 1 );
  EXPECT_EQ( dst(0,0), 1 );
  EXPECT_EQ( dst(0,1), 3 );
  EXPECT_EQ( dst(1,0), 2 );
  EXPECT_EQ( dst(1,1), 4 );

  // Test the traits
  ASSERT_TRUE( is_of_type<double>( cnv(0,0) ) );
  ASSERT_TRUE( is_of_type<double>( *(cnv.origin()) ) );
  ASSERT_FALSE( bool_trait<IsResizable>( cnv ) );
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( cnv ) );
  ASSERT_FALSE( bool_trait<IsFloatingPointIndexable>( cnv ) );
  ASSERT_TRUE( bool_trait<IsImageView>( cnv ) );
}

TEST( Convolution, SeparableView_2x0 ) {
  ImageView<double> src(2,2); src(0,0)=1; src(1,0)=2; src(0,1)=3; src(1,1)=4;
  std::vector<double> krnx; krnx.push_back(1); krnx.push_back(-1);
  std::vector<double> krny;
  SeparableConvolutionView<ImageView<double>,double,ZeroEdgeExtension> cnv( src, krnx, krny );

  // Test individual pixel access
  EXPECT_EQ( cnv.cols(), 2 );
  EXPECT_EQ( cnv.rows(), 2 );
  EXPECT_EQ( cnv.planes(), 1 );
  EXPECT_EQ( cnv(0,0), 1 );
  EXPECT_EQ( cnv(1,0), 1 );
  EXPECT_EQ( cnv(0,1), 3 );
  EXPECT_EQ( cnv(1,1), 1 );

  // Test rasterization
  ImageView<double> dst = cnv;
  EXPECT_EQ( dst.cols(), 2 );
  EXPECT_EQ( dst.rows(), 2 );
  EXPECT_EQ( dst.planes(), 1 );
  EXPECT_EQ( dst(0,0), 1 );
  EXPECT_EQ( dst(1,0), 1 );
  EXPECT_EQ( dst(0,1), 3 );
  EXPECT_EQ( dst(1,1), 1 );

  // Test the traits
  ASSERT_TRUE( is_of_type<double>( cnv(0,0) ) );
  ASSERT_TRUE( is_of_type<double>( *(cnv.origin()) ) );
  ASSERT_FALSE( bool_trait<IsResizable>( cnv ) );
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( cnv ) );
  ASSERT_FALSE( bool_trait<IsFloatingPointIndexable>( cnv ) );
  ASSERT_TRUE( bool_trait<IsImageView>( cnv ) );
}

TEST( Convolution, SeparableView_Compound ) {
  ImageView<PixelGray<float32> > src(2,2); src(0,0)=1; src(1,0)=0.2; src(0,1)=0.3; src(1,1)=0.4;
  std::vector<double> krn; krn.push_back(1); krn.push_back(-1);
  SeparableConvolutionView<ImageView<PixelGray<float32> >,double,ZeroEdgeExtension> cnv( src, krn, krn );
  ASSERT_TRUE( is_of_type<PixelGray<float32> >( cnv(0,0) ) );
}

// This unit test catches a bug in the separable convolution code
// that was causing the shift due to a crop operation to be applied
// twice to image view operations that included two layers of edge
// extension or convolution.
TEST( Convolution, Prerasterize ) {
  ImageView<float> test_image(1024,1024);
  fill(test_image,1.0);

  ImageViewRef<float> two_edge_extend_operations = edge_extend(gaussian_filter(channel_cast<float>(channels_to_planes(test_image)),1.5), ZeroEdgeExtension());

  BBox2i bbox(100,0,1024,1024);
  ImageView<float> right_buf = crop( two_edge_extend_operations, bbox );

  EXPECT_EQ(right_buf(1000,100), 0.0);
  EXPECT_EQ(right_buf(900,100), 1.0);
}
