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


// TestInterpolation.h
#include <gtest/gtest_VW.h>

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

TEST( Interpolation, DISABLED_PassThrough ) {
  typedef InterpolationView<FloatingView<float>, BilinearInterpolation> bilinear_view;
  typedef InterpolationView<FloatingView<float>, BicubicInterpolation> bicubic_view;
  typedef InterpolationView<FloatingView<float>, NearestPixelInterpolation> near_view;
  // If underlining image is floating point accessible, then
  // InterpolationView should pass through all traits and access.
  EXPECT_TRUE( bool_trait<IsFloatingPointIndexable>(
          bilinear_view(FloatingView<float>(5,5),BilinearInterpolation())) );
  EXPECT_TRUE( bool_trait<IsMultiplyAccessible>(
          bilinear_view(FloatingView<float>(5,5),BilinearInterpolation())) );
  EXPECT_TRUE( bool_trait<IsFloatingPointIndexable>(
          bicubic_view(FloatingView<float>(5,5),BicubicInterpolation())) );
  EXPECT_TRUE( bool_trait<IsMultiplyAccessible>(
          bicubic_view(FloatingView<float>(5,5),BicubicInterpolation())) );
  EXPECT_TRUE( bool_trait<IsFloatingPointIndexable>(
          near_view(FloatingView<float>(5,5),NearestPixelInterpolation())) );
  EXPECT_TRUE( bool_trait<IsMultiplyAccessible>(
          near_view(FloatingView<float>(5,5),NearestPixelInterpolation())) );

  // Verify that InterpolationView really isn't doing anything to mess
  // this up. However doesn't work for Bilinear and Bicubic.
  EXPECT_EQ( 0.25, (near_view(FloatingView<float>(5,5),
                              NearestPixelInterpolation())(0.5,0.5)) );

  // Verify that prerasterization is not happening for interpolations
  // of already floating point accessible images.
  near_view::prerasterize_type pre =
    near_view(FloatingView<float>(5,5),
              NearestPixelInterpolation()).prerasterize(BBox2i(0,0,5,5));
  EXPECT_EQ( 0.25, pre(0.5,0.5) );
}
