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


// TestImageMath.h
#include <gtest/gtest_VW.h>

#include <vw/config.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/AlgorithmFunctions.h>
#include <test/Helpers.h>

using namespace vw;

template <typename T>
class MaskedViewTest : public ::testing::Test {
public:
  typedef T Px;
  typedef typename PixelChannelType<Px>::type ChT;
  typedef PixelMask<T> MPx;

  MaskedViewTest() {}

  virtual void SetUp() {
    a.set_size(2,2);
    b.set_size(2,2);

    a(0,0) = MPx(1); a(1,0) = MPx(2);
    a(0,1) = MPx(3); a(1,1) = MPx(4);
    b(0,0) = MPx(10); b(1,0) = MPx(5);
    b(0,1) = MPx(7); b(1,1) = MPx(9);
    a(0,1).invalidate();
    a(1,1).invalidate();
    b(1,0).invalidate();
    b(1,1).invalidate();
  }

  ImageView<MPx> a, b;
};

typedef MaskedViewTest<uint8> MaskedViewU8Test;
typedef MaskedViewTest<PixelGray<float> > MaskedViewGrayTest;

typedef ::testing::Types<uint8, int16, uint16, float, double, PixelGray<uint8>, PixelGray<int16>, PixelGray<uint16>, PixelGray<float>, PixelRGB<uint8>, PixelRGB<int16>, PixelRGB<uint16>, PixelRGB<float> > MyTypes;
TYPED_TEST_CASE( MaskedViewTest, MyTypes );

TYPED_TEST( MaskedViewTest, create_mask ) {
  ImageView<typename TestFixture::Px> non(2,1);
  non(0,0) = 10;
  non(1,0) = 0;
  ImageView<typename TestFixture::MPx> c = create_mask(non);
  EXPECT_EQ( typename TestFixture::Px(10), c(0,0).child() );
  EXPECT_EQ( typename TestFixture::Px(0), c(1,0).child() );
  EXPECT_TRUE( is_valid(c(0,0)) );
  EXPECT_FALSE( is_valid(c(1,0)) );

  non(1,0) = 20;
  c = create_mask(non,20);
  EXPECT_EQ( typename TestFixture::Px(10), c(0,0).child() );
  EXPECT_EQ( typename TestFixture::Px(0), c(1,0).child() );
  EXPECT_TRUE( is_valid(c(0,0)) );
  EXPECT_FALSE( is_valid(c(1,0)) );
}

TYPED_TEST( MaskedViewTest, apply_mask ) {
  ImageView<typename TestFixture::Px> c = apply_mask(this->a);
  EXPECT_EQ( 2, c.cols() ); EXPECT_EQ( 2, c.rows() );
  EXPECT_EQ( typename TestFixture::Px(1), c(0,0) );
  EXPECT_EQ( typename TestFixture::Px(2), c(1,0) );
  EXPECT_EQ( typename TestFixture::Px(0), c(0,1) );
  EXPECT_EQ( typename TestFixture::Px(0), c(1,1) );
  ImageView<typename TestFixture::Px> d = apply_mask(this->a,200);
  EXPECT_EQ( typename TestFixture::Px(1), d(0,0) );
  EXPECT_EQ( typename TestFixture::Px(2), d(1,0) );
  EXPECT_EQ( typename TestFixture::Px(200), d(0,1) );
  EXPECT_EQ( typename TestFixture::Px(200), d(1,1) );
}

TYPED_TEST( MaskedViewTest, copy_mask ) {
  ImageView<typename TestFixture::Px> in(2,2);
  in(0,0) = typename TestFixture::Px(9); in(1,0) = typename TestFixture::Px(8);
  in(0,1) = typename TestFixture::Px(1); in(1,1) = typename TestFixture::Px(4);
  ImageView<typename TestFixture::MPx> c = copy_mask(in,this->a);
  EXPECT_EQ( 2, c.cols() ); EXPECT_EQ( 2, c.rows() );
  EXPECT_TRUE( is_valid(c(0,0)) );
  EXPECT_TRUE( is_valid(c(1,0)) );
  EXPECT_FALSE(is_valid(c(0,1)) );
  EXPECT_FALSE(is_valid(c(1,1)) );
  EXPECT_EQ( typename TestFixture::Px(9), c(0,0).child() );
  EXPECT_EQ( typename TestFixture::Px(8), c(1,0).child() );
  EXPECT_EQ( typename TestFixture::Px(1), c(0,1).child() );
  EXPECT_EQ( typename TestFixture::Px(4), c(1,1).child() );

  // Test a copy mask where the input is a different channel type
  ImageView< PixelRGB<uint8> > in_rgbu8(2,2);
  in_rgbu8(0,0) = PixelRGB<uint8>(9); in_rgbu8(1,0) = PixelRGB<uint8>(8);
  in_rgbu8(0,1) = PixelRGB<uint8>(1); in_rgbu8(1,1) = PixelRGB<uint8>(4);
  ImageView<PixelMask<PixelRGB<uint8> > > out_rgbu8 = copy_mask( in_rgbu8, this->a );
  EXPECT_EQ( 2, out_rgbu8.cols() ); EXPECT_EQ( 2, out_rgbu8.rows() );
  EXPECT_TRUE( is_valid(out_rgbu8(0,0)) );
  EXPECT_TRUE( is_valid(out_rgbu8(1,0)) );
  EXPECT_FALSE(is_valid(out_rgbu8(0,1)) );
  EXPECT_FALSE(is_valid(out_rgbu8(1,1)) );

  // Test a copy_mask where the mask is a different type than the
  // input.
  ImageView< PixelMask<uint8> > mask(2,2);
  mask(0,0) = PixelMask<uint8>(); mask(1,0) = PixelMask<uint8>(2);
  mask(0,1) = PixelMask<uint8>(255); mask(1,1) = PixelMask<uint8>();
  c = copy_mask(in, mask);
  EXPECT_EQ( 2, c.cols() ); EXPECT_EQ( 2, c.rows() );
  EXPECT_FALSE( is_valid(c(0,0)) );
  EXPECT_TRUE( is_valid(c(1,0)) );
  EXPECT_TRUE( is_valid(c(0,1)) );
  EXPECT_FALSE(is_valid(c(1,1)) );
  EXPECT_EQ( typename TestFixture::Px(9), c(0,0).child() );
  EXPECT_EQ( typename TestFixture::Px(8), c(1,0).child() );
  EXPECT_EQ( typename TestFixture::Px(1), c(0,1).child() );
  EXPECT_EQ( typename TestFixture::Px(4), c(1,1).child() );
}

TEST_F( MaskedViewGrayTest, mask_to_alpha ) {
  typedef PixelWithAlpha<Px>::type APx;
  ImageView<APx> c = mask_to_alpha(a);
  EXPECT_EQ( 2, c.cols() ); EXPECT_EQ( 2, c.rows() );
  EXPECT_PIXEL_EQ( APx(1,1), c(0,0) );
  EXPECT_PIXEL_EQ( APx(2,1), c(1,0) );
  EXPECT_PIXEL_EQ( APx(0,0), c(0,1) );
  EXPECT_PIXEL_EQ( APx(0,0), c(1,1) );
}

TEST_F( MaskedViewGrayTest, alpha_to_mask ) {
  typedef PixelWithAlpha<Px>::type APx;
  ImageView<APx> c(2,2);
  c(0,0) = APx(0,0); c(1,0) = APx(2,1);
  c(0,1) = APx(3,0.25); c(1,1) = APx(4,0.75);
  a = alpha_to_mask(c);
  EXPECT_EQ( 0, a(0,0) );
  EXPECT_EQ( 2, a(1,0) );
  EXPECT_EQ( 3, a(0,1) );
  EXPECT_EQ( 4, a(1,1) );
  EXPECT_FALSE( is_valid(a(0,0)) );
  EXPECT_TRUE( is_valid(a(1,0)) );
  EXPECT_TRUE( is_valid(a(0,1)) );
  EXPECT_TRUE( is_valid(a(1,1)) );
}

TYPED_TEST( MaskedViewTest, edge_mask ) {
  ImageView<typename TestFixture::Px> c(4,4);
  fill( c, typename TestFixture::Px(1) );
  c(0,0) = typename TestFixture::Px(); c(0,1) = typename TestFixture::Px();
  c(0,2) = typename TestFixture::Px(); c(2,1) = typename TestFixture::Px();
  ImageView<typename TestFixture::MPx> d = edge_mask( c );
  EXPECT_FALSE( is_valid(d(0,0)) );
  EXPECT_FALSE( is_valid(d(0,2)) );
  EXPECT_TRUE( is_valid(d(2,1)) );
  EXPECT_TRUE( is_valid(d(2,2)) );
  EXPECT_FALSE( is_valid(d(3,3)) );
  EXPECT_FALSE( is_valid(d(0,3)) );
  EXPECT_FALSE( is_valid(d(3,0)) );
}

TYPED_TEST( MaskedViewTest, validate_mask ) {
  ImageView<typename TestFixture::MPx> c = validate_mask(this->a);
  EXPECT_EQ( 2, c.cols() ); EXPECT_EQ( 2, c.rows() );
  EXPECT_TRUE( is_valid(c(0,0)) );
  EXPECT_TRUE( is_valid(c(1,0)) );
  EXPECT_TRUE( is_valid(c(0,1)) );
  EXPECT_TRUE( is_valid(c(1,0)) );
  EXPECT_EQ( typename TestFixture::Px(1), c(0,0).child() );
  EXPECT_EQ( typename TestFixture::Px(4), c(1,1).child() );
}

TYPED_TEST( MaskedViewTest, invalidate_mask ) {
  ImageView<typename TestFixture::MPx> c = invalidate_mask(this->a);
  EXPECT_EQ( 2, c.cols() ); EXPECT_EQ( 2, c.rows() );
  EXPECT_FALSE( is_valid(c(0,0)) );
  EXPECT_FALSE( is_valid(c(1,0)) );
  EXPECT_FALSE( is_valid(c(0,1)) );
  EXPECT_FALSE( is_valid(c(1,0)) );
  EXPECT_EQ( typename TestFixture::Px(1), c(0,0).child() );
  EXPECT_EQ( typename TestFixture::Px(4), c(1,1).child() );
}

TYPED_TEST( MaskedViewTest, union_mask ) {
  ImageView<typename TestFixture::MPx> c = union_mask(this->a,this->b);
  EXPECT_EQ( 2, c.cols() ); EXPECT_EQ( 2, c.rows() );
  EXPECT_FALSE( is_valid(c(1,1)) );
  EXPECT_TRUE(  is_valid(c(0,0)) );
  EXPECT_TRUE(  is_valid(c(0,1)) );
  EXPECT_TRUE(  is_valid(c(1,0)) );
  EXPECT_EQ( typename TestFixture::Px(1), c(0,0).child() );
  EXPECT_EQ( typename TestFixture::Px(2), c(1,0).child() );
  EXPECT_EQ( typename TestFixture::Px(3), c(0,1).child() );
  EXPECT_EQ( typename TestFixture::Px(4), c(1,1).child() );
}

TYPED_TEST( MaskedViewTest, intersect_mask ) {
  ImageView<typename TestFixture::MPx> c =
    intersect_mask(this->a,this->b);
  EXPECT_EQ( 2, c.cols() ); EXPECT_EQ( 2, c.rows() );
  EXPECT_FALSE( is_valid(c(1,1)) );
  EXPECT_FALSE( is_valid(c(1,0)) );
  EXPECT_FALSE( is_valid(c(0,1)) );
  EXPECT_TRUE(  is_valid(c(0,0)) );
  EXPECT_EQ( typename TestFixture::Px(1), c(0,0).child() );
  EXPECT_EQ( typename TestFixture::Px(2), c(1,0).child() );
  EXPECT_EQ( typename TestFixture::Px(3), c(0,1).child() );
  EXPECT_EQ( typename TestFixture::Px(4), c(1,1).child() );
}

TYPED_TEST( MaskedViewTest, invert_mask ) {
  this->b = invert_mask(this->a);
  EXPECT_EQ( 2, this->b.cols() ); EXPECT_EQ( 2, this->b.rows() );
  EXPECT_NE( is_valid(this->a(0,0)), is_valid(this->b(0,0)) );
  EXPECT_NE( is_valid(this->a(1,0)), is_valid(this->b(1,0)) );
  EXPECT_NE( is_valid(this->a(0,1)), is_valid(this->b(0,1)) );
  EXPECT_NE( is_valid(this->a(1,1)), is_valid(this->b(1,1)) );
}

TYPED_TEST( MaskedViewTest, create_apply_mask ) {
  this->a(0,1).validate();
  this->a(1,1).validate();

  EXPECT_EQ( typename TestFixture::Px(1), this->a(0,0).child() );
  EXPECT_EQ( typename TestFixture::Px(2), this->a(1,0).child() );
  EXPECT_EQ( typename TestFixture::Px(3), this->a(0,1).child() );
  EXPECT_EQ( typename TestFixture::Px(4), this->a(1,1).child() );

  ImageView<typename TestFixture::Px> c =
    apply_mask(create_mask(apply_mask(this->a),3));
  EXPECT_EQ( typename TestFixture::Px(1), c(0,0) );
  EXPECT_EQ( typename TestFixture::Px(2), c(1,0) );
  EXPECT_EQ( typename TestFixture::Px(0), c(0,1) );
  EXPECT_EQ( typename TestFixture::Px(4), c(1,1) );

  ImageView<typename TestFixture::MPx> d =
    create_mask(apply_mask(this->a),3);
  EXPECT_EQ( typename TestFixture::Px(1), d(0,0).child() );
  EXPECT_EQ( typename TestFixture::Px(2), d(1,0).child() );
  EXPECT_EQ( typename TestFixture::Px(), d(0,1).child() );
  EXPECT_EQ( typename TestFixture::Px(4), d(1,1).child() );

  ImageView<typename TestFixture::Px> e =
    apply_mask(this->a);
  EXPECT_EQ( typename TestFixture::Px(1), e(0,0) );
  EXPECT_EQ( typename TestFixture::Px(2), e(1,0) );
  EXPECT_EQ( typename TestFixture::Px(3), e(0,1) );
  EXPECT_EQ( typename TestFixture::Px(4), e(1,1) );
}
