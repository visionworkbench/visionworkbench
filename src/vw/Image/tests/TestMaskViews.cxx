// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

// TestImageMath.h
#include <gtest/gtest.h>

#include <vw/config.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/Algorithms.h>
#include <test/Helpers.h>

using namespace vw;

template <class PixelT>
class MaskedViewTest : public ::testing::Test {
protected:
  typedef PixelT Px;
  typedef typename PixelChannelType<Px>::type ChT;
  typedef PixelMask<PixelT> MPx;

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

TEST_F( MaskedViewU8Test, create_mask ) {
  ImageView<Px> non(2,1);
  non(0,0) = 10;
  non(1,0) = 0;
  ImageView<MPx> c = create_mask(non);
  EXPECT_EQ( 10, c(0,0) );
  EXPECT_EQ( 0, c(1,0) );
  EXPECT_TRUE( is_valid(c(0,0)) );
  EXPECT_FALSE( is_valid(c(1,0)) );

  non(1,0) = 20;
  c = create_mask(non,20);
  EXPECT_EQ( 10, c(0,0) );
  EXPECT_EQ( 0, c(1,0) );
  EXPECT_TRUE( is_valid(c(0,0)) );
  EXPECT_FALSE( is_valid(c(1,0)) );
}

TEST_F( MaskedViewGrayTest, create_mask ) {
  ImageView<Px> non(2,1);
  non(0,0) = 10;
  non(1,0) = 0;
  ImageView<MPx> c = create_mask(non);
  EXPECT_EQ( 10, c(0,0) );
  EXPECT_EQ( 0, c(1,0) );
  EXPECT_TRUE( is_valid(c(0,0)) );
  EXPECT_FALSE( is_valid(c(1,0)) );

  non(1,0) = 20;
  c = create_mask(non,20);
  EXPECT_EQ( 10, c(0,0) );
  EXPECT_EQ( 0, c(1,0) );
  EXPECT_TRUE( is_valid(c(0,0)) );
  EXPECT_FALSE( is_valid(c(1,0)) );
}

TEST_F( MaskedViewU8Test, apply_mask ) {
  ImageView<Px> c = apply_mask(a);
  EXPECT_EQ( 2, c.cols() ); EXPECT_EQ( 2, c.rows() );
  EXPECT_EQ( 1, c(0,0) );
  EXPECT_EQ( 2, c(1,0) );
  EXPECT_EQ( 0, c(0,1) );
  EXPECT_EQ( 0, c(1,1) );
  ImageView<Px> d = apply_mask(a,200);
  EXPECT_EQ( 1, d(0,0) );
  EXPECT_EQ( 2, d(1,0) );
  EXPECT_EQ( 200, d(0,1) );
  EXPECT_EQ( 200, d(1,1) );
}

TEST_F( MaskedViewGrayTest, apply_mask ) {
  ImageView<Px> c = apply_mask(a);
  EXPECT_EQ( 2, c.cols() ); EXPECT_EQ( 2, c.rows() );
  EXPECT_EQ( 1, c(0,0) );
  EXPECT_EQ( 2, c(1,0) );
  EXPECT_EQ( 0, c(0,1) );
  EXPECT_EQ( 0, c(1,1) );
  ImageView<Px> d = apply_mask(a,200);
  EXPECT_EQ( 1, d(0,0) );
  EXPECT_EQ( 2, d(1,0) );
  EXPECT_EQ( 200, d(0,1) );
  EXPECT_EQ( 200, d(1,1) );
}

TEST_F( MaskedViewU8Test, copy_mask ) {
  ImageView<Px> in(2,2);
  in(0,0) = Px(9); in(1,0) = Px(8);
  in(0,1) = Px(1); in(1,1) = Px(4);
  ImageView<MPx> c = copy_mask(in,a);
  EXPECT_EQ( 2, c.cols() ); EXPECT_EQ( 2, c.rows() );
  EXPECT_TRUE( is_valid(c(0,0)) );
  EXPECT_TRUE( is_valid(c(1,0)) );
  EXPECT_FALSE(is_valid(c(0,1)) );
  EXPECT_FALSE(is_valid(c(1,1)) );
  EXPECT_EQ( 9, c(0,0) );
  EXPECT_EQ( 8, c(1,0) );
  EXPECT_EQ( 1, c(0,1) );
  EXPECT_EQ( 4, c(1,1) );
}

TEST_F( MaskedViewGrayTest, copy_mask ) {
  ImageView<Px> in(2,2);
  in(0,0) = Px(9); in(1,0) = Px(8);
  in(0,1) = Px(1); in(1,1) = Px(4);
  ImageView<MPx> c = copy_mask(in,a);
  EXPECT_EQ( 2, c.cols() ); EXPECT_EQ( 2, c.rows() );
  EXPECT_TRUE( is_valid(c(0,0)) );
  EXPECT_TRUE( is_valid(c(1,0)) );
  EXPECT_FALSE(is_valid(c(0,1)) );
  EXPECT_FALSE(is_valid(c(1,1)) );
  EXPECT_EQ( 9, c(0,0) );
  EXPECT_EQ( 8, c(1,0) );
  EXPECT_EQ( 1, c(0,1) );
  EXPECT_EQ( 4, c(1,1) );
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

TEST_F( MaskedViewU8Test, edge_mask ) {
  ImageView<Px> c(4,4);
  fill( c, Px(1) );
  c(0,0) = Px(); c(0,1) = Px(); c(0,2) = Px(); c(2,1) = Px();
  ImageView<MPx> d = edge_mask( c );
  EXPECT_FALSE( is_valid(d(0,0)) );
  EXPECT_FALSE( is_valid(d(0,2)) );
  EXPECT_TRUE( is_valid(d(2,1)) );
  EXPECT_TRUE( is_valid(d(2,2)) );
  EXPECT_FALSE( is_valid(d(3,3)) );
  EXPECT_FALSE( is_valid(d(0,3)) );
  EXPECT_FALSE( is_valid(d(3,0)) );
}

TEST_F( MaskedViewGrayTest, edge_mask ) {
  ImageView<Px> c(4,4);
  fill( c, Px(1) );
  c(0,0) = Px(); c(0,1) = Px(); c(0,2) = Px(); c(2,1) = Px();
  ImageView<MPx> d = edge_mask( c );
  EXPECT_FALSE( is_valid(d(0,0)) );
  EXPECT_FALSE( is_valid(d(0,2)) );
  EXPECT_TRUE( is_valid(d(2,1)) );
  EXPECT_TRUE( is_valid(d(2,2)) );
  EXPECT_FALSE( is_valid(d(3,3)) );
  EXPECT_FALSE( is_valid(d(0,3)) );
  EXPECT_FALSE( is_valid(d(3,0)) );
}

TEST_F( MaskedViewU8Test, validate_mask ) {
  ImageView<MPx> c = validate_mask(a);
  EXPECT_EQ( 2, c.cols() ); EXPECT_EQ( 2, c.rows() );
  EXPECT_TRUE( is_valid(c(0,0)) );
  EXPECT_TRUE( is_valid(c(1,0)) );
  EXPECT_TRUE( is_valid(c(0,1)) );
  EXPECT_TRUE( is_valid(c(1,0)) );
  EXPECT_EQ( 1, c(0,0) );
  EXPECT_EQ( 4, c(1,1) );
}

TEST_F( MaskedViewGrayTest, validate_mask ) {
  ImageView<MPx> c = validate_mask(a);
  EXPECT_EQ( 2, c.cols() ); EXPECT_EQ( 2, c.rows() );
  EXPECT_TRUE( is_valid(c(0,0)) );
  EXPECT_TRUE( is_valid(c(1,0)) );
  EXPECT_TRUE( is_valid(c(0,1)) );
  EXPECT_TRUE( is_valid(c(1,0)) );
  EXPECT_EQ( 1, c(0,0) );
  EXPECT_EQ( 4, c(1,1) );
}

TEST_F( MaskedViewU8Test, invalidate_mask ) {
  ImageView<MPx> c = invalidate_mask(a);
  EXPECT_EQ( 2, c.cols() ); EXPECT_EQ( 2, c.rows() );
  EXPECT_FALSE( is_valid(c(0,0)) );
  EXPECT_FALSE( is_valid(c(1,0)) );
  EXPECT_FALSE( is_valid(c(0,1)) );
  EXPECT_FALSE( is_valid(c(1,0)) );
  EXPECT_EQ( 1, c(0,0) );
  EXPECT_EQ( 4, c(1,1) );
}

TEST_F( MaskedViewGrayTest, invalidate_mask ) {
  ImageView<MPx> c = invalidate_mask(a);
  EXPECT_EQ( 2, c.cols() ); EXPECT_EQ( 2, c.rows() );
  EXPECT_FALSE( is_valid(c(0,0)) );
  EXPECT_FALSE( is_valid(c(1,0)) );
  EXPECT_FALSE( is_valid(c(0,1)) );
  EXPECT_FALSE( is_valid(c(1,0)) );
  EXPECT_EQ( 1, c(0,0) );
  EXPECT_EQ( 4, c(1,1) );
}

TEST_F( MaskedViewU8Test, union_mask ) {
  ImageView<MPx> c = union_mask(a,b);
  EXPECT_EQ( 2, c.cols() ); EXPECT_EQ( 2, c.rows() );
  EXPECT_FALSE( is_valid(c(1,1)) );
  EXPECT_TRUE(  is_valid(c(0,0)) );
  EXPECT_TRUE(  is_valid(c(0,1)) );
  EXPECT_TRUE(  is_valid(c(1,0)) );
  EXPECT_EQ( 1, c(0,0) );
  EXPECT_EQ( 2, c(1,0) );
  EXPECT_EQ( 3, c(0,1) );
  EXPECT_EQ( 4, c(1,1) );
}

TEST_F( MaskedViewGrayTest, union_mask ) {
  ImageView<MPx> c = union_mask(a,b);
  EXPECT_EQ( 2, c.cols() ); EXPECT_EQ( 2, c.rows() );
  EXPECT_FALSE( is_valid(c(1,1)) );
  EXPECT_TRUE(  is_valid(c(0,0)) );
  EXPECT_TRUE(  is_valid(c(0,1)) );
  EXPECT_TRUE(  is_valid(c(1,0)) );
  EXPECT_EQ( 1, c(0,0) );
  EXPECT_EQ( 2, c(1,0) );
  EXPECT_EQ( 3, c(0,1) );
  EXPECT_EQ( 4, c(1,1) );
}

TEST_F( MaskedViewU8Test, intersect_mask ) {
  ImageView<MPx> c = intersect_mask(a,b);
  EXPECT_EQ( 2, c.cols() ); EXPECT_EQ( 2, c.rows() );
  EXPECT_FALSE( is_valid(c(1,1)) );
  EXPECT_FALSE( is_valid(c(1,0)) );
  EXPECT_FALSE( is_valid(c(0,1)) );
  EXPECT_TRUE(  is_valid(c(0,0)) );
  EXPECT_EQ( 1, c(0,0) );
  EXPECT_EQ( 2, c(1,0) );
  EXPECT_EQ( 3, c(0,1) );
  EXPECT_EQ( 4, c(1,1) );
}

TEST_F( MaskedViewGrayTest, intersect_mask ) {
  ImageView<MPx> c = intersect_mask(a,b);
  EXPECT_EQ( 2, c.cols() ); EXPECT_EQ( 2, c.rows() );
  EXPECT_FALSE( is_valid(c(1,1)) );
  EXPECT_FALSE( is_valid(c(1,0)) );
  EXPECT_FALSE( is_valid(c(0,1)) );
  EXPECT_TRUE(  is_valid(c(0,0)) );
  EXPECT_EQ( 1, c(0,0) );
  EXPECT_EQ( 2, c(1,0) );
  EXPECT_EQ( 3, c(0,1) );
  EXPECT_EQ( 4, c(1,1) );
}

TEST_F( MaskedViewU8Test, invert_mask ) {
  b = invert_mask(a);
  EXPECT_EQ( 2, b.cols() ); EXPECT_EQ( 2, b.rows() );
  EXPECT_NE( is_valid(a(0,0)), is_valid(b(0,0)) );
  EXPECT_NE( is_valid(a(1,0)), is_valid(b(1,0)) );
  EXPECT_NE( is_valid(a(0,1)), is_valid(b(0,1)) );
  EXPECT_NE( is_valid(a(1,1)), is_valid(b(1,1)) );
}

TEST_F( MaskedViewGrayTest, invert_mask ) {
  b = invert_mask(a);
  EXPECT_EQ( 2, b.cols() ); EXPECT_EQ( 2, b.rows() );
  EXPECT_NE( is_valid(a(0,0)), is_valid(b(0,0)) );
  EXPECT_NE( is_valid(a(1,0)), is_valid(b(1,0)) );
  EXPECT_NE( is_valid(a(0,1)), is_valid(b(0,1)) );
  EXPECT_NE( is_valid(a(1,1)), is_valid(b(1,1)) );
}
