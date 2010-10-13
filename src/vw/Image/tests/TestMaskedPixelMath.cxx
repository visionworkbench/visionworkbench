// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestMaskedPixelMath.h
#include <gtest/gtest.h>

#include <vw/Image/PixelMath.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/Interpolation.h>

#include <test/Helpers.h>

using namespace vw;

template <template<class> class TraitT, class T>
static bool bool_trait( T const& /*arg*/ ) {
  return TraitT<T>::value;
}

template <class PixelT>
class MaskedBinaryTest : public ::testing::Test {
protected:
  typedef PixelT Px;
  typedef typename PixelChannelType<Px>::type ChT;
  typedef PixelMask<PixelT> MPx;
  MPx a, b;

  MaskedBinaryTest() {}

  virtual void SetUp() {
    int ci = 1, cj = 2;
    a = MPx();
    b = MPx();
    a.validate();
    b.validate();
    for ( size_t i = 0; i < CompoundNumChannels<PixelT>::value; i++ ) {
      a[i] = ci;
      b[i] = cj;
      ci++;
      cj *= 2;
    }
  }

  void test_validness() {
    MPx c = a+b;
    MPx i;

    // Check boolean trait
    ASSERT_TRUE( bool_trait<IsMasked>( c ) );
    ASSERT_TRUE( bool_trait<IsMasked>( i ) );
    ASSERT_FALSE( bool_trait<IsMasked>( Px() ) );
    ASSERT_FALSE( bool_trait<IsMasked>( ChT() ) );

    // Check validity
    EXPECT_TRUE( is_valid(c) );
    EXPECT_FALSE( is_valid(i) );

    // Check all applications of sum operator
    c = a + i;
    EXPECT_FALSE( is_valid(c) );
    c = i + a;
    EXPECT_FALSE( is_valid(c) );
    c = i + i;
    EXPECT_FALSE( is_valid(c) );
    // Against Scalar
    c = a + 6;
    EXPECT_TRUE( is_valid(c) );
    c = 7 + b;
    EXPECT_TRUE( is_valid(c) );
    c = i + 6;
    EXPECT_FALSE( is_valid(c) );
    c = 2 + i;
    EXPECT_FALSE( is_valid(c) );
    c = a;
    c += 2;
    EXPECT_TRUE( is_valid(c) );
    c = i;
    c += 2;
    EXPECT_FALSE( is_valid(c) );
    // Against Masked Scalar
    c = a + PixelMask<ChT>(6);
    EXPECT_TRUE( is_valid(c) );
    c = PixelMask<ChT>(6) + b;
    EXPECT_TRUE( is_valid(c) );
    c = i + PixelMask<ChT>(3);
    EXPECT_FALSE( is_valid(c) );
    c = PixelMask<ChT>(3) + i;
    EXPECT_FALSE( is_valid(c) );
    c = a;
    c += PixelMask<ChT>(2);
    EXPECT_TRUE( is_valid(c) );
    c = i;
    c += PixelMask<ChT>(3);
    EXPECT_FALSE( is_valid(c) );

  }
};

typedef MaskedBinaryTest<PixelRGB<uint8> > MaskedRGBu8;
TEST_F( MaskedRGBu8, Arithmetic ) {
  test_validness();

  // Checking values
  MPx c = a+b;
  MPx i;
  EXPECT_PIXEL_EQ( c.child(), Px(3,6,11) );
  c = a + i;
  EXPECT_PIXEL_EQ( c.child(), Px(1,2,3) );
  c = i + i;
  EXPECT_PIXEL_EQ( c.child(), Px() );
  c = a - b;
  EXPECT_PIXEL_EQ( c.child(), Px(-1,-2,-5) );
  c = a - i;
  EXPECT_PIXEL_EQ( c.child(), Px(1,2,3) );
  c = i - a;
  EXPECT_PIXEL_EQ( c.child(), Px(-1,-2,-3) );
  c = i - i;
  EXPECT_PIXEL_EQ( c.child(), Px() );
  c = a * b;
  EXPECT_PIXEL_EQ( c.child(), Px( 2, 8, 24 ) );
  c = a * i;
  EXPECT_PIXEL_EQ( c.child(), Px() );
  c = i * i;
  EXPECT_PIXEL_EQ( c.child(), Px() );
  c = a / b;
  EXPECT_PIXEL_EQ( c.child(), Px(0.5,0.5,3.0/8.0) );
  c = i / b;
  EXPECT_PIXEL_EQ( c.child(), Px() );
  //  c = b / i; // THROWS FLOATING EXCEPTION (divide by zero)
  //  EXPECT_PIXEL_EQ( c.child(), Px() );
}

typedef MaskedBinaryTest<PixelRGB<float32> > MaskedRGBf32;
TEST_F( MaskedRGBf32, Arithmetic ) {
  test_validness();

  // Checking values
  MPx c = a+b;
  MPx i;
  EXPECT_PIXEL_EQ( c.child(), Px(3,6,11) );
  c = a + i;
  EXPECT_PIXEL_EQ( c.child(), Px(1,2,3) );
  c = i + i;
  EXPECT_PIXEL_EQ( c.child(), Px() );
  c = a - b;
  EXPECT_PIXEL_EQ( c.child(), Px(-1,-2,-5) );
  c = a - i;
  EXPECT_PIXEL_EQ( c.child(), Px(1,2,3) );
  c = i - a;
  EXPECT_PIXEL_EQ( c.child(), Px(-1,-2,-3) );
  c = i - i;
  EXPECT_PIXEL_EQ( c.child(), Px() );
  c = a * b;
  EXPECT_PIXEL_EQ( c.child(), Px( 2, 8, 24 ) );
  c = a * i;
  EXPECT_PIXEL_EQ( c.child(), Px() );
  c = i * i;
  EXPECT_PIXEL_EQ( c.child(), Px() );
  c = a / b;
  EXPECT_PIXEL_EQ( c.child(), Px(0.5,0.5,3.0/8.0) );
  c = i / b;
  EXPECT_PIXEL_EQ( c.child(), Px() );
}

typedef MaskedBinaryTest<uint8> MaskedScalaru8;
TEST_F( MaskedScalaru8, Arithmetic ) {
  test_validness();

  // Checking values
  MPx c = a+b;
  MPx i;
  EXPECT_EQ( c.child(), 3 );
  c = a + i;
  EXPECT_EQ( c.child(), 1 );
  c = i + i;
  EXPECT_EQ( c.child(), 0 );
  c = a - b;
  EXPECT_EQ( c.child(), 255 ); // WRAP
  c = a - i;
  EXPECT_EQ( c.child(), 1 );
  c = i - a;
  EXPECT_EQ( c.child(), 255 ); // WRAP
  c = i - i;
  EXPECT_EQ( c.child(), 0 );
  c = a * b;
  EXPECT_EQ( c.child(), 2 );
  c = a * i;
  EXPECT_EQ( c.child(), 0 );
  c = i * a;
  EXPECT_EQ( c.child(), 0 );
  c = a / b;
  EXPECT_EQ( c.child(), 0 );
  c = i / b;
  EXPECT_EQ( c.child(), 0 );
}

typedef MaskedBinaryTest<float32> MaskedScalarf32;
TEST_F( MaskedScalarf32, Arithmetic ) {
  test_validness();

  // Checking values
  MPx c = a+b;
  MPx i;
  EXPECT_EQ( c.child(), 3 );
  c = a + i;
  EXPECT_EQ( c.child(), 1 );
  c = i + i;
  EXPECT_EQ( c.child(), 0 );
  c = a - b;
  EXPECT_EQ( c.child(), -1 );
  c = a - i;
  EXPECT_EQ( c.child(), 1 );
  c = i - a;
  EXPECT_EQ( c.child(), -1 );
  c = i - i;
  EXPECT_EQ( c.child(), 0 );
  c = a * b;
  EXPECT_EQ( c.child(), 2 );
  c = a * i;
  EXPECT_EQ( c.child(), 0 );
  c = i * a;
  EXPECT_EQ( c.child(), 0 );
  c = a / b;
  EXPECT_EQ( c.child(), 0.5 );
  c = i / b;
  EXPECT_EQ( c.child(), 0 );
}

TEST( MaskedPixelMath, PixelMaskInterpolation ) {
  typedef PixelMask<PixelGray<float> > px_type;
  ImageView<px_type > test(2,2);
  test(0,0) = px_type(0);
  test(0,1) = px_type(0);
  test(0,0).invalidate();
  test(0,1).invalidate();
  test(1,0) = px_type(255);
  test(1,0).invalidate();
  test(1,1) = px_type(255);

  InterpolationView<EdgeExtensionView<ImageView<px_type>, ConstantEdgeExtension>, BilinearInterpolation> interp_test(edge_extend(test, ConstantEdgeExtension()));

  EXPECT_FALSE( is_valid( interp_test(0,0) ) );
  EXPECT_FALSE( is_valid( interp_test(0,1) ) );
  EXPECT_FALSE( is_valid( interp_test(1,0) ) );

  EXPECT_TRUE( is_valid( interp_test(1,1) ) );

  EXPECT_FALSE( is_valid( interp_test(0.5,0.5) ) );
  EXPECT_FALSE( is_valid( interp_test(0.5,0) ) );
  EXPECT_FALSE( is_valid( interp_test(0,0.5) ) );
}

TEST( MaskedPixelMath, MixedTypes ) {
  typedef PixelRGB<float> CPx;
  typedef PixelMask<CPx> MCPx;
  typedef PixelMask<float> MSx;

  // This is not possible in the current frame work
  //
  // We could allow these kind of operations in the future, but
  // it negates the reason for apply_mask and create_mask. It
  // would also cause a lot of boiler plating code.

  /*
  {
    CPx a(1, 2, 4);
    MSx b(3);
    MSx i;
    // Compound + MScalar
    CPx c = a + b;
    EXPECT_PIXEL_EQ( c, CPx(4,5,7) );
    c = a + i;
    EXPECT_PIXEL_EQ( c, CPx(1,2,4) );
    // MScalar + Compound
    c = b + a;
    EXPECT_PIXEL_EQ( c, CPx(4,5,7) );
    c = i + a;
    EXPECT_PIXEL_EQ( c, CPx(1,2,4) );
    // Compound += MScalar
    c += b;
    EXPECT_PIXEL_EQ( c, CPx(4,5,7) );
    c += i;
    EXPECT_PIXEL_EQ( c, CPx(4,5,7) );
  }
  */
  // MCompound + Scalar
  // Scalar + MCompound
  // MCompound += Scalar
  // MScalar + Scalar
  // Scalar + MScalar
  // MScalar += Scalar
  // Scalar += MScalar

}
