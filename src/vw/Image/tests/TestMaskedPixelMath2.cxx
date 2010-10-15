// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestMaskedPixelMath.h
#include <gtest/gtest.h>
#include <vw/Image.h>
#include <test/Helpers.h>

using namespace vw;
using boost::enable_if;
using boost::disable_if;

template <template<class> class TraitT, class T>
static bool bool_trait( T const& /*arg*/ ) {
  return TraitT<T>::value;
}

// Assignment, this in ensures a constant is applied to all elements
template <class T, class V>
void assignment( PixelMask<T> & a, V value ) {
  for ( size_t i = 0; i < CompoundNumChannels<PixelMask<T> >::value - 1; i++ )
    a[i] = value;
}

// Construction, this in ensures a constant is applied to all elements
template <class T, class V>
typename enable_if<IsCompound<T>, T>::type construct( V value ) {
  T a;
  for ( size_t i = 0; i < CompoundNumChannels<T>::value; i++ )
    a[i] = value;
  return a;
}
template <class T, class V>
typename disable_if<IsCompound<T>, T>::type construct( V value ) {
  return T(value);
}


template <class PixelT>
struct MaskedSum : public ::testing::Test {
  typedef PixelT Px;
  typedef typename PixelChannelType<Px>::type ChT;
  typedef PixelMask<PixelT> MPx;
  MPx Av, Ai, Bv, Bi;

  void SetUp() {
    assignment( Ai, 2 );
    assignment( Av, 2 );
    assignment( Bi, 100 );
    assignment( Bv, 100 );
    Av.validate();
    Bv.validate();
  }

  void test() {
    // Test traits
    ASSERT_TRUE(  bool_trait<IsMasked>( Av ) );
    ASSERT_TRUE(  bool_trait<IsMasked>( Ai ) );
    ASSERT_FALSE( bool_trait<IsMasked>( Px() ) );
    ASSERT_FALSE( bool_trait<IsMasked>( ChT() ) );
    ASSERT_TRUE(  bool_trait<IsMasked>( MPx() ) );

    // Test Operations
    MPx F = Av + Bv;
    EXPECT_TRUE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(102), F.child() );
    F = Ai + Bv;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(102), F.child() );
    F = Av + Bi;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(102), F.child() );
    F = Ai + Bi;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(102), F.child() );
    F = Av + 2;
    EXPECT_TRUE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(4), F.child() );
    F = Ai + 2;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(4), F.child() );
    F = 2 + Bv;
    EXPECT_TRUE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(102), F.child() );
    F = 2 + Bi;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(102), F.child() );
    F = Av;
    F += Bv;
    EXPECT_TRUE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(102), F.child() );
    F = Av;
    F += Bi;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(102), F.child() );
    F = Ai;
    F += Bv;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(102), F.child() );
    F = Ai;
    F += Bi;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(102), F.child() );
    F = Av;
    F += 2;
    EXPECT_TRUE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(4), F.child() );
    F = Ai;
    F += 2;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(4), F.child() );
  }
};

template <class PixelT>
struct MaskedDifference : public ::testing::Test {
  typedef PixelT Px;
  typedef typename PixelChannelType<Px>::type ChT;
  typedef PixelMask<PixelT> MPx;
  MPx Av, Ai, Bv, Bi;

  void SetUp() {
    assignment( Ai, 100 );
    assignment( Av, 100 );
    assignment( Bi, 50 );
    assignment( Bv, 50 );
    Av.validate();
    Bv.validate();
  }

  void test() {
    // Test Operations
    MPx F = Av - Bv;
    EXPECT_TRUE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
    F = Ai - Bv;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
    F = Av - Bi;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
    F = Ai - Bi;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
    F = Av - 2;
    EXPECT_TRUE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(98), F.child() );
    F = Ai - 2;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(98), F.child() );
    F = 100 - Bv;
    EXPECT_TRUE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
    F = 100 - Bi;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
    F = Av;
    F -= Bv;
    EXPECT_TRUE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
    F = Av;
    F -= Bi;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
    F = Ai;
    F -= Bv;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
    F = Ai;
    F -= Bi;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
    F = Av;
    F -= 2;
    EXPECT_TRUE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(98), F.child() );
    F = Ai;
    F -= 2;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(98), F.child() );
  }
};

template <class PixelT>
struct MaskedProduct : public ::testing::Test {
  typedef PixelT Px;
  typedef typename PixelChannelType<Px>::type ChT;
  typedef PixelMask<PixelT> MPx;
  MPx Av, Ai, Bv, Bi;

  void SetUp() {
    assignment( Ai, 100 );
    assignment( Av, 100 );
    assignment( Bi, 2 );
    assignment( Bv, 2 );
    Av.validate();
    Bv.validate();
  }

  void test() {
    // Test Operations
    MPx F = Av * Bv;
    EXPECT_TRUE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(200), F.child() );
    F = Ai * Bv;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(200), F.child() );
    F = Av * Bi;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(200), F.child() );
    F = Ai * Bi;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(200), F.child() );
    F = Av * 2;
    EXPECT_TRUE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(200), F.child() );
    F = Ai * 2;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(200), F.child() );
    F = 100 * Bv;
    EXPECT_TRUE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(200), F.child() );
    F = 100 * Bi;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(200), F.child() );
    F = Av;
    F *= Bv;
    EXPECT_TRUE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(200), F.child() );
    F = Av;
    F *= Bi;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(200), F.child() );
    F = Ai;
    F *= Bv;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(200), F.child() );
    F = Ai;
    F *= Bi;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(200), F.child() );
    F = Av;
    F *= 2;
    EXPECT_TRUE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(200), F.child() );
    F = Ai;
    F *= 2;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(200), F.child() );
  }
};

template <class PixelT>
struct MaskedQuotient : public ::testing::Test {
  typedef PixelT Px;
  typedef typename PixelChannelType<Px>::type ChT;
  typedef PixelMask<PixelT> MPx;
  MPx Av, Ai, Bv, Bi;

  void SetUp() {
    assignment( Ai, 100 );
    assignment( Av, 100 );
    assignment( Bi, 2 );
    assignment( Bv, 2 );
    Av.validate();
    Bv.validate();
  }

  void test() {
    // Test Operations
    MPx F = Av / Bv;
    EXPECT_TRUE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
    F = Ai / Bv;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
    F = Av / Bi;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
    F = Ai / Bi;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
    F = Av / 2;
    EXPECT_TRUE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
    F = Ai / 2;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
    F = 100 / Bv;
    EXPECT_TRUE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
    F = 100 / Bi;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
    F = Av;
    F /= Bv;
    EXPECT_TRUE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
    F = Av;
    F /= Bi;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
    F = Ai;
    F /= Bv;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
    F = Ai;
    F /= Bi;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
    F = Av;
    F /= 2;
    EXPECT_TRUE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
    F = Ai;
    F /= 2;
    EXPECT_FALSE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(50), F.child() );
  }
};

// Instantiate tests
typedef MaskedSum<PixelGray<uint8> > MaskedSumGrayU8;
TEST_F( MaskedSumGrayU8, PixelMath ) { test(); }
typedef MaskedSum<PixelGray<float> > MaskedSumGrayF32;
TEST_F( MaskedSumGrayF32, PixelMath ) { test(); }
typedef MaskedSum<PixelRGB<uint8> >  MaskedSumRGBU8;
TEST_F( MaskedSumRGBU8, PixelMath ) { test(); }
typedef MaskedSum<PixelRGB<float> >  MaskedSumRGBF32;
TEST_F( MaskedSumRGBF32, PixelMath ) { test(); }
typedef MaskedSum<uint8> MaskedSumU8;
TEST_F( MaskedSumU8, PixelMath ) { test(); }
typedef MaskedSum<float> MaskedSumF32;
TEST_F( MaskedSumF32, PixelMath ) { test(); }
typedef MaskedSum<Vector2f> MaskedSumVec2f;
TEST_F( MaskedSumVec2f, PixelMath ) { test(); }
typedef MaskedSum<Vector2i> MaskedSumVec2i;
TEST_F( MaskedSumVec2i, PixelMath ) { test(); }

typedef MaskedDifference<PixelGray<uint8> > MaskedDifferenceGrayU8;
TEST_F( MaskedDifferenceGrayU8, PixelMath ) { test(); }
typedef MaskedDifference<PixelGray<float> > MaskedDifferenceGrayF32;
TEST_F( MaskedDifferenceGrayF32, PixelMath ) { test(); }
typedef MaskedDifference<PixelRGB<uint8> >  MaskedDifferenceRGBU8;
TEST_F( MaskedDifferenceRGBU8, PixelMath ) { test(); }
typedef MaskedDifference<PixelRGB<float> >  MaskedDifferenceRGBF32;
TEST_F( MaskedDifferenceRGBF32, PixelMath ) { test(); }
typedef MaskedDifference<uint8> MaskedDifferenceU8;
TEST_F( MaskedDifferenceU8, PixelMath ) { test(); }
typedef MaskedDifference<float> MaskedDifferenceF32;
TEST_F( MaskedDifferenceF32, PixelMath ) { test(); }
typedef MaskedDifference<Vector2f> MaskedDifferenceVec2f;
TEST_F( MaskedDifferenceVec2f, PixelMath ) { test(); }
typedef MaskedDifference<Vector2i> MaskedDifferenceVec2i;
TEST_F( MaskedDifferenceVec2i, PixelMath ) { test(); }

typedef MaskedProduct<PixelGray<uint8> > MaskedProductGrayU8;
TEST_F( MaskedProductGrayU8, PixelMath ) { test(); }
typedef MaskedProduct<PixelGray<float> > MaskedProductGrayF32;
TEST_F( MaskedProductGrayF32, PixelMath ) { test(); }
typedef MaskedProduct<PixelRGB<uint8> >  MaskedProductRGBU8;
TEST_F( MaskedProductRGBU8, PixelMath ) { test(); }
typedef MaskedProduct<PixelRGB<float> >  MaskedProductRGBF32;
TEST_F( MaskedProductRGBF32, PixelMath ) { test(); }
typedef MaskedProduct<uint8> MaskedProductU8;
TEST_F( MaskedProductU8, PixelMath ) { test(); }
typedef MaskedProduct<float> MaskedProductF32;
TEST_F( MaskedProductF32, PixelMath ) { test(); }
typedef MaskedProduct<Vector2f> MaskedProductVec2f;
TEST_F( MaskedProductVec2f, PixelMath ) { test(); }
typedef MaskedProduct<Vector2i> MaskedProductVec2i;
TEST_F( MaskedProductVec2i, PixelMath ) { test(); }

typedef MaskedQuotient<PixelGray<uint8> > MaskedQuotientGrayU8;
TEST_F( MaskedQuotientGrayU8, PixelMath ) { test(); }
typedef MaskedQuotient<PixelGray<float> > MaskedQuotientGrayF32;
TEST_F( MaskedQuotientGrayF32, PixelMath ) { test(); }
typedef MaskedQuotient<PixelRGB<uint8> >  MaskedQuotientRGBU8;
TEST_F( MaskedQuotientRGBU8, PixelMath ) { test(); }
typedef MaskedQuotient<PixelRGB<float> >  MaskedQuotientRGBF32;
TEST_F( MaskedQuotientRGBF32, PixelMath ) { test(); }
typedef MaskedQuotient<uint8> MaskedQuotientU8;
TEST_F( MaskedQuotientU8, PixelMath ) { test(); }
typedef MaskedQuotient<float> MaskedQuotientF32;
TEST_F( MaskedQuotientF32, PixelMath ) { test(); }
typedef MaskedQuotient<Vector2f> MaskedQuotientVec2f;
TEST_F( MaskedQuotientVec2f, PixelMath ) { test(); }
typedef MaskedQuotient<Vector2i> MaskedQuotientVec2i;
TEST_F( MaskedQuotientVec2i, PixelMath ) { test(); }

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
