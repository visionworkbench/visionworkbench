// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

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
  typedef PixelMask<Px> MPx;
  ImageView<MPx> Av, Ai, Bv, Bi;

  void SetUp() {
    Av.set_size(1,1);
    assignment( Av(0,0), 2 );
    Ai.set_size(1,1);
    assignment( Ai(0,0), 2 );
    Bv.set_size(1,1);
    assignment( Bv(0,0), 100 );
    Bi.set_size(1,1);
    assignment( Bi(0,0), 100 );
    Av(0,0).validate();
    Bv(0,0).validate();
  }

  void test() {
    // Test traits
    ASSERT_TRUE(  bool_trait<IsMasked>( Av(0,0) ) );
    ASSERT_TRUE(  bool_trait<IsMasked>( Ai(0,0) ) );

    // Test Operations
    ImageView<MPx> F = Av + Bv;
    EXPECT_TRUE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(102), F(0,0).child() );
    F = Ai + Bv;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(102), F(0,0).child() );
    F = Av + Bi;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(102), F(0,0).child() );
    F = Ai + Bi;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(102), F(0,0).child() );
    F = Av + 2;
    EXPECT_TRUE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(4), F(0,0).child() );
    F = Ai + 2;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(4), F(0,0).child() );
    F = 2 + Bv;
    EXPECT_TRUE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(102), F(0,0).child()  );
    F = 2 + Bi;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(102), F(0,0).child() );
    F = copy(Av);
    F += Bv;
    EXPECT_TRUE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(102), F(0,0).child() );
    F = copy(Av);
    F += Bi;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(102), F(0,0).child() );
    F = copy(Ai);
    F += Bv;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(102), F(0,0).child() );
    F = copy(Ai);
    F += Bi;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(102), F(0,0).child() );
    F = copy(Av);
    F += 2;
    EXPECT_TRUE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(4), F(0,0).child() );
    F = copy(Ai);
    F += 2;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(4), F(0,0).child() );
  }
};

template <class PixelT>
struct MaskedDifference : public ::testing::Test {
  typedef PixelT Px;
  typedef typename PixelChannelType<Px>::type ChT;
  typedef PixelMask<Px> MPx;
  ImageView<MPx> Av, Ai, Bv, Bi;

  void SetUp() {
    Ai.set_size(1,1);
    assignment( Ai(0,0), 100 );
    Av.set_size(1,1);
    assignment( Av(0,0), 100 );
    Bi.set_size(1,1);
    assignment( Bi(0,0), 50 );
    Bv.set_size(1,1);
    assignment( Bv(0,0), 50 );
    Av(0,0).validate();
    Bv(0,0).validate();
  }

  void test() {
    // Test traits
    ASSERT_TRUE(  bool_trait<IsMasked>( Av(0,0) ) );
    ASSERT_TRUE(  bool_trait<IsMasked>( Ai(0,0) ) );

    // Test Operations
    ImageView<MPx> F = Av - Bv;
    EXPECT_TRUE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
    F = Ai - Bv;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
    F = Av - Bi;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
    F = Ai - Bi;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
    F = Av - 2;
    EXPECT_TRUE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(98), F(0,0).child() );
    F = Ai - 2;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(98), F(0,0).child() );
    F = 100 - Bv;
    EXPECT_TRUE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
    F = 100 - Bi;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
    F = copy(Av);
    F -= Bv;
    EXPECT_TRUE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
    F = copy(Av);
    F -= Bi;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
    F = copy(Ai);
    F -= Bv;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
    F = copy(Ai);
    F -= Bi;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
    F = copy(Av);
    F -= 2;
    EXPECT_TRUE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(98), F(0,0).child() );
    F = copy(Ai);
    F -= 2;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(98), F(0,0).child() );
  }
};

template <class PixelT>
struct MaskedProduct : public ::testing::Test {
  typedef PixelT Px;
  typedef typename PixelChannelType<Px>::type ChT;
  typedef PixelMask<Px> MPx;
  ImageView<MPx> Av, Ai, Bv, Bi;

  void SetUp() {
    Av.set_size(1,1);
    assignment( Av(0,0), 100 );
    Ai.set_size(1,1);
    assignment( Ai(0,0), 100 );
    Bv.set_size(1,1);
    assignment( Bv(0,0), 2 );
    Bi.set_size(1,1);
    assignment( Bi(0,0), 2 );
    Av(0,0).validate();
    Bv(0,0).validate();
  }

  void test() {
    // Test traits
    ASSERT_TRUE(  bool_trait<IsMasked>( Av(0,0) ) );
    ASSERT_TRUE(  bool_trait<IsMasked>( Ai(0,0) ) );

    // Test Operations
    ImageView<MPx> F = Av * Bv;
    EXPECT_TRUE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
    F = Ai * Bv;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
    F = Av * Bi;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
    F = Ai * Bi;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
    F = Av * 2;
    EXPECT_TRUE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
    F = Ai * 2;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
    F = 100 * Bv;
    EXPECT_TRUE( is_valid(F) );
    EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
    F = 100 * Bi;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
    F = copy(Av);
    F *= Bv;
    EXPECT_TRUE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
    F = copy(Av);
    F *= Bi;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
    F = copy(Ai);
    F *= Bv;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
    F = copy(Ai);
    F *= Bi;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
    F = copy(Av);
    F *= 2;
    EXPECT_TRUE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
    F = copy(Ai);
    F *= 2;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
  }
};

template <class PixelT>
struct MaskedQuotient : public ::testing::Test {
  typedef PixelT Px;
  typedef typename PixelChannelType<Px>::type ChT;
  typedef PixelMask<Px> MPx;
  ImageView<MPx> Av, Ai, Bv, Bi;

  void SetUp() {
    Av.set_size(1,1);
    assignment( Av(0,0), 100 );
    Ai.set_size(1,1);
    assignment( Ai(0,0), 100 );
    Bv.set_size(1,1);
    assignment( Bv(0,0), 2 );
    Bi.set_size(1,1);
    assignment( Bi(0,0), 2 );
    Av(0,0).validate();
    Bv(0,0).validate();
  }

  void test() {
    // Test Operations
    ImageView<MPx> F = Av / Bv;
    EXPECT_TRUE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
    F = Ai / Bv;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
    F = Av / Bi;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
    F = Ai / Bi;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
    F = Av / 2;
    EXPECT_TRUE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
    F = Ai / 2;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
    F = 100 / Bv;
    EXPECT_TRUE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
    F = 100 / Bi;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
    F = copy(Av);
    F /= Bv;
    EXPECT_TRUE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
    F = copy(Av);
    F /= Bi;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
    F = copy(Ai);
    F /= Bv;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
    F = copy(Ai);
    F /= Bi;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
    F = copy(Av);
    F /= 2;
    EXPECT_TRUE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
    F = copy(Ai);
    F /= 2;
    EXPECT_FALSE( is_valid(F(0,0)) );
    EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  }
};

// Instantiate tests
typedef MaskedSum<PixelGray<uint8> > MaskedSumGrayU8;
TEST_F( MaskedSumGrayU8, ImageMath ) { test(); }
typedef MaskedSum<PixelGray<float> > MaskedSumGrayF32;
TEST_F( MaskedSumGrayF32, ImageMath ) { test(); }
typedef MaskedSum<PixelRGB<uint8> >  MaskedSumRGBU8;
TEST_F( MaskedSumRGBU8, ImageMath ) { test(); }
typedef MaskedSum<PixelRGB<float> >  MaskedSumRGBF32;
TEST_F( MaskedSumRGBF32, ImageMath ) { test(); }
typedef MaskedSum<uint8> MaskedSumU8;
TEST_F( MaskedSumU8, ImageMath ) { test(); }
typedef MaskedSum<float> MaskedSumF32;
TEST_F( MaskedSumF32, ImageMath ) { test(); }
typedef MaskedSum<Vector2f> MaskedSumVec2f;
TEST_F( MaskedSumVec2f, ImageMath ) { test(); }
typedef MaskedSum<Vector2i> MaskedSumVec2i;
TEST_F( MaskedSumVec2i, ImageMath ) { test(); }

typedef MaskedDifference<PixelGray<uint8> > MaskedDifferenceGrayU8;
TEST_F( MaskedDifferenceGrayU8, ImageMath ) { test(); }
typedef MaskedDifference<PixelGray<float> > MaskedDifferenceGrayF32;
TEST_F( MaskedDifferenceGrayF32, ImageMath ) { test(); }
typedef MaskedDifference<PixelRGB<uint8> >  MaskedDifferenceRGBU8;
TEST_F( MaskedDifferenceRGBU8, ImageMath ) { test(); }
typedef MaskedDifference<PixelRGB<float> >  MaskedDifferenceRGBF32;
TEST_F( MaskedDifferenceRGBF32, ImageMath ) { test(); }
typedef MaskedDifference<uint8> MaskedDifferenceU8;
TEST_F( MaskedDifferenceU8, ImageMath ) { test(); }
typedef MaskedDifference<float> MaskedDifferenceF32;
TEST_F( MaskedDifferenceF32, ImageMath ) { test(); }
typedef MaskedDifference<Vector2f> MaskedDifferenceVec2f;
TEST_F( MaskedDifferenceVec2f, ImageMath ) { test(); }
typedef MaskedDifference<Vector2i> MaskedDifferenceVec2i;
TEST_F( MaskedDifferenceVec2i, ImageMath ) { test(); }

typedef MaskedProduct<PixelGray<uint8> > MaskedProductGrayU8;
TEST_F( MaskedProductGrayU8, ImageMath ) { test(); }
typedef MaskedProduct<PixelGray<float> > MaskedProductGrayF32;
TEST_F( MaskedProductGrayF32, ImageMath ) { test(); }
typedef MaskedProduct<PixelRGB<uint8> >  MaskedProductRGBU8;
TEST_F( MaskedProductRGBU8, ImageMath ) { test(); }
typedef MaskedProduct<PixelRGB<float> >  MaskedProductRGBF32;
TEST_F( MaskedProductRGBF32, ImageMath ) { test(); }
typedef MaskedProduct<uint8> MaskedProductU8;
TEST_F( MaskedProductU8, ImageMath ) { test(); }
typedef MaskedProduct<float> MaskedProductF32;
TEST_F( MaskedProductF32, ImageMath ) { test(); }
typedef MaskedProduct<Vector2f> MaskedProductVec2f;
TEST_F( MaskedProductVec2f, ImageMath ) { test(); }
typedef MaskedProduct<Vector2i> MaskedProductVec2i;
TEST_F( MaskedProductVec2i, ImageMath ) { test(); }

typedef MaskedQuotient<PixelGray<uint8> > MaskedQuotientGrayU8;
TEST_F( MaskedQuotientGrayU8, ImageMath ) { test(); }
typedef MaskedQuotient<PixelGray<float> > MaskedQuotientGrayF32;
TEST_F( MaskedQuotientGrayF32, ImageMath ) { test(); }
typedef MaskedQuotient<PixelRGB<uint8> >  MaskedQuotientRGBU8;
TEST_F( MaskedQuotientRGBU8, ImageMath ) { test(); }
typedef MaskedQuotient<PixelRGB<float> >  MaskedQuotientRGBF32;
TEST_F( MaskedQuotientRGBF32, ImageMath ) { test(); }
typedef MaskedQuotient<uint8> MaskedQuotientU8;
TEST_F( MaskedQuotientU8, ImageMath ) { test(); }
typedef MaskedQuotient<float> MaskedQuotientF32;
TEST_F( MaskedQuotientF32, ImageMath ) { test(); }
typedef MaskedQuotient<Vector2f> MaskedQuotientVec2f;
TEST_F( MaskedQuotientVec2f, ImageMath ) { test(); }
typedef MaskedQuotient<Vector2i> MaskedQuotientVec2i;
TEST_F( MaskedQuotientVec2i, ImageMath ) { test(); }
