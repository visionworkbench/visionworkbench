// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

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

typedef ::testing::Types<PixelGray<uint8>, PixelGray<float>, PixelRGB<uint8>, PixelRGB<float>, uint8, float, Vector2f, Vector2i> MyTypes;

template <class PixelT>
struct MaskedImageMath : public ::testing::Test {
  typedef PixelT Px;
  typedef typename PixelChannelType<Px>::type ChT;
  typedef PixelMask<Px> MPx;
  ImageView<MPx> Av, Ai, Bv, Bi;

  void SetUp() {
    Av.set_size(1,1);
    Ai.set_size(1,1);
    Bv.set_size(1,1);
    Bi.set_size(1,1);
  }
};

TYPED_TEST_CASE( MaskedImageMath, MyTypes );
TYPED_TEST( MaskedImageMath, Sum ) {

  assignment( this->Av(0,0), 2 );
  assignment( this->Ai(0,0), 2 );
  assignment( this->Bv(0,0), 100 );
  assignment( this->Bi(0,0), 100 );
  this->Av(0,0).validate();
  this->Bv(0,0).validate();

  // Test traits
  ASSERT_TRUE(  bool_trait<IsMasked>( this->Av(0,0) ) );
  ASSERT_TRUE(  bool_trait<IsMasked>( this->Ai(0,0) ) );

  // Test Operations
  typedef typename TestFixture::MPx MPx;
  typedef typename TestFixture::Px Px;
  ImageView<MPx> F = this->Av + this->Bv;
  EXPECT_TRUE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(102), F(0,0).child() );
  F = this->Ai + this->Bv;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(102), F(0,0).child() );
  F = this->Av + this->Bi;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(102), F(0,0).child() );
  F = this->Ai + this->Bi;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(102), F(0,0).child() );
  F = this->Av + 2;
  EXPECT_TRUE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(4), F(0,0).child() );
  F = this->Ai + 2;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(4), F(0,0).child() );
  F = 2 + this->Bv;
  EXPECT_TRUE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(102), F(0,0).child()  );
  F = 2 + this->Bi;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(102), F(0,0).child() );
  F = copy(this->Av);
  F += this->Bv;
  EXPECT_TRUE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(102), F(0,0).child() );
  F = copy(this->Av);
  F += this->Bi;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(102), F(0,0).child() );
  F = copy(this->Ai);
  F += this->Bv;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(102), F(0,0).child() );
  F = copy(this->Ai);
  F += this->Bi;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(102), F(0,0).child() );
  F = copy(this->Av);
  F += 2;
  EXPECT_TRUE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(4), F(0,0).child() );
  F = copy(this->Ai);
  F += 2;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(4), F(0,0).child() );
}

TYPED_TEST( MaskedImageMath, Difference ) {

  this->Ai.set_size(1,1);
  assignment( this->Ai(0,0), 100 );
  this->Av.set_size(1,1);
  assignment( this->Av(0,0), 100 );
  this->Bi.set_size(1,1);
  assignment( this->Bi(0,0), 50 );
  this->Bv.set_size(1,1);
  assignment( this->Bv(0,0), 50 );
  this->Av(0,0).validate();
  this->Bv(0,0).validate();

  // Test traits
  ASSERT_TRUE(  bool_trait<IsMasked>( this->Av(0,0) ) );
  ASSERT_TRUE(  bool_trait<IsMasked>( this->Ai(0,0) ) );

  // Test Operations
  typedef typename TestFixture::MPx MPx;
  typedef typename TestFixture::Px Px;
  ImageView<MPx> F = this->Av - this->Bv;
  EXPECT_TRUE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  F = this->Ai - this->Bv;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  F = this->Av - this->Bi;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  F = this->Ai - this->Bi;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  F = this->Av - 2;
  EXPECT_TRUE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(98), F(0,0).child() );
  F = this->Ai - 2;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(98), F(0,0).child() );
  F = 100 - this->Bv;
  EXPECT_TRUE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  F = 100 - this->Bi;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  F = copy(this->Av);
  F -= this->Bv;
  EXPECT_TRUE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  F = copy(this->Av);
  F -= this->Bi;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  F = copy(this->Ai);
  F -= this->Bv;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  F = copy(this->Ai);
  F -= this->Bi;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  F = copy(this->Av);
  F -= 2;
  EXPECT_TRUE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(98), F(0,0).child() );
  F = copy(this->Ai);
  F -= 2;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(98), F(0,0).child() );
}

TYPED_TEST( MaskedImageMath, Product ) {
  assignment( this->Av(0,0), 100 );
  assignment( this->Ai(0,0), 100 );
  assignment( this->Bv(0,0), 2 );
  assignment( this->Bi(0,0), 2 );
  this->Av(0,0).validate();
  this->Bv(0,0).validate();

  // Test traits
  ASSERT_TRUE(  bool_trait<IsMasked>( this->Av(0,0) ) );
  ASSERT_TRUE(  bool_trait<IsMasked>( this->Ai(0,0) ) );

  // Test Operations
  typedef typename TestFixture::MPx MPx;
  typedef typename TestFixture::Px Px;
  ImageView<MPx> F = this->Av * this->Bv;
  EXPECT_TRUE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
  F = this->Ai * this->Bv;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
  F = this->Av * this->Bi;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
  F = this->Ai * this->Bi;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
  F = this->Av * 2;
  EXPECT_TRUE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
  F = this->Ai * 2;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
  F = 100 * this->Bv;
  EXPECT_TRUE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
  F = 100 * this->Bi;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
  F = copy(this->Av);
  F *= this->Bv;
  EXPECT_TRUE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
  F = copy(this->Av);
  F *= this->Bi;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
  F = copy(this->Ai);
  F *= this->Bv;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
  F = copy(this->Ai);
  F *= this->Bi;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
  F = copy(this->Av);
  F *= 2;
  EXPECT_TRUE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
  F = copy(this->Ai);
  F *= 2;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(200), F(0,0).child() );
}

TYPED_TEST( MaskedImageMath, Quotient ) {
  assignment( this->Av(0,0), 100 );
  assignment( this->Ai(0,0), 100 );
  assignment( this->Bv(0,0), 2 );
  assignment( this->Bi(0,0), 2 );
  this->Av(0,0).validate();
  this->Bv(0,0).validate();

  // Test Operations
  typedef typename TestFixture::MPx MPx;
  typedef typename TestFixture::Px Px;
  ImageView<MPx> F = this->Av / this->Bv;
  EXPECT_TRUE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  F = this->Ai / this->Bv;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  F = this->Av / this->Bi;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  F = this->Ai / this->Bi;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  F = this->Av / 2;
  EXPECT_TRUE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  F = this->Ai / 2;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  F = 100 / this->Bv;
  EXPECT_TRUE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  F = 100 / this->Bi;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  F = copy(this->Av);
  F /= this->Bv;
  EXPECT_TRUE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  F = copy(this->Av);
  F /= this->Bi;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  F = copy(this->Ai);
  F /= this->Bv;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  F = copy(this->Ai);
  F /= this->Bi;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  F = copy(this->Av);
  F /= 2;
  EXPECT_TRUE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
  F = copy(this->Ai);
  F /= 2;
  EXPECT_FALSE( is_valid(F(0,0)) );
  EXPECT_VW_EQ( construct<Px>(50), F(0,0).child() );
}
