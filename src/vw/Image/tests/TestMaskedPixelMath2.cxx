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


// TestMaskedPixelMath.h
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/PixelMath.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/Interpolation.h>
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
struct MaskedPixelMath : public ::testing::Test {
  typedef PixelT Px;
  typedef typename PixelChannelType<Px>::type ChT;
  typedef PixelMask<PixelT> MPx;
  MPx Av, Ai, Bv, Bi;
};

TYPED_TEST_CASE( MaskedPixelMath, MyTypes );

TYPED_TEST( MaskedPixelMath, Sum ) {
  assignment( this->Ai, 2 );
  assignment( this->Av, 2 );
  assignment( this->Bi, 100 );
  assignment( this->Bv, 100 );
  this->Av.validate();
  this->Bv.validate();

  // Test traits
  typedef typename TestFixture::MPx MPx;
  typedef typename TestFixture::Px Px;
  typedef typename TestFixture::ChT ChT;
  ASSERT_TRUE(  bool_trait<IsMasked>( this->Av ) );
  ASSERT_TRUE(  bool_trait<IsMasked>( this->Ai ) );
  ASSERT_FALSE( bool_trait<IsMasked>( Px() ) );
  ASSERT_FALSE( bool_trait<IsMasked>( ChT() ) );
  ASSERT_TRUE(  bool_trait<IsMasked>( MPx() ) );

  // Test Operations
  MPx F = this->Av + this->Bv;
  EXPECT_TRUE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(102), F.child() );
  F = this->Ai + this->Bv;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(102), F.child() );
  F = this->Av + this->Bi;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(102), F.child() );
  F = this->Ai + this->Bi;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(102), F.child() );
  F = this->Av + 2;
  EXPECT_TRUE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(4), F.child() );
  F = this->Ai + 2;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(4), F.child() );
  F = 2 + this->Bv;
  EXPECT_TRUE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(102), F.child() );
  F = 2 + this->Bi;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(102), F.child() );
  F = this->Av;
  F += this->Bv;
  EXPECT_TRUE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(102), F.child() );
  F = this->Av;
  F += this->Bi;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(102), F.child() );
  F = this->Ai;
  F += this->Bv;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(102), F.child() );
  F = this->Ai;
  F += this->Bi;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(102), F.child() );
  F = this->Av;
  F += 2;
  EXPECT_TRUE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(4), F.child() );
  F = this->Ai;
  F += 2;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(4), F.child() );
}

TYPED_TEST( MaskedPixelMath, Difference ) {
  assignment( this->Ai, 100 );
  assignment( this->Av, 100 );
  assignment( this->Bi, 50 );
  assignment( this->Bv, 50 );
  this->Av.validate();
  this->Bv.validate();

  // Test Operations
  typedef typename TestFixture::MPx MPx;
  typedef typename TestFixture::Px Px;
  MPx F = this->Av - this->Bv;
  EXPECT_TRUE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
  F = this->Ai - this->Bv;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
  F = this->Av - this->Bi;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
  F = this->Ai - this->Bi;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
  F = this->Av - 2;
  EXPECT_TRUE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(98), F.child() );
  F = this->Ai - 2;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(98), F.child() );
  F = 100 - this->Bv;
  EXPECT_TRUE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
  F = 100 - this->Bi;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
  F = this->Av;
  F -= this->Bv;
  EXPECT_TRUE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
  F = this->Av;
  F -= this->Bi;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
  F = this->Ai;
  F -= this->Bv;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
  F = this->Ai;
  F -= this->Bi;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
  F = this->Av;
  F -= 2;
  EXPECT_TRUE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(98), F.child() );
  F = this->Ai;
  F -= 2;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(98), F.child() );
}

TYPED_TEST( MaskedPixelMath, Product ) {
  assignment( this->Ai, 100 );
  assignment( this->Av, 100 );
  assignment( this->Bi, 2 );
  assignment( this->Bv, 2 );
  this->Av.validate();
  this->Bv.validate();

  // Test Operations
  typedef typename TestFixture::MPx MPx;
  typedef typename TestFixture::Px Px;
  MPx F = this->Av * this->Bv;
  EXPECT_TRUE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(200), F.child() );
  F = this->Ai * this->Bv;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(200), F.child() );
  F = this->Av * this->Bi;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(200), F.child() );
  F = this->Ai * this->Bi;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(200), F.child() );
  F = this->Av * 2;
  EXPECT_TRUE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(200), F.child() );
  F = this->Ai * 2;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(200), F.child() );
  F = 100 * this->Bv;
  EXPECT_TRUE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(200), F.child() );
  F = 100 * this->Bi;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(200), F.child() );
  F = this->Av;
  F *= this->Bv;
  EXPECT_TRUE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(200), F.child() );
  F = this->Av;
  F *= this->Bi;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(200), F.child() );
  F = this->Ai;
  F *= this->Bv;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(200), F.child() );
  F = this->Ai;
  F *= this->Bi;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(200), F.child() );
  F = this->Av;
  F *= 2;
  EXPECT_TRUE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(200), F.child() );
  F = this->Ai;
  F *= 2;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(200), F.child() );
}

TYPED_TEST( MaskedPixelMath, Quotient ) {
  assignment( this->Ai, 100 );
  assignment( this->Av, 100 );
  assignment( this->Bi, 2 );
  assignment( this->Bv, 2 );
  this->Av.validate();
  this->Bv.validate();

  // Test Operations
  typedef typename TestFixture::MPx MPx;
  typedef typename TestFixture::Px Px;
  MPx F = this->Av / this->Bv;
  EXPECT_TRUE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
  F = this->Ai / this->Bv;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
  F = this->Av / this->Bi;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
  F = this->Ai / this->Bi;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
  F = this->Av / 2;
  EXPECT_TRUE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
  F = this->Ai / 2;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
  F = 100 / this->Bv;
  EXPECT_TRUE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
  F = 100 / this->Bi;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
  F = this->Av;
  F /= this->Bv;
  EXPECT_TRUE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
  F = this->Av;
  F /= this->Bi;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
  F = this->Ai;
  F /= this->Bv;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
  F = this->Ai;
  F /= this->Bi;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
  F = this->Av;
  F /= 2;
  EXPECT_TRUE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
  F = this->Ai;
  F /= 2;
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(50), F.child() );
}

TYPED_TEST( MaskedPixelMath, Abs ) {
  assignment( this->Ai, 2 );
  assignment( this->Av, 2 );
  assignment( this->Bi, 100 );
  assignment( this->Bv, 100 );
  this->Av.validate();
  this->Bv.validate();

  EXPECT_TRUE( is_valid( this->Av ) );
  EXPECT_FALSE( is_valid( this->Ai ) );
  EXPECT_TRUE( is_valid( this->Bv ) );
  EXPECT_FALSE( is_valid( this->Bi ) );

  // Test traits
  typedef typename TestFixture::MPx MPx;
  typedef typename TestFixture::Px Px;
  typedef typename TestFixture::ChT ChT;
  ASSERT_TRUE(  bool_trait<IsMasked>( this->Av ) );
  ASSERT_TRUE(  bool_trait<IsMasked>( this->Ai ) );
  ASSERT_FALSE( bool_trait<IsMasked>( Px() ) );
  ASSERT_FALSE( bool_trait<IsMasked>( ChT() ) );
  ASSERT_TRUE(  bool_trait<IsMasked>( MPx() ) );

  MPx F = abs(this->Av);
  EXPECT_TRUE( bool_trait<IsMasked>( abs(this->Av) ) );
  EXPECT_TRUE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(2), F.child() );
  F = abs(this->Ai);
  EXPECT_TRUE( bool_trait<IsMasked>( abs(this->Ai) ) );
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(2), F.child() );
  F = abs(this->Bv);
  EXPECT_TRUE( bool_trait<IsMasked>( abs(this->Bv) ) );
  EXPECT_TRUE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(100), F.child() );
  F = abs(this->Bi);
  EXPECT_TRUE( bool_trait<IsMasked>( abs(this->Bi) ) );
  EXPECT_FALSE( is_valid(F) );
  EXPECT_VW_EQ( construct<Px>(100), F.child() );
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
