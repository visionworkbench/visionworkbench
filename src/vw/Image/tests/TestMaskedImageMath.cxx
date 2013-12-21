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

#include <test/Helpers.h>
#include <vw/Core/CompoundTypes.h>      // for CompoundNumChannels, etc
#include <vw/Core/FundamentalTypes.h>   // for uint8
#include <vw/Math/Vector.h>             // for operator==, Vector2f, etc
#include <vw/Image/ImageMath.h>         // for operator*, operator+, etc
#include <vw/Image/ImageView.h>         // for ImageView
#include <vw/Image/Manipulation.h>      // for copy
#include <vw/Image/PixelMask.h>         // for is_valid, IsMasked, etc
#include <vw/Image/PixelMath.h>         // for operator*, operator+, etc
#include <vw/Image/PixelTypeInfo.h>     // for PixelChannelType
#include <vw/Image/PixelTypes.h>        // for operator<<, etc

#include <boost/utility/enable_if.hpp>  // for disable_if, enable_if

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
