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


// TestPixelMath.h
#include <gtest/gtest_VW.h>

#include <vw/config.h>
#include <vw/Image/PixelMath.h>
#include <vw/Image/PixelTypes.h>

#include <test/Helpers.h>

// Create a simple toy pixel type using the PixelMathBase
// class to test its default functionality in isolation.
namespace vw {
  template <class ChannelT>
  class ToyType : public PixelMathBase<ToyType<ChannelT> > {
    ChannelT m_val;
  public:
    ToyType() : m_val() {}
    ToyType( ChannelT val ) : m_val(val) {}
    ChannelT& operator[](int) { return m_val; }
    ChannelT const& operator[](int) const { return m_val; }
  };

VW_DECLARE_PIXEL_TYPE(ToyType,1);
}

using namespace vw;

template <class T1, class T2>
static bool is_of_type( T2 ) {
  return boost::is_same<T1,T2>::value;
}

#define TEST_UNARY_MATH_OPERATOR(op,arg,result)                         \
  do {                                                                  \
    ToyType<double> a(arg);                                             \
    EXPECT_NEAR( (op a)[0], (result), 0.00001 );                    \
    ASSERT_TRUE( is_of_type<ToyType<float> >( (op ToyType<float>(1)) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<double> >( (op ToyType<double>(1)) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<long double> >( (op ToyType<long double>(1)) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<int> >( (op ToyType<int>(1)) ) );     \
  } while(false)

#define TEST_BINARY_MATH_OPERATOR(op,arg1,arg2,result)                  \
  do {                                                                  \
    ToyType<double> a(arg1), b(arg2);                                   \
    EXPECT_NEAR( (a op b)[0], (result), 0.00001 );                  \
    EXPECT_NEAR( (a op double(arg2))[0], (result), 0.00001 );       \
    EXPECT_NEAR( (double(arg1) op b)[0], (result), 0.00001 );       \
    ASSERT_TRUE( is_of_type<ToyType<float> >( (ToyType<float>(1) op ToyType<float>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<double> >( (ToyType<double>(1) op ToyType<double>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<long double> >( (ToyType<long double>(1) op ToyType<long double>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<double> >( (ToyType<float>(1) op ToyType<double>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<long double> >( (ToyType<long double>(1) op ToyType<double>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<float> >( (ToyType<float>(1) op ToyType<int>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<int> >( (ToyType<int>(1) op ToyType<int>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<float> >( (ToyType<float>(1) op float(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<double> >( (ToyType<double>(1) op double(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<long double> >( (ToyType<long double>(1) op (long double)(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<double> >( (ToyType<float>(1) op double(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<long double> >( (ToyType<long double>(1) op double(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<float> >( (ToyType<float>(1) op int(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<int> >( (ToyType<int>(1) op int(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<float> >( (float(1) op ToyType<float>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<double> >( (double(1) op ToyType<double>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<long double> >( ((long double)(1) op ToyType<long double>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<double> >( (float(1) op ToyType<double>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<long double> >( ((long double)(1) op ToyType<double>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<float> >( (float(1) op ToyType<int>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<int> >( (int(1) op ToyType<int>(1) ) ) ); \
  } while(false)

TEST( PixelMath, Negation ) { TEST_UNARY_MATH_OPERATOR(-,1,-1); }
TEST( PixelMath, Sum ) { TEST_BINARY_MATH_OPERATOR(+,2,3,5); }
TEST( PixelMath, Difference ) { TEST_BINARY_MATH_OPERATOR(-,2,3,-1); }
TEST( PixelMath, Product ) { TEST_BINARY_MATH_OPERATOR(*,2,3,6); }
TEST( PixelMath, Quotient ) { TEST_BINARY_MATH_OPERATOR(/,3,2,1.5); }

#define TEST_UNARY_MATH_FUNCTION(name,arg,result)                       \
  do {                                                                  \
    ToyType<double> a(arg);                                             \
    EXPECT_NEAR( name(a)[0], (result), 0.00001 );                   \
    ASSERT_TRUE( is_of_type<ToyType<float> >( name(ToyType<float>(1)) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<double> >( name(ToyType<double>(1)) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<long double> >( name(ToyType<long double>(1)) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<double> >( name(ToyType<int>(1)) ) ); \
  } while(false)

#define TEST_BINARY_MATH_FUNCTION(name,arg1,arg2,result)                \
  do {                                                                  \
    ToyType<double> a(arg1), b(arg2);                                   \
    EXPECT_NEAR( name(a,b)[0], (result), 0.00001 );                 \
    EXPECT_NEAR( name(a,double(arg2))[0], (result), 0.00001 );      \
    EXPECT_NEAR( name(double(arg1),b)[0], (result), 0.00001 );      \
    ASSERT_TRUE( is_of_type<ToyType<float> >( name(ToyType<float>(1),ToyType<float>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<double> >( name(ToyType<double>(1),ToyType<double>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<long double> >( name(ToyType<long double>(1),ToyType<long double>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<double> >( name(ToyType<float>(1),ToyType<double>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<long double> >( name(ToyType<long double>(1),ToyType<double>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<float> >( name(ToyType<float>(1),ToyType<int>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<double> >( name(ToyType<int>(1),ToyType<int>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<float> >( name(ToyType<float>(1),float(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<double> >( name(ToyType<double>(1),double(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<long double> >( name(ToyType<long double>(1),(long double)(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<double> >( name(ToyType<float>(1),double(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<long double> >( name(ToyType<long double>(1),double(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<float> >( name(ToyType<float>(1),int(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<double> >( name(ToyType<int>(1),int(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<float> >( name(float(1),ToyType<float>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<double> >( name(double(1),ToyType<double>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<long double> >( name((long double)(1),ToyType<long double>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<double> >( name(float(1),ToyType<double>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<long double> >( name((long double)(1),ToyType<double>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<float> >( name(float(1),ToyType<int>(1) ) ) ); \
    ASSERT_TRUE( is_of_type<ToyType<double> >( name(int(1),ToyType<int>(1) ) ) ); \
  } while(false)

TEST( PixelMath, SelfAssignment ) {
  typedef PixelRGB<uint8> Px;
  Px a(1,2,3), b(2,3,4);
  a += 2;
  EXPECT_VW_EQ( a, Px(3,4,5) );
  a -= 2;
  EXPECT_VW_EQ( a, Px(1,2,3) );
  a *= 2;
  EXPECT_VW_EQ( a, Px(2,4,6) );
  a /= 2;
  EXPECT_VW_EQ( a, Px(1,2,3) );
  a += b;
  EXPECT_VW_EQ( a, Px(3,5,7) );
  a -= b;
  EXPECT_VW_EQ( a, Px(1,2,3) );
  a *= b;
  EXPECT_VW_EQ( a, Px(2,6,12) );
  a /= b;
  EXPECT_VW_EQ( a, Px(1,2,3) );
  (a += 2) += 2;
  EXPECT_VW_EQ( a, Px(5,6,7) );
  (a -= 2) -= 2;
  EXPECT_VW_EQ( a, Px(1,2,3) );
  (a *= 2) *= 2;
  EXPECT_VW_EQ( a, Px(4,8,12) );
  (a /= 2) /= 2;
  EXPECT_VW_EQ( a, Px(1,2,3) );
}

TEST( PixelMath, ACOS ) { TEST_UNARY_MATH_FUNCTION(acos,0.5,1.0472); }
TEST( PixelMath, ASIN ) { TEST_UNARY_MATH_FUNCTION(asin,0.5,0.523599); }
TEST( PixelMath, ATAN ) { TEST_UNARY_MATH_FUNCTION(atan,1.0,0.785398); }
TEST( PixelMath, COS ) { TEST_UNARY_MATH_FUNCTION(cos,1.0,0.540302); }
TEST( PixelMath, SIN ) { TEST_UNARY_MATH_FUNCTION(sin,1.0,0.841471); }
TEST( PixelMath, TAN ) { TEST_UNARY_MATH_FUNCTION(tan,1.0,1.55741); }
TEST( PixelMath, COSH ) { TEST_UNARY_MATH_FUNCTION(cosh,1.0,1.54308); }
TEST( PixelMath, SINH ) { TEST_UNARY_MATH_FUNCTION(sinh,1.0,1.1752); }
TEST( PixelMath, TANH ) { TEST_UNARY_MATH_FUNCTION(tanh,1.0,0.761594); }
TEST( PixelMath, EXP ) { TEST_UNARY_MATH_FUNCTION(exp,1.0,2.718281); }
TEST( PixelMath, LOG ) { TEST_UNARY_MATH_FUNCTION(log,2.0,0.693147); }
TEST( PixelMath, LOG10 ) { TEST_UNARY_MATH_FUNCTION(log10,2.0,0.30103); }
TEST( PixelMath, SQRT ) { TEST_UNARY_MATH_FUNCTION(sqrt,2.0,1.41421); }
TEST( PixelMath, CEIL ) { TEST_UNARY_MATH_FUNCTION(ceil,1.5,2.0);
  TEST_UNARY_MATH_FUNCTION(ceil,-1.5,-1.0); }
TEST( PixelMath, FLOOR ) { TEST_UNARY_MATH_FUNCTION(floor,1.5,1.0);
  TEST_UNARY_MATH_FUNCTION(floor,-1.5,-2.0); }

void test_atan2() { TEST_BINARY_MATH_FUNCTION(atan2,2.0,1.0,1.10715); }
void test_pow() { TEST_BINARY_MATH_FUNCTION(pow,3.0,2.0,9.0); }

#ifdef WIN32
#undef TEST_UNARY_MATH_FUNCTION
#define TEST_UNARY_MATH_FUNCTION(name,arg,result) do {} while(false)
#undef TEST_BINARY_MATH_FUNCTION
#define TEST_BINARY_MATH_FUNCTION(name,arg1,arg2,result) do {} while(false)
#endif

TEST( PixelMath, ACOSH ) { TEST_UNARY_MATH_FUNCTION(acosh,1.5,0.962424); }
TEST( PixelMath, ASINH ) { TEST_UNARY_MATH_FUNCTION(asinh,1.0,0.881374); }
TEST( PixelMath, ATANH ) { TEST_UNARY_MATH_FUNCTION(atanh,0.5,0.549306); }
TEST( PixelMath, EXP2 ) {
#ifdef VW_HAVE_EXP2
  TEST_UNARY_MATH_FUNCTION(exp2,1.0,2.0);
#endif
}
TEST( PixelMath, EXPM1 ) { TEST_UNARY_MATH_FUNCTION(expm1,1.0,1.718281); }
TEST( PixelMath, LOG2 ) {
#ifdef VW_HAVE_LOG2
  TEST_UNARY_MATH_FUNCTION(log2,2.0,1.0);
#endif
}
TEST( PixelMath, LOG1P ) { TEST_UNARY_MATH_FUNCTION(log1p,1.0,0.693147); }
TEST( PixelMath, CBRT ) { TEST_UNARY_MATH_FUNCTION(cbrt,2.0,1.25992); }
TEST( PixelMath, ERF ) { TEST_UNARY_MATH_FUNCTION(erf,1.0,0.842701); }
TEST( PixelMath, ERFC ) { TEST_UNARY_MATH_FUNCTION(erfc,1.0,0.157299); }
TEST( PixelMath, TGAMMA )  {
#ifdef VW_HAVE_TGAMMA
  TEST_UNARY_MATH_FUNCTION(tgamma,1.5,0.886227);
#endif
}
TEST( PixelMath, LGAMMA ) { TEST_UNARY_MATH_FUNCTION(lgamma,2.5,0.284683); }
TEST( PixelMath, ROUND ) { TEST_UNARY_MATH_FUNCTION(round,1.4,1.0);
  TEST_UNARY_MATH_FUNCTION(round,1.5,2.0); }
TEST( PixelMath, TRUNC ) { TEST_UNARY_MATH_FUNCTION(trunc,1.5,1.0);
  TEST_UNARY_MATH_FUNCTION(trunc,-1.5,-1.0); }

TEST( PixelMath, HYPOT ) { TEST_BINARY_MATH_FUNCTION(hypot,2.0,1.0,2.23607); }
TEST( PixelMath, COPYSIGN ) { TEST_BINARY_MATH_FUNCTION(copysign,3.0,-2.0,-3.0);
  TEST_BINARY_MATH_FUNCTION(copysign,3.0,2.0,3.0); }
TEST( PixelMath, FDIM ) { TEST_BINARY_MATH_FUNCTION(fdim,3.0,2.0,1.0);
  TEST_BINARY_MATH_FUNCTION(fdim,2.0,3.0,0.0); }

#define ASSERT_PRESERVED_TYPE( op ) \
  ASSERT_TRUE( is_of_type<ToyType<float> >( op(ToyType<float>(1)) ) ); \
  ASSERT_TRUE( is_of_type<ToyType<double> >( op(ToyType<double>(1)) ) ); \
  ASSERT_TRUE( is_of_type<ToyType<long double> >( op(ToyType<long double>(1)) ) ); \
  ASSERT_TRUE( is_of_type<ToyType<int> >( op(ToyType<int>(1)) ) ); \
  ASSERT_TRUE( is_of_type<ToyType<float> >( op(ToyType<std::complex<float> >(1)) ) ); \
  ASSERT_TRUE( is_of_type<ToyType<double> >( op(ToyType<std::complex<double> >(1)) ) ); \
  ASSERT_TRUE( is_of_type<ToyType<long double> >( op(ToyType<std::complex<long double> >(1)) ) ); \
  ASSERT_TRUE( is_of_type<ToyType<int> >( op(ToyType<std::complex<int> >(1)) ) );

TEST( PixelMath, Real ) {
  ToyType<double> ar(1.0);
  ToyType<std::complex<double> > ac(std::complex<double>(2,3));
  EXPECT_EQ( real(ToyType<double>(1.0))[0], 1.0 );
  EXPECT_EQ( real(ToyType<std::complex<double> >(std::complex<double>(2,3)))[0], 2.0 );
  ASSERT_PRESERVED_TYPE( real );
}

TEST( PixelMath, Imaginary ) {
  EXPECT_EQ( imag(ToyType<double>(1.0))[0], 0.0 );
  EXPECT_EQ( imag(ToyType<std::complex<double> >(std::complex<double>(2,3)))[0], 3.0 );
  ASSERT_PRESERVED_TYPE( imag );
}

TEST( PixelMath, Abs ) {
  EXPECT_EQ( abs(ToyType<double>(1.0))[0], 1.0 );
  EXPECT_EQ( abs(ToyType<double>(-1.0))[0], 1.0 );
  EXPECT_NEAR( abs(ToyType<std::complex<double> >(std::complex<double>(3,4)))[0], 5.0, 0.00001 );
  ASSERT_PRESERVED_TYPE( abs );
}

TEST( PixelMath, Conj ) {
  EXPECT_EQ( conj(ToyType<double>(1.0))[0], 1.0 );
  EXPECT_EQ( conj(ToyType<double>(-1.0))[0], -1.0 );
  EXPECT_EQ( conj(ToyType<std::complex<double> >(std::complex<double>(3,4)))[0], std::complex<double>(3,-4) );
  ASSERT_TRUE( is_of_type<ToyType<float> >( conj(ToyType<float>(1)) ) );
  ASSERT_TRUE( is_of_type<ToyType<double> >( conj(ToyType<double>(1)) ) );
  ASSERT_TRUE( is_of_type<ToyType<long double> >( conj(ToyType<long double>(1)) ) );
  ASSERT_TRUE( is_of_type<ToyType<int> >( conj(ToyType<int>(1)) ) );
  ASSERT_TRUE( is_of_type<ToyType<std::complex<float> > >( conj(ToyType<std::complex<float> >(1)) ) );
  ASSERT_TRUE( is_of_type<ToyType<std::complex<double> > >( conj(ToyType<std::complex<double> >(1)) ) );
  ASSERT_TRUE( is_of_type<ToyType<std::complex<long double> > >( conj(ToyType<std::complex<long double> >(1)) ) );
  ASSERT_TRUE( is_of_type<ToyType<std::complex<int> > >( conj(ToyType<std::complex<int> >(1)) ) );
}

TEST( PixelMath, Square ) {
  EXPECT_EQ( square(ToyType<double>(1.0))[0], 1.0 );
  EXPECT_EQ( square(ToyType<double>(-1.0))[0], 1.0 );
  EXPECT_EQ( square(ToyType<int32>(32))[0], 1024 );
  EXPECT_EQ( square(ToyType<uint8>(2u))[0], 4 );
  ASSERT_TRUE( is_of_type<ToyType<float> >( square(ToyType<float>(1)) ) );
  ASSERT_TRUE( is_of_type<ToyType<int32> >( square(ToyType<int32>(1)) ) );
}

TEST( PixelMath, AcosFunctor ) {
  ToyType<double> x(0.5);
  ASSERT_TRUE( is_of_type<ToyType<double> >( vw::math::ArgAcosFunctor()(x) ) );
  vw::math::ArgAcosFunctor f;  // Pulled out to workaround gcc 3.2 bug
  EXPECT_NEAR( f(x)[0], acos(x[0]), 1e-8 );
}

TEST( PixelMath, ArgArgHypotFunctor ) {
  ToyType<double> x(3), y(4);
  ASSERT_TRUE( is_of_type<ToyType<double> >( vw::math::ArgArgHypotFunctor()(x,y) ) );
  vw::math::ArgArgHypotFunctor f;  // Pulled out to workaround gcc 3.2 bug
  EXPECT_NEAR( f(x,y)[0], hypot(x[0],y[0]), 1e-8 );
}

TEST( PixelMath, ArgValHypotFunctor ) {
  ToyType<double> x(3), y(4);
  ASSERT_TRUE( is_of_type<ToyType<double> >( vw::math::ArgValHypotFunctor<ToyType<double> >(y)(x) ) );
  EXPECT_NEAR( vw::math::ArgValHypotFunctor<ToyType<double> >(y)(x)[0], hypot(x[0],y[0]), 1e-8 );
}

TEST( PixelMath, ValArgHypotFunctor ) {
  ToyType<double> x(3), y(4);
  ASSERT_TRUE( ( boost::is_same< ToyType<double>, boost::result_of<vw::math::ValArgHypotFunctor<ToyType<double> >(ToyType<double>)>::type >::value ) );
  ASSERT_TRUE( is_of_type<ToyType<double> >( vw::math::ValArgHypotFunctor<ToyType<double> >(x)(y) ) );
  EXPECT_NEAR( vw::math::ValArgHypotFunctor<ToyType<double> >(x)(y)[0], hypot(x[0],y[0]), 1e-8 );
}

