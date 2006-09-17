// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
//
// __END_LICENSE__

// TestPixelMath.h
#include <cxxtest/TestSuite.h>

#include <vw/Image/PixelMath.h>
#include <vw/Image/PixelTypeInfo.h>

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

class TestPixelMath : public CxxTest::TestSuite
{
public:

  template <class T1, class T2>
  static bool is_of_type( T2 ) {
    return boost::is_same<T1,T2>::value;
  }

#define TEST_UNARY_MATH_OPERATOR(op,arg,result)                         \
  do {                                                                  \
    ToyType<double> a(arg);                                             \
    TS_ASSERT_DELTA( (op a)[0], (result), 0.00001 );                    \
    TS_ASSERT( is_of_type<ToyType<float> >( (op ToyType<float>()) ) );  \
    TS_ASSERT( is_of_type<ToyType<double> >( (op ToyType<double>()) ) ); \
    TS_ASSERT( is_of_type<ToyType<long double> >( (op ToyType<long double>()) ) ); \
    TS_ASSERT( is_of_type<ToyType<int> >( (op ToyType<int>()) ) );      \
  } while(false)

#define TEST_BINARY_MATH_OPERATOR(op,arg1,arg2,result)                  \
  do {                                                                  \
    ToyType<double> a(arg1), b(arg2);                                   \
    TS_ASSERT_DELTA( (a op b)[0], (result), 0.00001 );                  \
    TS_ASSERT_DELTA( (a op double(arg2))[0], (result), 0.00001 );       \
    TS_ASSERT_DELTA( (double(arg1) op b)[0], (result), 0.00001 );       \
    TS_ASSERT( is_of_type<ToyType<float> >( (ToyType<float>() op ToyType<float>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<double> >( (ToyType<double>() op ToyType<double>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<long double> >( (ToyType<long double>() op ToyType<long double>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<double> >( (ToyType<float>() op ToyType<double>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<long double> >( (ToyType<long double>() op ToyType<double>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<float> >( (ToyType<float>() op ToyType<int>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<int> >( (ToyType<int>() op ToyType<int>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<float> >( (ToyType<float>() op float() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<double> >( (ToyType<double>() op double() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<long double> >( (ToyType<long double>() op (long double)(0) ) ) ); \
    TS_ASSERT( is_of_type<ToyType<double> >( (ToyType<float>() op double() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<long double> >( (ToyType<long double>() op double() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<float> >( (ToyType<float>() op int() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<int> >( (ToyType<int>() op int() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<float> >( (float() op ToyType<float>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<double> >( (double() op ToyType<double>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<long double> >( ((long double)(0) op ToyType<long double>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<double> >( (float() op ToyType<double>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<long double> >( ((long double)(0) op ToyType<double>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<float> >( (float() op ToyType<int>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<int> >( (int() op ToyType<int>() ) ) ); \
  } while(false)

  void test_netation() { TEST_UNARY_MATH_OPERATOR(-,1,-1); }
  void test_sum() { TEST_BINARY_MATH_OPERATOR(+,2,3,5); }
  void test_difference() { TEST_BINARY_MATH_OPERATOR(-,2,3,-1); }
  void test_product() { TEST_BINARY_MATH_OPERATOR(*,2,3,6); }
  void test_quotient() { TEST_BINARY_MATH_OPERATOR(/,3,2,1.5); }

#define TEST_UNARY_MATH_FUNCTION(name,arg,result)                       \
  do {                                                                  \
    ToyType<double> a(arg);                                             \
    TS_ASSERT_DELTA( name(a)[0], (result), 0.00001 );                   \
    TS_ASSERT( is_of_type<ToyType<float> >( name(ToyType<float>()) ) ); \
    TS_ASSERT( is_of_type<ToyType<double> >( name(ToyType<double>()) ) ); \
    TS_ASSERT( is_of_type<ToyType<long double> >( name(ToyType<long double>()) ) ); \
    TS_ASSERT( is_of_type<ToyType<double> >( name(ToyType<int>()) ) );  \
  } while(false)

#define TEST_BINARY_MATH_FUNCTION(name,arg1,arg2,result)                \
  do {                                                                  \
    ToyType<double> a(arg1), b(arg2);                                   \
    TS_ASSERT_DELTA( name(a,b)[0], (result), 0.00001 );                 \
    TS_ASSERT_DELTA( name(a,double(arg2))[0], (result), 0.00001 );      \
    TS_ASSERT_DELTA( name(double(arg1),b)[0], (result), 0.00001 );      \
    TS_ASSERT( is_of_type<ToyType<float> >( name(ToyType<float>(),ToyType<float>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<double> >( name(ToyType<double>(),ToyType<double>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<long double> >( name(ToyType<long double>(),ToyType<long double>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<double> >( name(ToyType<float>(),ToyType<double>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<long double> >( name(ToyType<long double>(),ToyType<double>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<float> >( name(ToyType<float>(),ToyType<int>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<double> >( name(ToyType<int>(),ToyType<int>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<float> >( name(ToyType<float>(),float() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<double> >( name(ToyType<double>(),double() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<long double> >( name(ToyType<long double>(),(long double)(0) ) ) ); \
    TS_ASSERT( is_of_type<ToyType<double> >( name(ToyType<float>(),double() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<long double> >( name(ToyType<long double>(),double() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<float> >( name(ToyType<float>(),int() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<double> >( name(ToyType<int>(),int() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<float> >( name(float(),ToyType<float>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<double> >( name(double(),ToyType<double>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<long double> >( name((long double)(0),ToyType<long double>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<double> >( name(float(),ToyType<double>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<long double> >( name((long double)(0),ToyType<double>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<float> >( name(float(),ToyType<int>() ) ) ); \
    TS_ASSERT( is_of_type<ToyType<double> >( name(int(),ToyType<int>() ) ) ); \
  } while(false)

  void test_acos() { TEST_UNARY_MATH_FUNCTION(acos,0.5,1.0472); }
  void test_asin() { TEST_UNARY_MATH_FUNCTION(asin,0.5,0.523599); }
  void test_atan() { TEST_UNARY_MATH_FUNCTION(atan,1.0,0.785398); }
  void test_cos() { TEST_UNARY_MATH_FUNCTION(cos,1.0,0.540302); }
  void test_sin() { TEST_UNARY_MATH_FUNCTION(sin,1.0,0.841471); }
  void test_tan() { TEST_UNARY_MATH_FUNCTION(tan,1.0,1.55741); }
  void test_cosh() { TEST_UNARY_MATH_FUNCTION(cosh,1.0,1.54308); }
  void test_sinh() { TEST_UNARY_MATH_FUNCTION(sinh,1.0,1.1752); }
  void test_tanh() { TEST_UNARY_MATH_FUNCTION(tanh,1.0,0.761594); }
  void test_exp() { TEST_UNARY_MATH_FUNCTION(exp,1.0,2.718281); }
  void test_log() { TEST_UNARY_MATH_FUNCTION(log,2.0,0.693147); }
  void test_log10() { TEST_UNARY_MATH_FUNCTION(log10,2.0,0.30103); }
  void test_sqrt() { TEST_UNARY_MATH_FUNCTION(sqrt,2.0,1.41421); }
  void test_ceil() { TEST_UNARY_MATH_FUNCTION(ceil,1.5,2.0);
                     TEST_UNARY_MATH_FUNCTION(ceil,-1.5,-1.0); }
  void test_floor() { TEST_UNARY_MATH_FUNCTION(floor,1.5,1.0);
                      TEST_UNARY_MATH_FUNCTION(floor,-1.5,-2.0); }

  void test_atan2() { TEST_BINARY_MATH_FUNCTION(atan2,2.0,1.0,1.10715); }
  void test_pow() { TEST_BINARY_MATH_FUNCTION(pow,3.0,2.0,9.0); }

#ifdef WIN32
#undef TEST_UNARY_MATH_FUNCTION
#define TEST_UNARY_MATH_FUNCTION(name,arg,result) do {} while(false)
#undef TEST_BINARY_MATH_FUNCTION
#define TEST_BINARY_MATH_FUNCTION(name,arg1,arg2,result) do {} while(false)
#endif

  void test_acosh() { TEST_UNARY_MATH_FUNCTION(acosh,1.5,0.962424); }
  void test_asinh() { TEST_UNARY_MATH_FUNCTION(asinh,1.0,0.881374); }
  void test_atanh() { TEST_UNARY_MATH_FUNCTION(atanh,0.5,0.549306); }
  void test_exp2() { TEST_UNARY_MATH_FUNCTION(exp2,1.0,2.0); }
  void test_expm1() { TEST_UNARY_MATH_FUNCTION(expm1,1.0,1.718281); }
  void test_log2() { TEST_UNARY_MATH_FUNCTION(log2,2.0,1.0); }
  void test_log1p() { TEST_UNARY_MATH_FUNCTION(log1p,1.0,0.693147); }
  void test_cbrt() { TEST_UNARY_MATH_FUNCTION(cbrt,2.0,1.25992); }
  void test_erf() { TEST_UNARY_MATH_FUNCTION(erf,1.0,0.842701); }
  void test_erfc() { TEST_UNARY_MATH_FUNCTION(erfc,1.0,0.157299); }
  void test_tgamma() { TEST_UNARY_MATH_FUNCTION(tgamma,1.5,0.886227); }
  void test_lgamma() { TEST_UNARY_MATH_FUNCTION(lgamma,2.5,0.284683); }
  void test_round() { TEST_UNARY_MATH_FUNCTION(round,1.4,1.0);
                      TEST_UNARY_MATH_FUNCTION(round,1.5,2.0); }
  void test_trunc() { TEST_UNARY_MATH_FUNCTION(trunc,1.5,1.0);
                      TEST_UNARY_MATH_FUNCTION(trunc,-1.5,-1.0); }

  void test_hypot() { TEST_BINARY_MATH_FUNCTION(hypot,2.0,1.0,2.23607); }
  void test_copysign() { TEST_BINARY_MATH_FUNCTION(copysign,3.0,-2.0,-3.0);
                         TEST_BINARY_MATH_FUNCTION(copysign,3.0,2.0,3.0); }
  void test_fdim() { TEST_BINARY_MATH_FUNCTION(fdim,3.0,2.0,1.0);
                     TEST_BINARY_MATH_FUNCTION(fdim,2.0,3.0,0.0); }

  void test_real() {
    ToyType<double> ar(1.0);
    ToyType<std::complex<double> > ac(std::complex<double>(2,3));
    TS_ASSERT_EQUALS( real(ToyType<double>(1.0))[0], 1.0 );
    TS_ASSERT_EQUALS( real(ToyType<std::complex<double> >(std::complex<double>(2,3)))[0], 2.0 );
    TS_ASSERT( is_of_type<ToyType<float> >( real(ToyType<float>()) ) );
    TS_ASSERT( is_of_type<ToyType<double> >( real(ToyType<double>()) ) );
    TS_ASSERT( is_of_type<ToyType<long double> >( real(ToyType<long double>()) ) );
    TS_ASSERT( is_of_type<ToyType<int> >( real(ToyType<int>()) ) );
    TS_ASSERT( is_of_type<ToyType<float> >( real(ToyType<std::complex<float> >()) ) );
    TS_ASSERT( is_of_type<ToyType<double> >( real(ToyType<std::complex<double> >()) ) );
    TS_ASSERT( is_of_type<ToyType<long double> >( real(ToyType<std::complex<long double> >()) ) );
    TS_ASSERT( is_of_type<ToyType<int> >( real(ToyType<std::complex<int> >()) ) );
  }

  void test_imag() {
    TS_ASSERT_EQUALS( imag(ToyType<double>(1.0))[0], 0.0 );
    TS_ASSERT_EQUALS( imag(ToyType<std::complex<double> >(std::complex<double>(2,3)))[0], 3.0 );
    TS_ASSERT( is_of_type<ToyType<float> >( imag(ToyType<float>()) ) );
    TS_ASSERT( is_of_type<ToyType<double> >( imag(ToyType<double>()) ) );
    TS_ASSERT( is_of_type<ToyType<long double> >( imag(ToyType<long double>()) ) );
    TS_ASSERT( is_of_type<ToyType<int> >( imag(ToyType<int>()) ) );
    TS_ASSERT( is_of_type<ToyType<float> >( imag(ToyType<std::complex<float> >()) ) );
    TS_ASSERT( is_of_type<ToyType<double> >( imag(ToyType<std::complex<double> >()) ) );
    TS_ASSERT( is_of_type<ToyType<long double> >( imag(ToyType<std::complex<long double> >()) ) );
    TS_ASSERT( is_of_type<ToyType<int> >( imag(ToyType<std::complex<int> >()) ) );
  }

  void test_abs() {
    TS_ASSERT_EQUALS( abs(ToyType<double>(1.0))[0], 1.0 );
    TS_ASSERT_EQUALS( abs(ToyType<double>(-1.0))[0], 1.0 );
    TS_ASSERT_DELTA( abs(ToyType<std::complex<double> >(std::complex<double>(3,4)))[0], 5.0, 0.00001 );
    TS_ASSERT( is_of_type<ToyType<float> >( abs(ToyType<float>()) ) );
    TS_ASSERT( is_of_type<ToyType<double> >( abs(ToyType<double>()) ) );
    TS_ASSERT( is_of_type<ToyType<long double> >( abs(ToyType<long double>()) ) );
    TS_ASSERT( is_of_type<ToyType<int> >( abs(ToyType<int>()) ) );
    TS_ASSERT( is_of_type<ToyType<float> >( abs(ToyType<std::complex<float> >()) ) );
    TS_ASSERT( is_of_type<ToyType<double> >( abs(ToyType<std::complex<double> >()) ) );
    TS_ASSERT( is_of_type<ToyType<long double> >( abs(ToyType<std::complex<long double> >()) ) );
  }

  void test_conj() {
    TS_ASSERT_EQUALS( conj(ToyType<double>(1.0))[0], 1.0 );
    TS_ASSERT_EQUALS( conj(ToyType<double>(-1.0))[0], -1.0 );
    TS_ASSERT_EQUALS( conj(ToyType<std::complex<double> >(std::complex<double>(3,4)))[0], std::complex<double>(3,-4) );
    TS_ASSERT( is_of_type<ToyType<float> >( conj(ToyType<float>()) ) );
    TS_ASSERT( is_of_type<ToyType<double> >( conj(ToyType<double>()) ) );
    TS_ASSERT( is_of_type<ToyType<long double> >( conj(ToyType<long double>()) ) );
    TS_ASSERT( is_of_type<ToyType<int> >( conj(ToyType<int>()) ) );
    TS_ASSERT( is_of_type<ToyType<std::complex<float> > >( conj(ToyType<std::complex<float> >()) ) );
    TS_ASSERT( is_of_type<ToyType<std::complex<double> > >( conj(ToyType<std::complex<double> >()) ) );
    TS_ASSERT( is_of_type<ToyType<std::complex<long double> > >( conj(ToyType<std::complex<long double> >()) ) );
    TS_ASSERT( is_of_type<ToyType<std::complex<int> > >( conj(ToyType<std::complex<int> >()) ) );
  }

};


