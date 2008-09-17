// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
//
// Copyright 2006 Carnegie Mellon University. All rights reserved.
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

// Math/TestFunctors.h
#include <cxxtest/TestSuite.h>

#include <vw/Math/Functors.h>

using namespace vw;
using namespace vw::math;

class TestFunctors : public CxxTest::TestSuite
{
public:

  template <class T1, class T2>
  static bool is_of_type( T2 ) {
    return boost::is_same<T1,T2>::value;
  }

  void test_real() {
    ArgRealFunctor f;
    TS_ASSERT_EQUALS( f(2.0), 2.0 );
    TS_ASSERT_EQUALS( f(std::complex<double>(2.0,3.0)), 2.0 );
    TS_ASSERT( is_of_type<float>( f(float()) ) );
    TS_ASSERT( is_of_type<double>( f(double()) ) );
    TS_ASSERT( is_of_type<long double>( f((long double)(0)) ) );
    TS_ASSERT( is_of_type<int>( f(int()) ) );
    TS_ASSERT( is_of_type<float>( f(std::complex<float> ()) ) );
    TS_ASSERT( is_of_type<double>( f(std::complex<double> ()) ) );
    TS_ASSERT( is_of_type<long double>( f(std::complex<long double> ()) ) );
    TS_ASSERT( is_of_type<int>( f(std::complex<int> ()) ) );
  }

  void test_imag() {
    ArgImagFunctor f;
    TS_ASSERT_EQUALS( f(2.0), 0.0 );
    TS_ASSERT_EQUALS( f(std::complex<double>(2.0,3.0)), 3.0 );
    TS_ASSERT( is_of_type<float>( f(float()) ) );
    TS_ASSERT( is_of_type<double>( f(double()) ) );
    TS_ASSERT( is_of_type<long double>( f((long double)(0)) ) );
    TS_ASSERT( is_of_type<int>( f(int()) ) );
    TS_ASSERT( is_of_type<float>( f(std::complex<float>()) ) );
    TS_ASSERT( is_of_type<double>( f(std::complex<double>()) ) );
    TS_ASSERT( is_of_type<long double>( f(std::complex<long double>()) ) );
    TS_ASSERT( is_of_type<int>( f(std::complex<int>()) ) );
  }

  void test_abs() {
    ArgAbsFunctor f;
    TS_ASSERT_EQUALS( f(2.0), 2.0 );
    TS_ASSERT_EQUALS( f(-2.0), 2.0 );
    TS_ASSERT_DELTA( f(std::complex<double>(3.0,4.0)), 5.0, 0.00001 );
    TS_ASSERT( is_of_type<float>( f(float()) ) );
    TS_ASSERT( is_of_type<double>( f(double()) ) );
    TS_ASSERT( is_of_type<long double>( f((long double)(0)) ) );
    TS_ASSERT( is_of_type<int>( f(int()) ) );
    TS_ASSERT( is_of_type<float>( f(std::complex<float>()) ) );
    TS_ASSERT( is_of_type<double>( f(std::complex<double>()) ) );
    TS_ASSERT( is_of_type<long double>( f(std::complex<long double>()) ) );
  }

  void test_conj() {
    ArgConjFunctor f;
    TS_ASSERT_EQUALS( f(2.0), 2.0 );
    TS_ASSERT_EQUALS( f(std::complex<double>(1.0,2.0)), std::complex<double>(1.0,-2.0) );
    TS_ASSERT( is_of_type<float>( f(float()) ) );
    TS_ASSERT( is_of_type<double>( f(double()) ) );
    TS_ASSERT( is_of_type<long double>( f((long double)(0)) ) );
    TS_ASSERT( is_of_type<int>( f(int()) ) );
    TS_ASSERT( is_of_type<std::complex<float> >( f(std::complex<float>()) ) );
    TS_ASSERT( is_of_type<std::complex<double> >( f(std::complex<double>()) ) );
    TS_ASSERT( is_of_type<std::complex<long double> >( f(std::complex<long double>()) ) );
    TS_ASSERT( is_of_type<std::complex<int> >( f(std::complex<int>()) ) );
  }

#define TEST_UNARY_MATH_FUNCTOR(func,arg,result)                        \
  do {                                                                  \
    Arg##func##Functor f;                                               \
    TS_ASSERT_DELTA( f((float)(arg)), (result), 0.00001 );              \
    TS_ASSERT_DELTA( f((double)(arg)), (result), 0.00001 );             \
    TS_ASSERT_DELTA( f((long double)(arg)), (result), 0.00001 );        \
    TS_ASSERT( is_of_type<float>(f((float)(arg))) );                    \
    TS_ASSERT( is_of_type<double>(f((double)(arg))) );                  \
    TS_ASSERT( is_of_type<long double>(f((long double)(arg))) );        \
    TS_ASSERT( is_of_type<double>(f((int)(arg))) );                     \
  } while(false)

#define TEST_BINARY_MATH_FUNCTOR(func,arg1,arg2,result)                 \
  do {                                                                  \
    ArgArg##func##Functor f;                                            \
    TS_ASSERT_DELTA( f((float)(arg1),(float)(arg2)), (result), 0.00001 ); \
    TS_ASSERT_DELTA( f((double)(arg1),(double)(arg2)), (result), 0.00001 ); \
    TS_ASSERT_DELTA( f((long double)(arg1),(long double)(arg2)), (result), 0.00001 ); \
    TS_ASSERT( is_of_type<float>(f((float)(arg1),(float)(arg2))) );     \
    TS_ASSERT( is_of_type<double>(f((double)(arg1),(double)(arg2))) );  \
    TS_ASSERT( is_of_type<long double>(f((long double)(arg1),(long double)(arg2))) ); \
    TS_ASSERT( is_of_type<double>( f((float)(arg1),(double)(arg2))) );      \
    TS_ASSERT( is_of_type<long double>( f((float)(arg1),(long double)(arg2))) ); \
    TS_ASSERT( is_of_type<double>( f((int)(arg1),(int)(arg2))) );             \
    TS_ASSERT( is_of_type<float>( f((int)(arg1),(float)(arg2))) );            \
    TS_ASSERT_DELTA( ArgVal##func##Functor<float>(arg2)((float)(arg1)), (result), 0.00001 ); \
    TS_ASSERT_DELTA( ArgVal##func##Functor<double>(arg2)((double)(arg1)), (result), 0.00001 ); \
    TS_ASSERT_DELTA( ArgVal##func##Functor<long double>(arg2)((long double)(arg1)), (result), 0.00001 ); \
    TS_ASSERT( is_of_type<float>(ArgVal##func##Functor<float>(arg2)((float)(arg1))) );     \
    TS_ASSERT( is_of_type<double>(ArgVal##func##Functor<double>(arg2)((double)(arg1))) );  \
    TS_ASSERT( is_of_type<long double>(ArgVal##func##Functor<long double>(arg2)((long double)(arg1))) ); \
    TS_ASSERT( is_of_type<double>( ArgVal##func##Functor<double>(arg2)((float)(arg1))) ); \
    TS_ASSERT( is_of_type<long double>( ArgVal##func##Functor<long double>(arg2)((float)(arg1))) ); \
    TS_ASSERT( is_of_type<double>( ArgVal##func##Functor<int>((int)(arg2))((int)(arg1))) ); \
    TS_ASSERT( is_of_type<float>( ArgVal##func##Functor<float>(arg2)((int)(arg1))) );            \
    TS_ASSERT_DELTA( ValArg##func##Functor<float>(arg1)((float)(arg2)), (result), 0.00001 ); \
    TS_ASSERT_DELTA( ValArg##func##Functor<double>(arg1)((double)(arg2)), (result), 0.00001 ); \
    TS_ASSERT_DELTA( ValArg##func##Functor<long double>(arg1)((long double)(arg2)), (result), 0.00001 ); \
    TS_ASSERT( is_of_type<float>(ValArg##func##Functor<float>(arg1)((float)(arg2))) );     \
    TS_ASSERT( is_of_type<double>(ValArg##func##Functor<double>(arg1)((double)(arg2))) );  \
    TS_ASSERT( is_of_type<long double>(ValArg##func##Functor<long double>(arg1)((long double)(arg2))) ); \
    TS_ASSERT( is_of_type<double>( ValArg##func##Functor<float>(arg1)((double)(arg2))) ); \
    TS_ASSERT( is_of_type<long double>( ValArg##func##Functor<float>(arg1)((long double)(arg2))) ); \
    TS_ASSERT( is_of_type<double>( ValArg##func##Functor<int>((int)(arg1))((int)(arg2))) ); \
    TS_ASSERT( is_of_type<float>( ValArg##func##Functor<int>((int)arg1)((float)(arg2))) ); \
  } while(false)

  void testAcosFunctor() { TEST_UNARY_MATH_FUNCTOR(Acos,0.5,1.0472); }
  void testAsinFunctor() { TEST_UNARY_MATH_FUNCTOR(Asin,0.5,0.523599); }
  void testAtanFunctor() { TEST_UNARY_MATH_FUNCTOR(Atan,1.0,0.785398); }
  void testCosFunctor() { TEST_UNARY_MATH_FUNCTOR(Cos,1.0,0.540302); }
  void testSinFunctor() { TEST_UNARY_MATH_FUNCTOR(Sin,1.0,0.841471); }
  void testTanFunctor() { TEST_UNARY_MATH_FUNCTOR(Tan,1.0,1.55741); }
  void testCoshFunctor() { TEST_UNARY_MATH_FUNCTOR(Cosh,1.0,1.54308); }
  void testSinhFunctor() { TEST_UNARY_MATH_FUNCTOR(Sinh,1.0,1.1752); }
  void testTanhFunctor() { TEST_UNARY_MATH_FUNCTOR(Tanh,1.0,0.761594); }
  void testExpFunctor() { TEST_UNARY_MATH_FUNCTOR(Exp,1.0,2.718281); }
  void testLogFunctor() { TEST_UNARY_MATH_FUNCTOR(Log,2.0,0.693147); }
  void testLog10Functor() { TEST_UNARY_MATH_FUNCTOR(Log10,2.0,0.30103); }
  void testSqrtFunctor() { TEST_UNARY_MATH_FUNCTOR(Sqrt,2.0,1.41421); }
  void testCeilFunctor() { TEST_UNARY_MATH_FUNCTOR(Ceil,1.5,2.0);
                           TEST_UNARY_MATH_FUNCTOR(Ceil,-1.5,-1.0); }
  void testFloorFunctor() { TEST_UNARY_MATH_FUNCTOR(Floor,1.5,1.0);
                            TEST_UNARY_MATH_FUNCTOR(Floor,-1.5,-2.0); }

  void testAtan2Functor() { TEST_BINARY_MATH_FUNCTOR(Atan2,2.0,1.0,1.10715); }
  void testPowFunctor() { TEST_BINARY_MATH_FUNCTOR(Pow,3.0,2.0,9.0); }

#ifndef WIN32
  void testAcoshFunctor() { TEST_UNARY_MATH_FUNCTOR(Acosh,1.5,0.962424); }
  void testAsinhFunctor() { TEST_UNARY_MATH_FUNCTOR(Asinh,1.0,0.881374); }
  void testAtanhFunctor() { TEST_UNARY_MATH_FUNCTOR(Atanh,0.5,0.549306); }
  void testExp2Functor() {
#ifndef __FreeBSD__
    TEST_UNARY_MATH_FUNCTOR(Exp2,1.0,2.0);
#endif
  }
  void testLog2Functor() {
#ifndef __FreeBSD__
    TEST_UNARY_MATH_FUNCTOR(Log2,2.0,1.0);
#endif
  }
  void testTgammaFunctor() {
#ifndef __FreeBSD__
    TEST_UNARY_MATH_FUNCTOR(Tgamma,1.5,0.886227);
#endif
  }
  void testExpm1Functor() { TEST_UNARY_MATH_FUNCTOR(Expm1,1.0,1.718281); }
  void testLog1pFunctor() { TEST_UNARY_MATH_FUNCTOR(Log1p,1.0,0.693147); }
  void testCbrtFunctor() { TEST_UNARY_MATH_FUNCTOR(Cbrt,2.0,1.25992); }
  void testErfFunctor() { TEST_UNARY_MATH_FUNCTOR(Erf,1.0,0.842701); }
  void testErfcFunctor() { TEST_UNARY_MATH_FUNCTOR(Erfc,1.0,0.157299); }
  void testLgammaFunctor() { TEST_UNARY_MATH_FUNCTOR(Lgamma,2.5,0.284683); }
  void testRoundFunctor() { TEST_UNARY_MATH_FUNCTOR(Round,1.4,1.0);
                            TEST_UNARY_MATH_FUNCTOR(Round,1.5,2.0); }
  void testTruncFunctor() { TEST_UNARY_MATH_FUNCTOR(Trunc,1.5,1.0);
                            TEST_UNARY_MATH_FUNCTOR(Trunc,-1.5,-1.0); }

  void testHypotFunctor() { TEST_BINARY_MATH_FUNCTOR(Hypot,2.0,1.0,2.23607); }
  void testCopysignFunctor() { TEST_BINARY_MATH_FUNCTOR(Copysign,3.0,-2.0,-3.0);
                               TEST_BINARY_MATH_FUNCTOR(Copysign,3.0,2.0,3.0); }
  void testFdimFunctor() { TEST_BINARY_MATH_FUNCTOR(Fdim,3.0,2.0,1.0);
                           TEST_BINARY_MATH_FUNCTOR(Fdim,2.0,3.0,0.0); }
#endif
};
