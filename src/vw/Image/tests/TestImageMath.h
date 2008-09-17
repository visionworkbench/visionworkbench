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

// TestImageMath.h
#include <cxxtest/TestSuite.h>

#include <vw/Image/ImageMath.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/Manipulation.h>

using namespace std;
using namespace vw;

class TestImageMath : public CxxTest::TestSuite
{
public:

  template <template<class> class TraitT, class T>
  static bool bool_trait( T const& arg ) {
    return TraitT<T>::value;
  }

  template <class T1, class T2>
  static bool has_pixel_type( T2 ) {
    return boost::is_same<T1,typename T2::pixel_type>::value;
  }

  void testNegation()
  {
    ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;

    ImageView<double> im2 = -im1;
    TS_ASSERT_EQUALS( im2(0,0), -1 );
    TS_ASSERT_EQUALS( im2(1,0), -2 );
    TS_ASSERT_EQUALS( im2(0,1), -3 );
    TS_ASSERT_EQUALS( im2(1,1), -4 );

    // Test the traits
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( -im1 ) );
  }

  void testSum()
  {
    ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
    ImageView<double> im2(2,2); im2(0,0)=0; im2(1,0)=3; im2(0,1)=2; im2(1,1)=9;

    ImageView<double> im3 = im1 + im2;
    TS_ASSERT_EQUALS( im3(0,0), 1 );
    TS_ASSERT_EQUALS( im3(1,0), 5 );
    TS_ASSERT_EQUALS( im3(0,1), 5 );
    TS_ASSERT_EQUALS( im3(1,1), 13 );

    im3 = im1 + 2;
    TS_ASSERT_EQUALS( im3(0,0), 3 );
    TS_ASSERT_EQUALS( im3(1,0), 4 );
    TS_ASSERT_EQUALS( im3(0,1), 5 );
    TS_ASSERT_EQUALS( im3(1,1), 6 );

    im3 = 1 + im1;
    TS_ASSERT_EQUALS( im3(0,0), 2 );
    TS_ASSERT_EQUALS( im3(1,0), 3 );
    TS_ASSERT_EQUALS( im3(0,1), 4 );
    TS_ASSERT_EQUALS( im3(1,1), 5 );

    im3 += im1;
    TS_ASSERT_EQUALS( im3(0,0), 3 );
    TS_ASSERT_EQUALS( im3(1,0), 5 );
    TS_ASSERT_EQUALS( im3(0,1), 7 );
    TS_ASSERT_EQUALS( im3(1,1), 9 );

    im3 += 2;
    TS_ASSERT_EQUALS( im3(0,0), 5 );
    TS_ASSERT_EQUALS( im3(1,0), 7 );
    TS_ASSERT_EQUALS( im3(0,1), 9 );
    TS_ASSERT_EQUALS( im3(1,1), 11 );

    (im3 += 2) += 2;
    TS_ASSERT_EQUALS( im3(0,0), 9 );
    TS_ASSERT_EQUALS( im3(1,0), 11 );
    TS_ASSERT_EQUALS( im3(0,1), 13 );
    TS_ASSERT_EQUALS( im3(1,1), 15 );

    crop(im3,0,1,2,1) -= 2;
    TS_ASSERT_EQUALS( im3(0,0), 9 );
    TS_ASSERT_EQUALS( im3(1,0), 11 );
    TS_ASSERT_EQUALS( im3(0,1), 11 );
    TS_ASSERT_EQUALS( im3(1,1), 13 );

    // Test the traits
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( im1 + im2 ) );
  }

  void testDifference()
  {
    ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
    ImageView<double> im2(2,2); im2(0,0)=0; im2(1,0)=3; im2(0,1)=2; im2(1,1)=9;

    ImageView<double> im3 = im1 - im2;
    TS_ASSERT_EQUALS( im3(0,0), 1 );
    TS_ASSERT_EQUALS( im3(1,0), -1 );
    TS_ASSERT_EQUALS( im3(0,1), 1 );
    TS_ASSERT_EQUALS( im3(1,1), -5 );

    im3 = im1 - 2;
    TS_ASSERT_EQUALS( im3(0,0), -1 );
    TS_ASSERT_EQUALS( im3(1,0), 0 );
    TS_ASSERT_EQUALS( im3(0,1), 1 );
    TS_ASSERT_EQUALS( im3(1,1), 2 );

    im3 = 1 - im1;
    TS_ASSERT_EQUALS( im3(0,0), 0 );
    TS_ASSERT_EQUALS( im3(1,0), -1 );
    TS_ASSERT_EQUALS( im3(0,1), -2 );
    TS_ASSERT_EQUALS( im3(1,1), -3 );

    im3 -= im1;
    TS_ASSERT_EQUALS( im3(0,0), -1 );
    TS_ASSERT_EQUALS( im3(1,0), -3 );
    TS_ASSERT_EQUALS( im3(0,1), -5 );
    TS_ASSERT_EQUALS( im3(1,1), -7 );

    im3 -= 2;
    TS_ASSERT_EQUALS( im3(0,0), -3 );
    TS_ASSERT_EQUALS( im3(1,0), -5 );
    TS_ASSERT_EQUALS( im3(0,1), -7 );
    TS_ASSERT_EQUALS( im3(1,1), -9 );

    (im3 -= 2 ) -= 2;
    TS_ASSERT_EQUALS( im3(0,0), -7 );
    TS_ASSERT_EQUALS( im3(1,0), -9 );
    TS_ASSERT_EQUALS( im3(0,1), -11 );
    TS_ASSERT_EQUALS( im3(1,1), -13 );

    // Test the traits
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( im1 - im2 ) );
  }

  void testProduct()
  {
    ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
    ImageView<double> im2(2,2); im2(0,0)=0; im2(1,0)=3; im2(0,1)=2; im2(1,1)=9;

    ImageView<double> im3 = im1 * im2;
    TS_ASSERT_EQUALS( im3(0,0), 0 );
    TS_ASSERT_EQUALS( im3(1,0), 6 );
    TS_ASSERT_EQUALS( im3(0,1), 6 );
    TS_ASSERT_EQUALS( im3(1,1), 36 );

    im3 = im1 * 2;
    TS_ASSERT_EQUALS( im3(0,0), 2 );
    TS_ASSERT_EQUALS( im3(1,0), 4 );
    TS_ASSERT_EQUALS( im3(0,1), 6 );
    TS_ASSERT_EQUALS( im3(1,1), 8 );

    im3 = 3 * im1;
    TS_ASSERT_EQUALS( im3(0,0), 3 );
    TS_ASSERT_EQUALS( im3(1,0), 6 );
    TS_ASSERT_EQUALS( im3(0,1), 9 );
    TS_ASSERT_EQUALS( im3(1,1), 12 );

    im3 *= im1;
    TS_ASSERT_EQUALS( im3(0,0), 3 );
    TS_ASSERT_EQUALS( im3(1,0), 12 );
    TS_ASSERT_EQUALS( im3(0,1), 27 );
    TS_ASSERT_EQUALS( im3(1,1), 48 );

    im3 *= 0.5;
    TS_ASSERT_EQUALS( im3(0,0), 1.5 );
    TS_ASSERT_EQUALS( im3(1,0), 6 );
    TS_ASSERT_EQUALS( im3(0,1), 13.5 );
    TS_ASSERT_EQUALS( im3(1,1), 24 );

    (im3 *= 2) *= 2;

    TS_ASSERT_EQUALS( im3(0,0), 6 );
    TS_ASSERT_EQUALS( im3(1,0), 24 );
    TS_ASSERT_EQUALS( im3(0,1), 54 );
    TS_ASSERT_EQUALS( im3(1,1), 96 );

    // Test the traits
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( im1 * im2 ) );
  }

  void testQuotient()
  {
    ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
    ImageView<double> im2(2,2); im2(0,0)=1; im2(1,0)=4; im2(0,1)=2; im2(1,1)=8;

    ImageView<double> im3 = im1 / im2;
    TS_ASSERT_EQUALS( im3(0,0), 1 );
    TS_ASSERT_EQUALS( im3(1,0), 0.5 );
    TS_ASSERT_EQUALS( im3(0,1), 1.5 );
    TS_ASSERT_EQUALS( im3(1,1), 0.5 );

    im3 = im1 / 2;
    TS_ASSERT_EQUALS( im3(0,0), 0.5 );
    TS_ASSERT_EQUALS( im3(1,0), 1 );
    TS_ASSERT_EQUALS( im3(0,1), 1.5 );
    TS_ASSERT_EQUALS( im3(1,1), 2 );

    im3 = 2 / im2;
    TS_ASSERT_EQUALS( im3(0,0), 2 );
    TS_ASSERT_EQUALS( im3(1,0), 0.5 );
    TS_ASSERT_EQUALS( im3(0,1), 1 );
    TS_ASSERT_EQUALS( im3(1,1), 0.25 );

    im3 /= im2;
    TS_ASSERT_EQUALS( im3(0,0), 2 );
    TS_ASSERT_EQUALS( im3(1,0), 0.125 );
    TS_ASSERT_EQUALS( im3(0,1), 0.5 );
    TS_ASSERT_EQUALS( im3(1,1), 0.03125 );

    im3 /= 0.5;
    TS_ASSERT_EQUALS( im3(0,0), 4 );
    TS_ASSERT_EQUALS( im3(1,0), 0.25 );
    TS_ASSERT_EQUALS( im3(0,1), 1 );
    TS_ASSERT_EQUALS( im3(1,1), 0.0625 );

    (im3 /= 0.5) /= 0.5;
    TS_ASSERT_EQUALS( im3(0,0), 16 );
    TS_ASSERT_EQUALS( im3(1,0), 1 );
    TS_ASSERT_EQUALS( im3(0,1), 4 );
    TS_ASSERT_EQUALS( im3(1,1), 0.25 );

    // Test the traits
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( im1 / im2 ) );
  }

  void test_real() {
    ImageView<std::complex<double> > im(1,1); im(0,0)=std::complex<double>(1.0,2.0);
    TS_ASSERT_EQUALS( real(im).cols(), im.cols()  );
    TS_ASSERT_EQUALS( real(im).rows(), im.rows()  );
    TS_ASSERT_EQUALS( real(im).planes(), im.planes()  );
    TS_ASSERT_EQUALS( real(im)(0,0), 1.0 );
    TS_ASSERT( has_pixel_type<float>( real(ImageView<float>()) ) );
    TS_ASSERT( has_pixel_type<double>( real(ImageView<double>()) ) );
    TS_ASSERT( has_pixel_type<long double>( real(ImageView<long double>()) ) );
    TS_ASSERT( has_pixel_type<int>( real(ImageView<int>()) ) );
    TS_ASSERT( has_pixel_type<float>( real(ImageView<std::complex<float> >()) ) );
    TS_ASSERT( has_pixel_type<double>( real(ImageView<std::complex<double> >()) ) );
    TS_ASSERT( has_pixel_type<long double>( real(ImageView<std::complex<long double> >()) ) );
    TS_ASSERT( has_pixel_type<int>( real(ImageView<std::complex<int> >()) ) );
  }

  void test_imag() {
    ImageView<std::complex<double> > im(1,1);
    im(0,0)=std::complex<double>(1.0,2.0);
    TS_ASSERT_EQUALS( imag(im).cols(), im.cols()  );
    TS_ASSERT_EQUALS( imag(im).rows(), im.rows()  );
    TS_ASSERT_EQUALS( imag(im).planes(), im.planes()  );
    TS_ASSERT_EQUALS( imag(im)(0,0), 2.0 );
    TS_ASSERT( has_pixel_type<float>( imag(ImageView<float>()) ) );
    TS_ASSERT( has_pixel_type<double>( imag(ImageView<double>()) ) );
    TS_ASSERT( has_pixel_type<long double>( imag(ImageView<long double>()) ) );
    TS_ASSERT( has_pixel_type<int>( imag(ImageView<int>()) ) );
    TS_ASSERT( has_pixel_type<float>( imag(ImageView<std::complex<float> >()) ) );
    TS_ASSERT( has_pixel_type<double>( imag(ImageView<std::complex<double> >()) ) );
    TS_ASSERT( has_pixel_type<long double>( imag(ImageView<std::complex<long double> >()) ) );
    TS_ASSERT( has_pixel_type<int>( imag(ImageView<std::complex<int> >()) ) );
  }

  void test_abs() {
    ImageView<std::complex<double> > im(1,1); im(0,0)=std::complex<double>(1.0,2.0);
    TS_ASSERT_EQUALS( abs(im).cols(), im.cols()  );
    TS_ASSERT_EQUALS( abs(im).rows(), im.rows()  );
    TS_ASSERT_EQUALS( abs(im).planes(), im.planes()  );
    TS_ASSERT_DELTA( abs(im)(0,0), 2.23607, 0.00001 );
    TS_ASSERT( has_pixel_type<float>( abs(ImageView<float>()) ) );
    TS_ASSERT( has_pixel_type<double>( abs(ImageView<double>()) ) );
    TS_ASSERT( has_pixel_type<long double>( abs(ImageView<long double>()) ) );
    TS_ASSERT( has_pixel_type<int>( abs(ImageView<int>()) ) );
    TS_ASSERT( has_pixel_type<float>( abs(ImageView<std::complex<float> >()) ) );
    TS_ASSERT( has_pixel_type<double>( abs(ImageView<std::complex<double> >()) ) );
    TS_ASSERT( has_pixel_type<long double>( abs(ImageView<std::complex<long double> >()) ) );
  }

  void test_conj() {
    ImageView<std::complex<double> > im(1,1); im(0,0)=std::complex<double>(1.0,2.0);
    TS_ASSERT_EQUALS( conj(im).cols(), im.cols()  );
    TS_ASSERT_EQUALS( conj(im).rows(), im.rows()  );
    TS_ASSERT_EQUALS( conj(im).planes(), im.planes()  );
    TS_ASSERT_EQUALS( conj(im)(0,0).real(), 1.0 );
    TS_ASSERT_EQUALS( conj(im)(0,0).imag(), -2.0 );
    TS_ASSERT( has_pixel_type<float>( conj(ImageView<float>()) ) );
    TS_ASSERT( has_pixel_type<double>( conj(ImageView<double>()) ) );
    TS_ASSERT( has_pixel_type<long double>( conj(ImageView<long double>()) ) );
    TS_ASSERT( has_pixel_type<int>( conj(ImageView<int>()) ) );
    TS_ASSERT( has_pixel_type<std::complex<float> >( conj(ImageView<std::complex<float> >()) ) );
    TS_ASSERT( has_pixel_type<std::complex<double> >( conj(ImageView<std::complex<double> >()) ) );
    TS_ASSERT( has_pixel_type<std::complex< long double> >( conj(ImageView<std::complex<long double> >()) ) );
    TS_ASSERT( has_pixel_type<std::complex<int> >( conj(ImageView<std::complex<int> >()) ) );
  }

#define TEST_UNARY_MATH_FUNCTION(func,arg,result)                       \
  do {                                                                  \
    ImageView<double> im(1,1); im(0,0)=(arg);                           \
    TS_ASSERT_EQUALS( func(im).cols(), im.cols() );                     \
    TS_ASSERT_EQUALS( func(im).rows(), im.rows() );                     \
    TS_ASSERT_EQUALS( func(im).planes(), im.planes() );                 \
    TS_ASSERT_DELTA( func(im)(0,0), result, 0.00001 );                  \
    TS_ASSERT( has_pixel_type<float>( func(ImageView<float>()) ) );     \
    TS_ASSERT( has_pixel_type<double>( func(ImageView<double>()) ) );   \
    TS_ASSERT( has_pixel_type<long double>( func(ImageView<long double>()) ) ); \
    TS_ASSERT( has_pixel_type<double>( func(ImageView<int>()) ) );      \
  } while(false)

#define TEST_BINARY_MATH_FUNCTION(func,arg1,arg2,result)                \
  do {                                                                  \
    ImageView<double> im1(1,1); im1(0,0)=(arg1);                        \
    ImageView<double> im2(1,1); im2(0,0)=(arg2);                        \
    TS_ASSERT_EQUALS( func(im1,im2).cols(), im1.cols() );               \
    TS_ASSERT_EQUALS( func(im1,im2).rows(), im1.rows() );               \
    TS_ASSERT_EQUALS( func(im1,im2).planes(), im1.planes() );           \
    TS_ASSERT_DELTA( func(im1,im2)(0,0), result, 0.00001 );             \
    TS_ASSERT_EQUALS( func(im1,arg2).cols(), im1.cols() );              \
    TS_ASSERT_EQUALS( func(im1,arg2).rows(), im1.rows() );              \
    TS_ASSERT_EQUALS( func(im1,arg2).planes(), im1.planes() );          \
    TS_ASSERT_DELTA( func(im1,arg2)(0,0), result, 0.00001 );            \
    TS_ASSERT_EQUALS( func(arg1,im2).cols(), im1.cols() );              \
    TS_ASSERT_EQUALS( func(arg1,im2).rows(), im1.rows() );              \
    TS_ASSERT_EQUALS( func(arg1,im2).planes(), im1.planes() );          \
    TS_ASSERT_DELTA( func(arg1,im2)(0,0), result, 0.00001 );            \
    TS_ASSERT( has_pixel_type<float>( func(ImageView<float>(),ImageView<float>()) ) ); \
    TS_ASSERT( has_pixel_type<double>( func(ImageView<double>(),ImageView<double>()) ) ); \
    TS_ASSERT( has_pixel_type<long double>( func(ImageView<long double>(),ImageView<long double>()) ) ); \
    TS_ASSERT( has_pixel_type<double>( func(ImageView<int>(),ImageView<int>()) ) );  \
    TS_ASSERT( has_pixel_type<float>( func(ImageView<int>(),ImageView<float>()) ) ); \
    TS_ASSERT( has_pixel_type<double>( func(ImageView<float>(),ImageView<double>()) ) ); \
    TS_ASSERT( has_pixel_type<long double>( func(ImageView<double>(),ImageView<long double>()) ) ); \
    TS_ASSERT( has_pixel_type<double>( func(ImageView<char>(),ImageView<int>()) ) );  \
    TS_ASSERT( has_pixel_type<float>( func(ImageView<float>(),float()) ) ); \
    TS_ASSERT( has_pixel_type<double>( func(ImageView<double>(),double()) ) ); \
    TS_ASSERT( has_pixel_type<double>( func(ImageView<int>(),int()) ) );  \
    TS_ASSERT( has_pixel_type<float>( func(ImageView<int>(),float()) ) ); \
    TS_ASSERT( has_pixel_type<double>( func(ImageView<float>(),double()) ) ); \
    TS_ASSERT( has_pixel_type<long double>( func(ImageView<long double>(),double()) ) ); \
    TS_ASSERT( has_pixel_type<double>( func(ImageView<char>(),int()) ) );  \
    TS_ASSERT( has_pixel_type<float>( func(float(),ImageView<float>()) ) ); \
    TS_ASSERT( has_pixel_type<double>( func(double(),ImageView<double>()) ) ); \
    TS_ASSERT( has_pixel_type<double>( func(int(),ImageView<int>()) ) );  \
    TS_ASSERT( has_pixel_type<float>( func(int(),ImageView<float>()) ) ); \
    TS_ASSERT( has_pixel_type<double>( func(float(),ImageView<double>()) ) ); \
    TS_ASSERT( has_pixel_type<long double>( func(double(),ImageView<long double>()) ) ); \
    TS_ASSERT( has_pixel_type<double>( func(char(),ImageView<int>()) ) );  \
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

#ifndef WIN32
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
#endif

}; // class TestImageMath
