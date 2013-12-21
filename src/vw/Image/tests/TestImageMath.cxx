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


// TestImageMath.h
#include <gtest/gtest_VW.h>

#include <vw/config.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/PixelMask.h>

#include <boost/type_traits/is_same.hpp>

using namespace vw;

template <template<class> class TraitT, class T>
static bool bool_trait( T const& /*arg*/ ) {
  return TraitT<T>::value;
}

template <class T1, class T2>
static bool has_pixel_type( T2 ) {
  return boost::is_same<T1,typename T2::pixel_type>::value;
}

// BINARY OPERATION TEST
template <class PixelT>
class ImageBinaryMathTest : public ::testing::Test {
protected:
  ImageBinaryMathTest() {}

  virtual void SetUp() {
    image_a.set_size(2,2);
    image_a(0,0)=PixelT(1); image_a(1,0)=PixelT(2);
    image_a(0,1)=PixelT(3); image_a(1,1)=PixelT(4);
    image_b.set_size(2,2);
    image_b(0,0)=PixelT(0); image_b(1,0)=PixelT(3);
    image_b(0,1)=PixelT(2); image_b(1,1)=PixelT(9);
  }

  ImageView<PixelT> image_a, image_b;
};

#define IMAGE_EXPECT_EQ( image, a, b, c, d ) \
  EXPECT_EQ( image(0,0), a );  \
  EXPECT_EQ( image(1,0), b );  \
  EXPECT_EQ( image(0,1), c );  \
  EXPECT_EQ( image(1,1), d );

TEST( ImageMath, Negation ) {
  ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;

  ImageView<double> im2 = -im1;
  IMAGE_EXPECT_EQ( im2, -1, -2, -3, -4 );

  // Test the traits
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( -im1 ) );
}

typedef ImageBinaryMathTest<double> BinaryDoubleMath;
TEST_F( BinaryDoubleMath, Sum ) {
  ImageView<double> im3 = image_a + image_b;
  IMAGE_EXPECT_EQ( im3, 1, 5, 5, 13 );
  im3 = image_a + 2;
  IMAGE_EXPECT_EQ( im3, 3, 4, 5, 6 );
  im3 = 1 + image_a;
  IMAGE_EXPECT_EQ( im3, 2, 3, 4, 5 );
  im3 += image_a;
  IMAGE_EXPECT_EQ( im3, 3, 5, 7, 9 );
  im3 += 2;
  IMAGE_EXPECT_EQ( im3, 5, 7, 9, 11 );
  (im3 += 2 ) += 2;
  IMAGE_EXPECT_EQ( im3, 9, 11, 13, 15 );
  crop(im3,0,1,2,1) += 2;
  IMAGE_EXPECT_EQ( im3, 9, 11, 15, 17 );

  // Test the traits
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( image_a + image_b ) );
}


TEST_F( BinaryDoubleMath, Difference ) {
  ImageView<double> im3 = image_a - image_b;
  IMAGE_EXPECT_EQ( im3, 1, -1, 1, -5 );
  im3 = image_a - 2;
  IMAGE_EXPECT_EQ( im3, -1, 0, 1, 2 );
  im3 = 1 - image_a;
  IMAGE_EXPECT_EQ( im3, 0, -1, -2, -3 );
  im3 -= image_a;
  IMAGE_EXPECT_EQ( im3, -1, -3, -5, -7 );
  im3 -= 2;
  IMAGE_EXPECT_EQ( im3, -3, -5, -7, -9 );
  (im3 -= 2 ) -= 2;
  IMAGE_EXPECT_EQ( im3, -7, -9, -11, -13 );

  // Test the traits
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( image_a - image_b ) );
}

TEST_F( BinaryDoubleMath, Product ) {
  ImageView<double> im3 = image_a * image_b;
  IMAGE_EXPECT_EQ( im3, 0, 6, 6, 36 );
  im3 = image_a * 2;
  IMAGE_EXPECT_EQ( im3, 2, 4, 6, 8 );
  im3 = 3 * image_a;
  IMAGE_EXPECT_EQ( im3, 3, 6, 9, 12 );
  im3 *= image_a;
  IMAGE_EXPECT_EQ( im3, 3, 12, 27, 48 );
  im3 *= 0.5;
  IMAGE_EXPECT_EQ( im3, 1.5, 6, 13.5, 24 );
  (im3 *= 2 ) *= 2;
  IMAGE_EXPECT_EQ( im3, 6, 24, 54, 96 );

  // Test the traits
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( image_a * image_b ) );
}

TEST_F( BinaryDoubleMath, Quotient ) {
  ImageView<double> im3 = image_a / image_b;
  IMAGE_EXPECT_EQ( im3, 0, (2.0/3.0), 1.5, (4.0/9.0) );
  im3 = image_a / 2;
  IMAGE_EXPECT_EQ( im3, 0.5, 1, 1.5, 2 );
  im3 = 2 / image_b;
  IMAGE_EXPECT_EQ( im3, 0, (2.0/3.0), 1, (2.0/9.0) );
  im3 /= image_b;
  IMAGE_EXPECT_EQ( im3, 0, (2.0/9.0), 0.5, (2.0/81) );
  im3 /= 0.5;
  IMAGE_EXPECT_EQ( im3, 0, (4.0/9.0), 1, (4.0/81) );
  (im3 /= 1.0/9 ) /= 1.0/9;
  IMAGE_EXPECT_EQ( im3, 0, 36, 81, 4 );

  // Test the traits
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( image_a / image_b ) );
}

typedef ImageBinaryMathTest<PixelMask<double> > BinaryMaskedDoubleMath;
TEST_F( BinaryMaskedDoubleMath, Sum ) {

  ImageView<PixelMask<double> > im3 = image_a + image_b;
  IMAGE_EXPECT_EQ( im3, 1, 5, 5, 13 );
  im3 = image_a + 2;
  IMAGE_EXPECT_EQ( im3, 3, 4, 5, 6 );
  im3 = 1 + image_a;
  IMAGE_EXPECT_EQ( im3, 2, 3, 4, 5 );
  im3 += image_a;
  IMAGE_EXPECT_EQ( im3, 3, 5, 7, 9 );
  im3 += 2;
  IMAGE_EXPECT_EQ( im3, 5, 7, 9, 11 );
  ( im3 += 2 ) -= 3;
  IMAGE_EXPECT_EQ( im3, 4, 6, 8, 10 );

  ASSERT_TRUE( is_valid( im3(0,0) ) );

  im3 += ImageView<PixelMask<double> >(2,2); // All pixel should be invalid.
  ASSERT_FALSE( is_valid( im3(0,0) ) );
}

TEST_F( BinaryMaskedDoubleMath, Quotient ) {
  ImageView<PixelMask<double> > im3 = image_a / image_b;
  IMAGE_EXPECT_EQ( im3, 0, 2.0/3.0, 1.5, 4.0/9.0 );
  im3 = image_a / 2;
  IMAGE_EXPECT_EQ( im3, 0.5, 1, 1.5, 2 );
  im3 = 2 / image_b;
  IMAGE_EXPECT_EQ( im3, 0, 2.0/3.0, 1, 2.0/9.0 );
  im3 /= image_b;
  IMAGE_EXPECT_EQ( im3, 0, 2.0/9.0, 0.5, 2.0/81 );
  im3 /= 0.5;
  IMAGE_EXPECT_EQ( im3, 0, 4.0/9.0, 1, 4.0/81 );
  (im3 /= 1.0/9 ) /= 1.0/9;
  IMAGE_EXPECT_EQ( im3, 0, 36, 81, 4 );

  // Test the traits
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( image_a / image_b ) );
}

#define ASSERT_PRESERVED_TYPE( op ) \
  ASSERT_TRUE( has_pixel_type<float>( op(ImageView<float>()) ) ); \
  ASSERT_TRUE( has_pixel_type<double>( op(ImageView<double>()) ) ); \
  ASSERT_TRUE( has_pixel_type<long double>( op(ImageView<long double>()) ) ); \
  ASSERT_TRUE( has_pixel_type<int>( op(ImageView<int>()) ) ); \
  ASSERT_TRUE( has_pixel_type<float>( op(ImageView<std::complex<float> >()) ) ); \
  ASSERT_TRUE( has_pixel_type<double>( op(ImageView<std::complex<double> >()) ) ); \
  ASSERT_TRUE( has_pixel_type<long double>( op(ImageView<std::complex<long double> >()) ) ); \
  ASSERT_TRUE( has_pixel_type<int>( op(ImageView<std::complex<int> >()) ) ); \
  ASSERT_TRUE( has_pixel_type<PixelGray<float> >( op(ImageView<PixelGray<float> >()) ) ); \
  ASSERT_TRUE( has_pixel_type<PixelGray<int> >( op(ImageView<PixelGray<int> >()) ) ); \
  ASSERT_TRUE( has_pixel_type<PixelGray<float> >( op(ImageView<PixelGray<std::complex<float> > >()) ) ); \
  ASSERT_TRUE( has_pixel_type<PixelRGB<float> >( op(ImageView<PixelRGB<float> >()) ) ); \
  ASSERT_TRUE( has_pixel_type<PixelRGB<int> >( op(ImageView<PixelRGB<int> >()) ) ); \
  ASSERT_TRUE( has_pixel_type<PixelRGB<float> >( op(ImageView<PixelRGB<std::complex<float> > >()) ) ); \

TEST( ImageMath, Real ) {
  ImageView<std::complex<double> > im(1,1); im(0,0)=std::complex<double>(1.0,2.0);
  EXPECT_EQ( real(im).cols(), im.cols()  );
  EXPECT_EQ( real(im).rows(), im.rows()  );
  EXPECT_EQ( real(im).planes(), im.planes()  );
  EXPECT_EQ( real(im)(0,0), 1.0 );
  ASSERT_PRESERVED_TYPE( real );
}

TEST( ImageMath, Imag ) {
  ImageView<std::complex<double> > im(1,1);
  im(0,0)=std::complex<double>(1.0,2.0);
  EXPECT_EQ( imag(im).cols(), im.cols()  );
  EXPECT_EQ( imag(im).rows(), im.rows()  );
  EXPECT_EQ( imag(im).planes(), im.planes()  );
  EXPECT_EQ( imag(im)(0,0), 2.0 );
  ASSERT_PRESERVED_TYPE( imag );
}

TEST( ImageMath, Abs ) {
  ImageView<std::complex<double> > im(1,1); im(0,0)=std::complex<double>(1.0,2.0);
  EXPECT_EQ( abs(im).cols(), im.cols()  );
  EXPECT_EQ( abs(im).rows(), im.rows()  );
  EXPECT_EQ( abs(im).planes(), im.planes()  );
  EXPECT_NEAR( abs(im)(0,0), 2.23607, 1e-5 );
  ASSERT_PRESERVED_TYPE( abs );
}

TEST( ImageMath, Conj ) {
  ImageView<std::complex<double> > im(1,1); im(0,0)=std::complex<double>(1.0,2.0);
  EXPECT_EQ( conj(im).cols(), im.cols()  );
  EXPECT_EQ( conj(im).rows(), im.rows()  );
  EXPECT_EQ( conj(im).planes(), im.planes()  );
  EXPECT_EQ( conj(im)(0,0).real(), 1.0 );
  EXPECT_EQ( conj(im)(0,0).imag(), -2.0 );


  ASSERT_TRUE( has_pixel_type<float>( conj(ImageView<float>()) ) );
  ASSERT_TRUE( has_pixel_type<double>( conj(ImageView<double>()) ) );
  ASSERT_TRUE( has_pixel_type<long double>( conj(ImageView<long double>()) ) );
  ASSERT_TRUE( has_pixel_type<int>( conj(ImageView<int>()) ) );
  ASSERT_TRUE( has_pixel_type<std::complex<float> >( conj(ImageView<std::complex<float> >()) ) );
  ASSERT_TRUE( has_pixel_type<std::complex<double> >( conj(ImageView<std::complex<double> >()) ) );
  ASSERT_TRUE( has_pixel_type<std::complex< long double> >( conj(ImageView<std::complex<long double> >()) ) );
  ASSERT_TRUE( has_pixel_type<std::complex<int> >( conj(ImageView<std::complex<int> >()) ) );
}

TEST( ImageMath, Square ) {
  ImageView<double> im(2,1);
  im(0,0) = 3; im(1,0) = -4;
  EXPECT_EQ( square(im).cols(), im.cols() );
  EXPECT_EQ( square(im).rows(), im.rows() );
  EXPECT_EQ( square(im).planes(), im.planes() );
  EXPECT_EQ( square(im)(0,0), 9 );
  EXPECT_EQ( square(im)(1,0), 16 );
  ASSERT_TRUE( has_pixel_type<float>( square(ImageView<float>()) ) );
  ASSERT_TRUE( has_pixel_type<double>( square(ImageView<double>()) ) );
  ASSERT_TRUE( has_pixel_type<int32>( square(ImageView<int32>()) ) );
  ASSERT_TRUE( has_pixel_type<PixelGray<float> >( square(ImageView<PixelGray<float> >()) ) );
  ASSERT_TRUE( has_pixel_type<PixelRGB<float> >( square(ImageView<PixelRGB<float> >()) ) );
  ASSERT_TRUE( has_pixel_type<PixelGray<int32> >( square(ImageView<PixelGray<int32> >()) ) );
  ASSERT_TRUE( has_pixel_type<PixelRGB<int32> >( square(ImageView<PixelRGB<int32> >()) ) );
}

TEST( ImageMath, MaskMath ) {
  ImageView<PixelMask<PixelGray<float> > > image(2,1);
  for ( int i = 0; i < 2; i++ )
    image(i,0) = .5;
  image(1,0).invalidate();

  ASSERT_TRUE( is_valid(image(0,0)) == true );
  ASSERT_TRUE( is_valid(image(1,0)) == false );

  for ( int i = 0; i < 2; i++ )
    ASSERT_TRUE( image(i,0) == .5 );

  image += .25;
  for ( int i = 0; i < 2; i++ )
    ASSERT_TRUE( image(i,0) == .75 );

  image /= 2;
  for ( int i = 0; i < 2; i++ )
    ASSERT_TRUE( image(i,0) == .375 );

  image *= 2;
  for ( int i = 0; i < 2; i++ )
    ASSERT_TRUE( image(i,0) == .75 );

  image -= .5;
  for ( int i = 0; i < 2; i++ )
    ASSERT_TRUE( image(i,0) == .25 );

  for ( int i = 0; i < 2; i++ )
    EXPECT_EQ( image(i,0), .25 );
}

#define TEST_UNARY_MATH_FUNCTION(func,arg,result)                       \
  do {                                                                  \
    ImageView<double> im(1,1); im(0,0)=(arg);                           \
    EXPECT_EQ( func(im).cols(), im.cols() );                     \
    EXPECT_EQ( func(im).rows(), im.rows() );                     \
    EXPECT_EQ( func(im).planes(), im.planes() );                 \
    EXPECT_NEAR( func(im)(0,0), result, 1e-5 );                  \
    ASSERT_TRUE( has_pixel_type<float>( func(ImageView<float>()) ) );     \
    ASSERT_TRUE( has_pixel_type<double>( func(ImageView<double>()) ) );   \
    ASSERT_TRUE( has_pixel_type<long double>( func(ImageView<long double>()) ) ); \
    ASSERT_TRUE( has_pixel_type<double>( func(ImageView<int>()) ) );      \
  } while(false)

#define TEST_BINARY_MATH_FUNCTION(func,arg1,arg2,result)                \
  do {                                                                  \
    ImageView<double> im1(1,1); im1(0,0)=(arg1);                        \
    ImageView<double> im2(1,1); im2(0,0)=(arg2);                        \
    EXPECT_EQ( func(im1,im2).cols(), im1.cols() );               \
    EXPECT_EQ( func(im1,im2).rows(), im1.rows() );               \
    EXPECT_EQ( func(im1,im2).planes(), im1.planes() );           \
    EXPECT_NEAR( func(im1,im2)(0,0), result, 1e-5 );             \
    EXPECT_EQ( func(im1,arg2).cols(), im1.cols() );              \
    EXPECT_EQ( func(im1,arg2).rows(), im1.rows() );              \
    EXPECT_EQ( func(im1,arg2).planes(), im1.planes() );          \
    EXPECT_NEAR( func(im1,arg2)(0,0), result, 1e-5 );            \
    EXPECT_EQ( func(arg1,im2).cols(), im1.cols() );              \
    EXPECT_EQ( func(arg1,im2).rows(), im1.rows() );              \
    EXPECT_EQ( func(arg1,im2).planes(), im1.planes() );          \
    EXPECT_NEAR( func(arg1,im2)(0,0), result, 1e-5 );            \
    ASSERT_TRUE( has_pixel_type<float>( func(ImageView<float>(),ImageView<float>()) ) ); \
    ASSERT_TRUE( has_pixel_type<double>( func(ImageView<double>(),ImageView<double>()) ) ); \
    ASSERT_TRUE( has_pixel_type<long double>( func(ImageView<long double>(),ImageView<long double>()) ) ); \
    ASSERT_TRUE( has_pixel_type<double>( func(ImageView<int>(),ImageView<int>()) ) );  \
    ASSERT_TRUE( has_pixel_type<float>( func(ImageView<int>(),ImageView<float>()) ) ); \
    ASSERT_TRUE( has_pixel_type<double>( func(ImageView<float>(),ImageView<double>()) ) ); \
    ASSERT_TRUE( has_pixel_type<long double>( func(ImageView<double>(),ImageView<long double>()) ) ); \
    ASSERT_TRUE( has_pixel_type<double>( func(ImageView<char>(),ImageView<int>()) ) );  \
    ASSERT_TRUE( has_pixel_type<float>( func(ImageView<float>(),float()) ) ); \
    ASSERT_TRUE( has_pixel_type<double>( func(ImageView<double>(),double()) ) ); \
    ASSERT_TRUE( has_pixel_type<double>( func(ImageView<int>(),int()) ) );  \
    ASSERT_TRUE( has_pixel_type<float>( func(ImageView<int>(),float()) ) ); \
    ASSERT_TRUE( has_pixel_type<double>( func(ImageView<float>(),double()) ) ); \
    ASSERT_TRUE( has_pixel_type<long double>( func(ImageView<long double>(),double()) ) ); \
    ASSERT_TRUE( has_pixel_type<double>( func(ImageView<char>(),int()) ) );  \
    ASSERT_TRUE( has_pixel_type<float>( func(float(),ImageView<float>()) ) ); \
    ASSERT_TRUE( has_pixel_type<double>( func(double(),ImageView<double>()) ) ); \
    ASSERT_TRUE( has_pixel_type<double>( func(int(),ImageView<int>()) ) );  \
    ASSERT_TRUE( has_pixel_type<float>( func(int(),ImageView<float>()) ) ); \
    ASSERT_TRUE( has_pixel_type<double>( func(float(),ImageView<double>()) ) ); \
    ASSERT_TRUE( has_pixel_type<long double>( func(double(),ImageView<long double>()) ) ); \
    ASSERT_TRUE( has_pixel_type<double>( func(char(),ImageView<int>()) ) );  \
  } while(false)

TEST( ImageMath, ACOS ) { TEST_UNARY_MATH_FUNCTION(acos,0.5,1.0472); }
TEST( ImageMath, ASIN ) { TEST_UNARY_MATH_FUNCTION(asin,0.5,0.523599); }
TEST( ImageMath, ATAN ) { TEST_UNARY_MATH_FUNCTION(atan,1.0,0.785398); }
TEST( ImageMath, COS ) { TEST_UNARY_MATH_FUNCTION(cos,1.0,0.540302); }
TEST( ImageMath, SIN ) { TEST_UNARY_MATH_FUNCTION(sin,1.0,0.841471); }
TEST( ImageMath, TAN ) { TEST_UNARY_MATH_FUNCTION(tan,1.0,1.55741); }
TEST( ImageMath, COSH ) { TEST_UNARY_MATH_FUNCTION(cosh,1.0,1.54308); }
TEST( ImageMath, SINH ) { TEST_UNARY_MATH_FUNCTION(sinh,1.0,1.1752); }
TEST( ImageMath, TANH ) { TEST_UNARY_MATH_FUNCTION(tanh,1.0,0.761594); }
TEST( ImageMath, EXP ) { TEST_UNARY_MATH_FUNCTION(exp,1.0,2.718281); }
TEST( ImageMath, LOG ) { TEST_UNARY_MATH_FUNCTION(log,2.0,0.693147); }
TEST( ImageMath, LOG10 ) { TEST_UNARY_MATH_FUNCTION(log10,2.0,0.30103); }
TEST( ImageMath, SQRT ) { TEST_UNARY_MATH_FUNCTION(sqrt,2.0,1.41421); }
TEST( ImageMath, CEIL ) { TEST_UNARY_MATH_FUNCTION(ceil,1.5,2.0);
  TEST_UNARY_MATH_FUNCTION(ceil,-1.5,-1.0); }
TEST( ImageMath, FLOOR ) { TEST_UNARY_MATH_FUNCTION(floor,1.5,1.0);
  TEST_UNARY_MATH_FUNCTION(floor,-1.5,-2.0); }
TEST( ImageMath, ATAN2 ) { TEST_BINARY_MATH_FUNCTION(atan2,2.0,1.0,1.10715); }
TEST( ImageMath, POW ) { TEST_BINARY_MATH_FUNCTION(pow,3.0,2.0,9.0); }

#ifndef WIN32
TEST( ImageMath, ACOSH ) { TEST_UNARY_MATH_FUNCTION(acosh,1.5,0.962424); }
TEST( ImageMath, ASINH ){ TEST_UNARY_MATH_FUNCTION(asinh,1.0,0.881374); }
TEST( ImageMath, ATANH ) { TEST_UNARY_MATH_FUNCTION(atanh,0.5,0.549306); }
TEST( ImageMath, EXP2 ) {
#ifdef VW_HAVE_EXP2
  TEST_UNARY_MATH_FUNCTION(exp2,1.0,2.0);
#endif
}
TEST( ImageMath, EXPM1 ) { TEST_UNARY_MATH_FUNCTION(expm1,1.0,1.718281); }
TEST( ImageMath, LOG2 ) {
#ifdef VW_HAVE_LOG2
  TEST_UNARY_MATH_FUNCTION(log2,2.0,1.0);
#endif
}
TEST( ImageMath, LOG1P ) { TEST_UNARY_MATH_FUNCTION(log1p,1.0,0.693147); }
TEST( ImageMath, CBRT ) { TEST_UNARY_MATH_FUNCTION(cbrt,2.0,1.25992); }
TEST( ImageMath, ERF ) { TEST_UNARY_MATH_FUNCTION(erf,1.0,0.842701); }
TEST( ImageMath, ERFC ) { TEST_UNARY_MATH_FUNCTION(erfc,1.0,0.157299); }
TEST( ImageMath, TGAMMA ) {
#ifdef VW_HAVE_TGAMMA
  TEST_UNARY_MATH_FUNCTION(tgamma,1.5,0.886227);
#endif
}
TEST( ImageMath, LGAMMA ) { TEST_UNARY_MATH_FUNCTION(lgamma,2.5,0.284683); }
TEST( ImageMath, ROUND ) { TEST_UNARY_MATH_FUNCTION(round,1.4,1.0);
  TEST_UNARY_MATH_FUNCTION(round,1.5,2.0); }
TEST( ImageMath, TRUNC ) { TEST_UNARY_MATH_FUNCTION(trunc,1.5,1.0);
  TEST_UNARY_MATH_FUNCTION(trunc,-1.5,-1.0); }

TEST( ImageMath, HYPOT ) { TEST_BINARY_MATH_FUNCTION(hypot,2.0,1.0,2.23607); }
TEST( ImageMath, COPYSIGN ) { TEST_BINARY_MATH_FUNCTION(copysign,3.0,-2.0,-3.0);
  TEST_BINARY_MATH_FUNCTION(copysign,3.0,2.0,3.0); }
TEST( ImageMath, FDIM ) { TEST_BINARY_MATH_FUNCTION(fdim,3.0,2.0,1.0);
  TEST_BINARY_MATH_FUNCTION(fdim,2.0,3.0,0.0); }
#endif
