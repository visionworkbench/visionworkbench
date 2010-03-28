// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestImageMath.h
#include <gtest/gtest.h>

#include <vw/config.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/Manipulation.h>

#include <vw/Image/PixelMask.h>
#include <vw/Image/MaskViews.h>

using namespace vw;

template <template<class> class TraitT, class T>
static bool bool_trait( T const& /*arg*/ ) {
  return TraitT<T>::value;
}

template <class T1, class T2>
static bool has_pixel_type( T2 ) {
  return boost::is_same<T1,typename T2::pixel_type>::value;
}

TEST( ImageMath, Negation ) {
  ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;

  ImageView<double> im2 = -im1;
  EXPECT_EQ( im2(0,0), -1 );
  EXPECT_EQ( im2(1,0), -2 );
  EXPECT_EQ( im2(0,1), -3 );
  EXPECT_EQ( im2(1,1), -4 );

  // Test the traits
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( -im1 ) );
}

TEST( ImageMath, Sum ) {
  ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
  ImageView<double> im2(2,2); im2(0,0)=0; im2(1,0)=3; im2(0,1)=2; im2(1,1)=9;

  ImageView<double> im3 = im1 + im2;
  EXPECT_EQ( im3(0,0), 1 );
  EXPECT_EQ( im3(1,0), 5 );
  EXPECT_EQ( im3(0,1), 5 );
  EXPECT_EQ( im3(1,1), 13 );

  im3 = im1 + 2;
  EXPECT_EQ( im3(0,0), 3 );
  EXPECT_EQ( im3(1,0), 4 );
  EXPECT_EQ( im3(0,1), 5 );
  EXPECT_EQ( im3(1,1), 6 );

  im3 = 1 + im1;
  EXPECT_EQ( im3(0,0), 2 );
  EXPECT_EQ( im3(1,0), 3 );
  EXPECT_EQ( im3(0,1), 4 );
  EXPECT_EQ( im3(1,1), 5 );

  im3 += im1;
  EXPECT_EQ( im3(0,0), 3 );
  EXPECT_EQ( im3(1,0), 5 );
  EXPECT_EQ( im3(0,1), 7 );
  EXPECT_EQ( im3(1,1), 9 );

  im3 += 2;
  EXPECT_EQ( im3(0,0), 5 );
  EXPECT_EQ( im3(1,0), 7 );
  EXPECT_EQ( im3(0,1), 9 );
  EXPECT_EQ( im3(1,1), 11 );

  (im3 += 2) += 2;
  EXPECT_EQ( im3(0,0), 9 );
  EXPECT_EQ( im3(1,0), 11 );
  EXPECT_EQ( im3(0,1), 13 );
  EXPECT_EQ( im3(1,1), 15 );

  crop(im3,0,1,2,1) -= 2;
  EXPECT_EQ( im3(0,0), 9 );
  EXPECT_EQ( im3(1,0), 11 );
  EXPECT_EQ( im3(0,1), 11 );
  EXPECT_EQ( im3(1,1), 13 );

  // Test the traits
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( im1 + im2 ) );
}

TEST( ImageMath, Difference ) {
  ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
  ImageView<double> im2(2,2); im2(0,0)=0; im2(1,0)=3; im2(0,1)=2; im2(1,1)=9;

  ImageView<double> im3 = im1 - im2;
  EXPECT_EQ( im3(0,0), 1 );
  EXPECT_EQ( im3(1,0), -1 );
  EXPECT_EQ( im3(0,1), 1 );
  EXPECT_EQ( im3(1,1), -5 );

  im3 = im1 - 2;
  EXPECT_EQ( im3(0,0), -1 );
  EXPECT_EQ( im3(1,0), 0 );
  EXPECT_EQ( im3(0,1), 1 );
  EXPECT_EQ( im3(1,1), 2 );

  im3 = 1 - im1;
  EXPECT_EQ( im3(0,0), 0 );
  EXPECT_EQ( im3(1,0), -1 );
  EXPECT_EQ( im3(0,1), -2 );
  EXPECT_EQ( im3(1,1), -3 );

  im3 -= im1;
  EXPECT_EQ( im3(0,0), -1 );
  EXPECT_EQ( im3(1,0), -3 );
  EXPECT_EQ( im3(0,1), -5 );
  EXPECT_EQ( im3(1,1), -7 );

  im3 -= 2;
  EXPECT_EQ( im3(0,0), -3 );
  EXPECT_EQ( im3(1,0), -5 );
  EXPECT_EQ( im3(0,1), -7 );
  EXPECT_EQ( im3(1,1), -9 );

  (im3 -= 2 ) -= 2;
  EXPECT_EQ( im3(0,0), -7 );
  EXPECT_EQ( im3(1,0), -9 );
  EXPECT_EQ( im3(0,1), -11 );
  EXPECT_EQ( im3(1,1), -13 );

  // Test the traits
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( im1 - im2 ) );
}

TEST( ImageMath, Product ) {
  ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
  ImageView<double> im2(2,2); im2(0,0)=0; im2(1,0)=3; im2(0,1)=2; im2(1,1)=9;

  ImageView<double> im3 = im1 * im2;
  EXPECT_EQ( im3(0,0), 0 );
  EXPECT_EQ( im3(1,0), 6 );
  EXPECT_EQ( im3(0,1), 6 );
  EXPECT_EQ( im3(1,1), 36 );

  im3 = im1 * 2;
  EXPECT_EQ( im3(0,0), 2 );
  EXPECT_EQ( im3(1,0), 4 );
  EXPECT_EQ( im3(0,1), 6 );
  EXPECT_EQ( im3(1,1), 8 );

  im3 = 3 * im1;
  EXPECT_EQ( im3(0,0), 3 );
  EXPECT_EQ( im3(1,0), 6 );
  EXPECT_EQ( im3(0,1), 9 );
  EXPECT_EQ( im3(1,1), 12 );

  im3 *= im1;
  EXPECT_EQ( im3(0,0), 3 );
  EXPECT_EQ( im3(1,0), 12 );
  EXPECT_EQ( im3(0,1), 27 );
  EXPECT_EQ( im3(1,1), 48 );

  im3 *= 0.5;
  EXPECT_EQ( im3(0,0), 1.5 );
  EXPECT_EQ( im3(1,0), 6 );
  EXPECT_EQ( im3(0,1), 13.5 );
  EXPECT_EQ( im3(1,1), 24 );

  (im3 *= 2) *= 2;

  EXPECT_EQ( im3(0,0), 6 );
  EXPECT_EQ( im3(1,0), 24 );
  EXPECT_EQ( im3(0,1), 54 );
  EXPECT_EQ( im3(1,1), 96 );

  // Test the traits
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( im1 * im2 ) );
}

TEST( ImageMath, Quotient ) {
  ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
  ImageView<double> im2(2,2); im2(0,0)=1; im2(1,0)=4; im2(0,1)=2; im2(1,1)=8;

  ImageView<double> im3 = im1 / im2;
  EXPECT_EQ( im3(0,0), 1 );
  EXPECT_EQ( im3(1,0), 0.5 );
  EXPECT_EQ( im3(0,1), 1.5 );
  EXPECT_EQ( im3(1,1), 0.5 );

  im3 = im1 / 2;
  EXPECT_EQ( im3(0,0), 0.5 );
  EXPECT_EQ( im3(1,0), 1 );
  EXPECT_EQ( im3(0,1), 1.5 );
  EXPECT_EQ( im3(1,1), 2 );

  im3 = 2 / im2;
  EXPECT_EQ( im3(0,0), 2 );
  EXPECT_EQ( im3(1,0), 0.5 );
  EXPECT_EQ( im3(0,1), 1 );
  EXPECT_EQ( im3(1,1), 0.25 );

  im3 /= im2;
  EXPECT_EQ( im3(0,0), 2 );
  EXPECT_EQ( im3(1,0), 0.125 );
  EXPECT_EQ( im3(0,1), 0.5 );
  EXPECT_EQ( im3(1,1), 0.03125 );

  im3 /= 0.5;
  EXPECT_EQ( im3(0,0), 4 );
  EXPECT_EQ( im3(1,0), 0.25 );
  EXPECT_EQ( im3(0,1), 1 );
  EXPECT_EQ( im3(1,1), 0.0625 );

  (im3 /= 0.5) /= 0.5;
  EXPECT_EQ( im3(0,0), 16 );
  EXPECT_EQ( im3(1,0), 1 );
  EXPECT_EQ( im3(0,1), 4 );
  EXPECT_EQ( im3(1,1), 0.25 );

  // Test the traits
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( im1 / im2 ) );
}

TEST( ImageMath, Real ) {
  ImageView<std::complex<double> > im(1,1); im(0,0)=std::complex<double>(1.0,2.0);
  EXPECT_EQ( real(im).cols(), im.cols()  );
  EXPECT_EQ( real(im).rows(), im.rows()  );
  EXPECT_EQ( real(im).planes(), im.planes()  );
  EXPECT_EQ( real(im)(0,0), 1.0 );
  ASSERT_TRUE( has_pixel_type<float>( real(ImageView<float>()) ) );
  ASSERT_TRUE( has_pixel_type<double>( real(ImageView<double>()) ) );
  ASSERT_TRUE( has_pixel_type<long double>( real(ImageView<long double>()) ) );
  ASSERT_TRUE( has_pixel_type<int>( real(ImageView<int>()) ) );
  ASSERT_TRUE( has_pixel_type<float>( real(ImageView<std::complex<float> >()) ) );
  ASSERT_TRUE( has_pixel_type<double>( real(ImageView<std::complex<double> >()) ) );
  ASSERT_TRUE( has_pixel_type<long double>( real(ImageView<std::complex<long double> >()) ) );
  ASSERT_TRUE( has_pixel_type<int>( real(ImageView<std::complex<int> >()) ) );
}

TEST( ImageMath, Imag ) {
  ImageView<std::complex<double> > im(1,1);
  im(0,0)=std::complex<double>(1.0,2.0);
  EXPECT_EQ( imag(im).cols(), im.cols()  );
  EXPECT_EQ( imag(im).rows(), im.rows()  );
  EXPECT_EQ( imag(im).planes(), im.planes()  );
  EXPECT_EQ( imag(im)(0,0), 2.0 );
  ASSERT_TRUE( has_pixel_type<float>( imag(ImageView<float>()) ) );
  ASSERT_TRUE( has_pixel_type<double>( imag(ImageView<double>()) ) );
  ASSERT_TRUE( has_pixel_type<long double>( imag(ImageView<long double>()) ) );
  ASSERT_TRUE( has_pixel_type<int>( imag(ImageView<int>()) ) );
  ASSERT_TRUE( has_pixel_type<float>( imag(ImageView<std::complex<float> >()) ) );
  ASSERT_TRUE( has_pixel_type<double>( imag(ImageView<std::complex<double> >()) ) );
  ASSERT_TRUE( has_pixel_type<long double>( imag(ImageView<std::complex<long double> >()) ) );
  ASSERT_TRUE( has_pixel_type<int>( imag(ImageView<std::complex<int> >()) ) );
}

TEST( ImageMath, Abs ) {
  ImageView<std::complex<double> > im(1,1); im(0,0)=std::complex<double>(1.0,2.0);
  EXPECT_EQ( abs(im).cols(), im.cols()  );
  EXPECT_EQ( abs(im).rows(), im.rows()  );
  EXPECT_EQ( abs(im).planes(), im.planes()  );
  EXPECT_NEAR( abs(im)(0,0), 2.23607, 1e-5 );
  ASSERT_TRUE( has_pixel_type<float>( abs(ImageView<float>()) ) );
  ASSERT_TRUE( has_pixel_type<double>( abs(ImageView<double>()) ) );
  ASSERT_TRUE( has_pixel_type<long double>( abs(ImageView<long double>()) ) );
  ASSERT_TRUE( has_pixel_type<int>( abs(ImageView<int>()) ) );
  ASSERT_TRUE( has_pixel_type<float>( abs(ImageView<std::complex<float> >()) ) );
  ASSERT_TRUE( has_pixel_type<double>( abs(ImageView<std::complex<double> >()) ) );
  ASSERT_TRUE( has_pixel_type<long double>( abs(ImageView<std::complex<long double> >()) ) );
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

  ASSERT_TRUE( is_valid(image(0,0)) );
  ASSERT_FALSE( is_valid(image(1,0)) );
  image = invert_mask( image );
  ASSERT_FALSE( is_valid(image(0,0)) );
  ASSERT_TRUE( is_valid(image(1,0)) );
  image = invert_mask( image );
  ASSERT_TRUE( is_valid(image(0,0)) );
  ASSERT_FALSE( is_valid(image(1,0)) );
  image = validate_mask( image );
  ASSERT_TRUE( is_valid(image(0,0)) );
  ASSERT_TRUE( is_valid(image(1,0)) );
  image = invalidate_mask( image );
  ASSERT_FALSE( is_valid(image(0,0)) );
  ASSERT_FALSE( is_valid(image(1,0)) );

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
