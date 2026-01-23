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


/// \file Image/ImageMath.h
///
/// Standard mathematical functions operating on images.
///
/// This header provides overloaded versions of many of the standard
/// C++ mathematical operators and functions to operate on images.
/// The binary functions can operate either on two images (with the
/// same pixel type and dimensions) or on an image and a scalar.  The
/// following operators are currently supported:
///
///  - <TT>-</TT> image
///  - image <TT>+</TT> image
///  - image <TT>+=</TT> image
///  - image <TT>+ </TT>scalar
///  - image <TT>+= </TT>scalar
///  - scalar <TT>+</TT> image
///  - image <TT>-</TT> image
///  - image <TT>-=</TT> image
///  - image <TT>- </TT>scalar
///  - image <TT>-= </TT>scalar
///  - scalar <TT>-</TT> image
///  - image <TT>*</TT> image
///  - image <TT>*=</TT> image
///  - image <TT>* </TT>scalar
///  - image <TT>*= </TT>scalar
///  - scalar <TT>*</TT> image
///  - image <TT>/</TT> image
///  - image <TT>/=</TT> image
///  - image <TT>/ </TT>scalar
///  - image <TT>/= </TT>scalar
///  - scalar <TT>/</TT> image
///
/// Each function returns either an UnaryPerPixelView or a
/// BinaryPerPixelView, as appropriate, so they can be combined
/// efficiently without introducing intermediate image buffers.  The
/// resulting view adopts the same pixel type and dimensions as the
/// source image, except where otherwise noted.  In order for the
/// image functions to work, the corresponding functions for the
/// underlying pixel types must exist.  We currently don't support the
/// comparison operators here.
///
#ifndef __VW_IMAGE_IMAGEMATH_H__
#define __VW_IMAGE_IMAGEMATH_H__

#include <vw/vw_config.h>
#include <vw/Core/Functors.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/PixelMask.h>

namespace vw {

  // *******************************************************************
  // Per-pixel image function definition macros
  // *******************************************************************

#define VW_IMAGE_MATH_UNARY_FUNCTION(func,ftor)       \
  template <class ImageT>                             \
  inline UnaryPerPixelView<ImageT,ftor> func( ImageViewBase<ImageT> const& image ) { \
    return UnaryPerPixelView<ImageT,ftor>( image.impl() );    \
  }

  // Used to get the pixel math operation result .. not necessarily
  // the math functor solution. This is useful for examples like
  // taking the abs() of a PixelGray and still getting back a
  // PixelGray .. instead of having PixelGray cast to a float and then
  // absoluted by the math functor.
#define VW_IMAGE_MATH_UNARY_PIXELMATH_FUNCTION(func, ftor)           \
  namespace detail {                                                 \
    struct Functorized##func {                                       \
      template <class Args> struct result;                           \
      template <class FuncT, class ValT>                             \
      struct result<FuncT(ValT)> {                                   \
        typedef typename CompoundResult<FuncT, ValT>::type type;     \
      };                                                             \
      template <class ValT>                                          \
      struct result<Functorized##func(ValT)> {                       \
        typedef typename boost::mpl::if_c<IsCompound<ValT>::value, typename CompoundResult<ftor,ValT>::type, typename boost::result_of<ftor(ValT)>::type >::type type; \
      };                                                             \
      template <class ValT>                                          \
      typename boost::enable_if< IsCompound< ValT >, typename result<Functorized##func(ValT)>::type>::type \
      inline operator()( ValT val ) const {                          \
        return vw::func( val );                                      \
      }                                                              \
      template <class ValT>                                          \
      typename boost::disable_if< IsCompound< ValT >, typename result<Functorized##func(ValT)>::type>::type \
      inline operator()( ValT val ) const {                          \
        return ftor()( val );                                        \
      }                                                              \
    };                                                               \
  }                                                                  \
  template <class ImageT>                                            \
  inline UnaryPerPixelView<ImageT, detail::Functorized##func> func( ImageViewBase<ImageT> const& image ) { \
    return UnaryPerPixelView<ImageT, detail::Functorized##func>( image.impl() ); \
  }

#define VW_IMAGE_MATH_BINARY_II_FUNCTION(func,ftor)                     \
  template <class Image1T, class Image2T>                               \
  BinaryPerPixelView<Image1T,Image2T,ftor>                              \
  inline func( ImageViewBase<Image1T> const& image1, ImageViewBase<Image2T> const& image2 ) { \
    return BinaryPerPixelView<Image1T,Image2T,ftor>( image1.impl(), image2.impl() );       \
  }

#define VW_IMAGE_MATH_BINARY_II_FUNCTION_ENABLE_IF( cond, func, ftor )  \
  template <class Image1T, class Image2T>                               \
  typename boost::enable_if_c< cond, BinaryPerPixelView<Image1T,Image2T,ftor> >::type \
  inline func( ImageViewBase<Image1T> const& image1, ImageViewBase<Image2T> const& image2 ) { \
    return BinaryPerPixelView<Image1T,Image2T,ftor>( image1.impl(), image2.impl() );       \
  }

#define VW_IMAGE_MATH_BINARY_IS_FUNCTION(func,ftor)                     \
  template <class ImageT, class ScalarT>                                \
  typename boost::disable_if< IsImageView<ScalarT>, UnaryPerPixelView<ImageT,ftor<ScalarT> > >::type \
  inline func( ImageViewBase<ImageT> const& image, ScalarT scalar ) {   \
    return UnaryPerPixelView<ImageT,ftor<ScalarT> >( image.impl(), ftor<ScalarT>(scalar) ); \
  }

#define VW_IMAGE_MATH_BINARY_IS_FUNCTION_ENABLE_IF( cond, func, ftor )  \
  template <class ImageT, class ScalarT>                                \
  typename boost::enable_if_c< cond && !IsImageView<ScalarT>::value,    \
    UnaryPerPixelView<ImageT,ftor<ScalarT> > >::type                    \
  inline func( ImageViewBase<ImageT> const& image, ScalarT scalar ) {   \
    return UnaryPerPixelView<ImageT,ftor<ScalarT> >( image.impl(), ftor<ScalarT>(scalar) ); \
  }

#define VW_IMAGE_MATH_BINARY_SI_FUNCTION(func,ftor)                     \
  template <class ImageT, class ScalarT>                                \
  typename boost::disable_if< IsImageView<ScalarT>, UnaryPerPixelView<ImageT,ftor<ScalarT> > >::type \
  inline func( ScalarT scalar, ImageViewBase<ImageT> const& image ) {   \
    return UnaryPerPixelView<ImageT,ftor<ScalarT> >( image.impl(), ftor<ScalarT>(scalar) ); \
  }

#define VW_IMAGE_MATH_BINARY_SI_FUNCTION_ENABLE_IF( cond, func, ftor )  \
  template <class ImageT, class ScalarT>                                \
  typename boost::enable_if_c< cond && !IsImageView<ScalarT>::value,    \
    UnaryPerPixelView<ImageT,ftor<ScalarT> > >::type                    \
  inline func( ScalarT scalar, ImageViewBase<ImageT> const& image ) {   \
    return UnaryPerPixelView<ImageT,ftor<ScalarT> >( image.impl(), ftor<ScalarT>(scalar) ); \
  }

#define VW_IMAGE_MATH_BINARY_IP_II_FUNCTION(func,ftor)                  \
  template <class Image1T, class Image2T>                               \
  inline Image1T& func( ImageViewBase<Image1T>& image1, ImageViewBase<Image2T> const& image2 ) { \
    for_each_pixel( image1, image2, ftor() );                           \
    return image1.impl();                                               \
  }                                                                     \
  template <class Image1T, class Image2T>                               \
  inline Image1T const& func( ImageViewBase<Image1T> const& image1, ImageViewBase<Image2T> const& image2 ) { \
    for_each_pixel( image1, image2, ftor() );                           \
    return image1.impl();                                               \
  }

#define VW_IMAGE_MATH_BINARY_IP_II_FUNCTION_ENABLE_IF( cond, func, ftor ) \
  template <class Image1T, class Image2T>                               \
  typename boost::enable_if_c< cond, Image1T& >::type                   \
  inline func( ImageViewBase<Image1T>& image1, ImageViewBase<Image2T> const& image2 ) { \
    for_each_pixel( image1, image2, ftor() );                           \
    return image1.impl();                                               \
  }                                                                     \
  template <class Image1T, class Image2T>                               \
  typename boost::enable_if_c< cond, Image1T const& >::type             \
  inline func( ImageViewBase<Image1T> const& image1, ImageViewBase<Image2T> const& image2 ) { \
    for_each_pixel( image1, image2, ftor() );                           \
    return image1.impl();                                               \
  }

#define VW_IMAGE_MATH_BINARY_IP_IS_FUNCTION(func,ftor)                  \
  template <class ImageT, class ScalarT>                                \
  typename boost::disable_if< IsImageView<ScalarT>, ImageT& >::type     \
  inline func( ImageViewBase<ImageT>& image, ScalarT scalar ) {         \
    for_each_pixel( image, ftor<ScalarT>(scalar) );                     \
    return image.impl();                                                \
  }                                                                     \
  template <class ImageT, class ScalarT>                                \
  typename boost::disable_if< IsImageView<ScalarT>, ImageT const& >::type     \
  inline func( ImageViewBase<ImageT> const& image, ScalarT scalar ) {         \
    for_each_pixel( image, ftor<ScalarT>(scalar) );                     \
    return image.impl();                                                \
  }

#define VW_IMAGE_MATH_BINARY_IP_IS_FUNCTION_ENABLE_IF( cond, func, ftor ) \
  template <class ImageT, class ScalarT>                                \
  typename boost::enable_if_c< cond && !IsImageView<ScalarT>::value, ImageT& >::type \
  inline func( ImageViewBase<ImageT>& image, ScalarT scalar ) {         \
    for_each_pixel( image, ftor<ScalarT>(scalar) );                     \
    return image.impl();                                                \
  }                                                                     \
  template <class ImageT, class ScalarT>                                \
  typename boost::enable_if_c< cond && !IsImageView<ScalarT>::value, ImageT const& >::type \
  inline func( ImageViewBase<ImageT> const& image, ScalarT scalar ) {   \
    for_each_pixel( image, ftor<ScalarT>(scalar) );                     \
    return image.impl();                                                \
  }

  // *******************************************************************
  // Default mathematical operator overlaods
  // *******************************************************************

  /// Negation of an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(operator-, vw::ArgNegationFunctor)

  /// Sum of two images.
  VW_IMAGE_MATH_BINARY_II_FUNCTION(operator +, vw::ArgArgSumFunctor)

  /// Sum of an image and a scalar.
  VW_IMAGE_MATH_BINARY_IS_FUNCTION(operator +, vw::ArgValSumFunctor)

  /// Sum of a scalar and an image.
  VW_IMAGE_MATH_BINARY_SI_FUNCTION(operator +, vw::ValArgSumFunctor)

  /// In-place sum of two images.
  VW_IMAGE_MATH_BINARY_IP_II_FUNCTION(operator +=, vw::ArgArgInPlaceSumFunctor)

  /// In-place sum of an image and a scalar.
  VW_IMAGE_MATH_BINARY_IP_IS_FUNCTION(operator +=, vw::ArgValInPlaceSumFunctor)

  /// Difference of two images.
  VW_IMAGE_MATH_BINARY_II_FUNCTION(operator -, vw::ArgArgDifferenceFunctor)

  /// Difference of an image and a scalar.
  VW_IMAGE_MATH_BINARY_IS_FUNCTION(operator -, vw::ArgValDifferenceFunctor)

  /// Difference of a scalar and an image.
  VW_IMAGE_MATH_BINARY_SI_FUNCTION(operator -, vw::ValArgDifferenceFunctor)

  /// In-place difference of two images.
  VW_IMAGE_MATH_BINARY_IP_II_FUNCTION(operator -=, vw::ArgArgInPlaceDifferenceFunctor)

  /// In-place difference of an image and a scalar.
  VW_IMAGE_MATH_BINARY_IP_IS_FUNCTION(operator -=, vw::ArgValInPlaceDifferenceFunctor)

  /// Product of two images.
  VW_IMAGE_MATH_BINARY_II_FUNCTION(operator *, vw::ArgArgProductFunctor)

  /// Product of an image and a scalar.
  VW_IMAGE_MATH_BINARY_IS_FUNCTION(operator *, vw::ArgValProductFunctor)

  /// Product of a scalar and an image.
  VW_IMAGE_MATH_BINARY_SI_FUNCTION(operator *, vw::ValArgProductFunctor)

  /// In-place product of two images.
  VW_IMAGE_MATH_BINARY_IP_II_FUNCTION(operator *=, vw::ArgArgInPlaceProductFunctor)

  /// In-place product of an image and a scalar.
  VW_IMAGE_MATH_BINARY_IP_IS_FUNCTION(operator *=, vw::ArgValInPlaceProductFunctor)

  /// Quotient of two images.
  VW_IMAGE_MATH_BINARY_II_FUNCTION_ENABLE_IF( !IsMasked<typename Image2T::pixel_type>::value, operator /, vw::ArgArgSafeQuotientFunctor )

  VW_IMAGE_MATH_BINARY_II_FUNCTION_ENABLE_IF( IsMasked<typename Image2T::pixel_type>::value, operator /, vw::ArgArgMaskedSafeQuotientFunctor )

  /// Quotient of an image and a scalar.
  VW_IMAGE_MATH_BINARY_IS_FUNCTION(operator /, vw::ArgValSafeQuotientFunctor)

  /// Quotient of a scalar and an image.
  VW_IMAGE_MATH_BINARY_SI_FUNCTION_ENABLE_IF( !IsMasked<typename ImageT::pixel_type>::value, operator /, vw::ValArgSafeQuotientFunctor )
  VW_IMAGE_MATH_BINARY_SI_FUNCTION_ENABLE_IF( IsMasked<typename ImageT::pixel_type>::value, operator /, vw::ValArgMaskedSafeQuotientFunctor )

  /// In-place quotient of two images.
  VW_IMAGE_MATH_BINARY_IP_II_FUNCTION_ENABLE_IF( !IsMasked<typename Image2T::pixel_type>::value, operator /=, vw::ArgArgInPlaceSafeQuotientFunctor )
  VW_IMAGE_MATH_BINARY_IP_II_FUNCTION_ENABLE_IF( IsMasked<typename Image2T::pixel_type>::value, operator /=, vw::ArgArgInPlaceMaskedSafeQuotientFunctor )

  /// In-place quotient of an image and a scalar.
  VW_IMAGE_MATH_BINARY_IP_IS_FUNCTION(operator /=, vw::ArgValInPlaceSafeQuotientFunctor)


  // *******************************************************************
  // Default mathematical function overlaods
  // *******************************************************************

  /// Computes the arccosine, \f$\cos^{-1} x\f$, of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(acos, vw::math::ArgAcosFunctor)

  /// Computes the arcsine, \f$\sin^{-1} x\f$, of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(asin, vw::math::ArgAsinFunctor)

  /// Computes the arctangent, \f$\tan^{-1} x\f$, of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(atan, vw::math::ArgAtanFunctor)

  /// Computes the cosine, \f$\cos x\f$, of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(cos, vw::math::ArgCosFunctor)

  /// Computes the sine, \f$\sin x\f$, of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(sin, vw::math::ArgSinFunctor)

  /// Computes the tangent, \f$\tan x\f$, of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(tan, vw::math::ArgTanFunctor)

  /// Computes the hyperbolic cosine, \f$\cosh x\f$, of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(cosh, vw::math::ArgCoshFunctor)

  /// Computes the hyperbolic sine, \f$\sinh x\f$, of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(sinh, vw::math::ArgSinhFunctor)

  /// Computes the hyperbolic tangent, \f$\tanh x\f$, of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(tanh, vw::math::ArgTanhFunctor)

  /// Computes the base-e exponential, \f$e^x\f$, of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(exp, vw::math::ArgExpFunctor)

  /// Computes the natural logarithm, \f$\ln x\f$, of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(log, vw::math::ArgLogFunctor)

  /// Computes the base-10 logarithm, \f$\log_{10} x\f$, of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(log10, vw::math::ArgLog10Functor)

  /// Computes the square root of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(sqrt, vw::math::ArgSqrtFunctor)

  /// Computes the smallest integer not less than each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(ceil, vw::math::ArgCeilFunctor)

  /// Computes the largest integer not greater than each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(floor, vw::math::ArgFloorFunctor)

  /// Computes the two-argument arctangent of each pixel from two images.
  VW_IMAGE_MATH_BINARY_II_FUNCTION(atan2, vw::math::ArgArgAtan2Functor)

  /// Computes the two-argument arctangent of each pixel of an image with a scalar.
  VW_IMAGE_MATH_BINARY_IS_FUNCTION(atan2, vw::math::ArgValAtan2Functor)

  /// Computes the two-argument arctangent of a scalar with each pixel of an image.
  VW_IMAGE_MATH_BINARY_SI_FUNCTION(atan2, vw::math::ValArgAtan2Functor)

  /// Computes the power function, \f$x^y\f$, of each pixel from two images.
  VW_IMAGE_MATH_BINARY_II_FUNCTION(pow, vw::math::ArgArgPowFunctor)

  /// Computes the given power of each pixel of an image.
  VW_IMAGE_MATH_BINARY_IS_FUNCTION(pow, vw::math::ArgValPowFunctor)

  /// Raises a scalar to the power of each pixel of an image.
  VW_IMAGE_MATH_BINARY_SI_FUNCTION(pow, vw::math::ValArgPowFunctor)

  /// Computes the hypotenuse, \f$\sqrt{x^2+y^2}\f$,  of each pixel from two images.
  VW_IMAGE_MATH_BINARY_II_FUNCTION(hypot, vw::math::ArgArgHypotFunctor)

  /// Computes the hypotenuse, \f$\sqrt{x^2+y^2}\f$, of each pixel of an image with a scalar.
  VW_IMAGE_MATH_BINARY_IS_FUNCTION(hypot, vw::math::ArgValHypotFunctor)

  /// Computes the hypotenuse, \f$\sqrt{x^2+y^2}\f$, of a scalar with each pixel of an image.
  VW_IMAGE_MATH_BINARY_SI_FUNCTION(hypot, vw::math::ValArgHypotFunctor)

#ifndef WIN32

  /// Computes the hyperbolic arccosine, \f$\cosh^{-1} x\f$, of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(acosh, vw::math::ArgAcoshFunctor)

  /// Computes the hyperbolic arcsine, \f$\sinh^{-1} x\f$, of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(asinh, vw::math::ArgAsinhFunctor)

  /// Computes the hyperbolic arctangent, \f$\tanh^{-1} x\f$, of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(atanh, vw::math::ArgAtanhFunctor)

#ifdef VW_HAVE_EXP2
  /// Computes the base-2 exponential, \f$2^x\f$, of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(exp2, vw::math::ArgExp2Functor)
#endif

#ifdef VW_HAVE_LOG2
  /// Computes the base-2 logarithm, \f$\log_{2} x\f$, of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(log2, vw::math::ArgLog2Functor)
#endif

#ifdef VW_HAVE_TGAMMA
  /// Computes the gamma function of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(tgamma, vw::math::ArgTgammaFunctor)
#endif

  /// Computes the base-e exponential minus one, \f$e^x-1\f$, of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(expm1, vw::math::ArgExpm1Functor)

  /// Computes the natural logarithm of one plus each pixel in an image, \f$\ln(1+x)\f$.
  VW_IMAGE_MATH_UNARY_FUNCTION(log1p, vw::math::ArgLog1pFunctor)

  /// Computes the cube root of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(cbrt, vw::math::ArgCbrtFunctor)

  /// Computes the error function of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(erf, vw::math::ArgErfFunctor)

  /// Computes the complementary error function of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(erfc, vw::math::ArgErfcFunctor)

  /// Computes the natural lograithm of the absolute value of the gamma function of each pixel in an image.
  VW_IMAGE_MATH_UNARY_FUNCTION(lgamma, vw::math::ArgLgammaFunctor)

  /// Computes the nearest integer to each pixel in an image, rounding half-way cases away from zero.
  VW_IMAGE_MATH_UNARY_FUNCTION(round, vw::math::ArgRoundFunctor)

  /// Takes the integer part of each pixel in an image by simple truncation.
  VW_IMAGE_MATH_UNARY_FUNCTION(trunc, vw::math::ArgTruncFunctor)

  /// Returns an image equal to first argument but with the sign of each pixel adjusted to match the corresponding sign in the second argument.
  VW_IMAGE_MATH_BINARY_II_FUNCTION(copysign, vw::math::ArgArgCopysignFunctor)

  /// Returns an image equal to first argument but with the sign of each pixel adjusted to match the sign of the second argument.
  VW_IMAGE_MATH_BINARY_IS_FUNCTION(copysign, vw::math::ArgValCopysignFunctor)

  /// Returns an image equal to second argument but with the magnitude of each pixel adjusted to match the magnitude of the first argument.
  VW_IMAGE_MATH_BINARY_SI_FUNCTION(copysign, vw::math::ValArgCopysignFunctor)

  /// Computes the positive difference of two images, <TT>(x>y)?(x-y):0</TT>.
  VW_IMAGE_MATH_BINARY_II_FUNCTION(fdim, vw::math::ArgArgFdimFunctor)

  /// Computes the positive difference between an image and a scalar, <TT>(x>y)?(x-y):0</TT>.
  VW_IMAGE_MATH_BINARY_IS_FUNCTION(fdim, vw::math::ArgValFdimFunctor)

  /// Computes the positive difference between a scalar and an image, <TT>(x>y)?(x-y):0</TT>.
  VW_IMAGE_MATH_BINARY_SI_FUNCTION(fdim, vw::math::ValArgFdimFunctor)

#endif // WIN32

  /// Takes the real part of each pixel in an image.  If the source
  /// image is complex then the resulting image has the corresponding
  /// real pixel type.  If the source image is real then this function
  /// effectively performs no operation.
  VW_IMAGE_MATH_UNARY_PIXELMATH_FUNCTION(real, vw::math::ArgRealFunctor)

  /// Takes the complex part of each pixel in an image.  If the source
  /// image is complex then the resulting image has the corresponding
  /// real pixel type.  If the source image is real then this function
  /// effectively returns an image of zeros.
  VW_IMAGE_MATH_UNARY_PIXELMATH_FUNCTION(imag, vw::math::ArgImagFunctor)

  /// Computes the absolute value of each pixel in an image.
  /// If the source image is complex then the resulting image has the
  /// corresponding real pixel type.
  VW_IMAGE_MATH_UNARY_PIXELMATH_FUNCTION(abs, vw::math::ArgAbsFunctor)

  /// Computes the complex conjugate of each pixel in an image.
  /// If the source image is real then this function effectively performs
  /// no operation.
  VW_IMAGE_MATH_UNARY_FUNCTION(conj, vw::math::ArgConjFunctor)

  /// Computes the square of each pixel in an image
  VW_IMAGE_MATH_UNARY_FUNCTION(square, vw::math::ArgSquareFunctor)

} // namespace vw

#endif // __VW_IMAGE_IMAGEMATH_H__
