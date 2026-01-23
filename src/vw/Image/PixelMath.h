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


/// \file PixelMath.h
///
/// Pixel type math base classes and supporting code.
///
#ifndef __VW_IMAGE_PIXELMATH_H__
#define __VW_IMAGE_PIXELMATH_H__

#include <vw/vw_config.h>
#include <boost/utility/enable_if.hpp>

#include <vw/Core/Functors.h>
#include <vw/Math/Functors.h>


namespace vw {

  // *******************************************************************
  // Pixel math base class
  // *******************************************************************

  /// This CRTP base class, and its related free operator overloads,
  /// provide some default pixel math behavior used by pixel types
  /// such as PixelRGB.  To use this class, simply derive your pixel
  /// type from it in the usual CRTP manner.  Your pixel type must
  /// provide operator[] and you must specialize the usual compound
  /// type traits: CompoundChannelType, CompoundNumChannels, and
  /// CompoundChannelCast.
  template <class DerivedT>
  class PixelMathBase {
  public:
    DerivedT      & impl()       { return static_cast<DerivedT      &>(*this); }
    DerivedT const& impl() const { return static_cast<DerivedT const&>(*this); }
  };


  // *******************************************************************
  // Per-channel pixel function definition macros
  // *******************************************************************

#define VW_PIXEL_MATH_UNARY_FUNCTION(func,ftor)       \
  template <class PixelT>                             \
  typename CompoundResult<ftor,PixelT>::type          \
  inline func( PixelMathBase<PixelT> const& pixel ) { \
    return compound_apply( ftor(), pixel.impl() );    \
  }

#define VW_PIXEL_MATH_BINARY_PP_FUNCTION(func,ftor)                     \
  template <class Pixel1T, class Pixel2T>                               \
  typename boost::enable_if< CompoundIsCompatible<Pixel1T,Pixel2T>, typename CompoundResult<ftor,Pixel1T,Pixel2T>::type >::type \
  inline func( PixelMathBase<Pixel1T> const& pixel1, PixelMathBase<Pixel2T> const& pixel2 ) { \
    return compound_apply(ftor(), pixel1.impl(), pixel2.impl() );       \
  }

#define VW_PIXEL_MATH_BINARY_PS_FUNCTION(func,ftor)                     \
  template <class PixelT, class ScalarT>                                \
  typename boost::enable_if< IsScalar<ScalarT>, typename CompoundResult<ftor<ScalarT>,PixelT>::type >::type \
  inline func( PixelMathBase<PixelT> const& pixel, ScalarT scalar ) {   \
    return compound_apply(ftor<ScalarT>(scalar), pixel.impl() );        \
  }

#define VW_PIXEL_MATH_BINARY_IP_FUNCTION(func,ftor)                     \
  template <class Pixel1T, class Pixel2T>                               \
  typename boost::enable_if< CompoundIsCompatible<Pixel1T,Pixel2T>, Pixel1T&>::type \
  inline func( PixelMathBase<Pixel1T>& pixel1, PixelMathBase<Pixel2T> const& pixel2 ) { \
    return compound_apply_in_place(ftor(), pixel1.impl(), pixel2.impl() ); \
  }

#define VW_PIXEL_MATH_BINARY_IS_FUNCTION(func,ftor)                     \
  template <class PixelT, class ScalarT>                                \
  typename boost::enable_if< IsScalar<ScalarT>, PixelT&>::type          \
  inline func( PixelMathBase<PixelT>& pixel, ScalarT scalar ) {         \
    return compound_apply_in_place(ftor<ScalarT>(scalar), pixel.impl() ); \
  }

#define VW_PIXEL_MATH_BINARY_SP_FUNCTION(func,ftor)                     \
  template <class PixelT, class ScalarT>                                \
  typename boost::enable_if< IsScalar<ScalarT>, typename CompoundResult<ftor<ScalarT>,PixelT>::type >::type \
  inline func( ScalarT scalar, PixelMathBase<PixelT> const& pixel ) {   \
    return compound_apply(ftor<ScalarT>(scalar), pixel.impl() );        \
  }

#define VW_PIXEL_MATH_STD_UNARY_FUNCTION(func,ftor)     \
  VW_PIXEL_MATH_UNARY_FUNCTION(func,ftor)               \
  using ::func;

#define VW_PIXEL_MATH_STD_BINARY_PP_FUNCTION(func,ftor) \
  VW_PIXEL_MATH_BINARY_PP_FUNCTION(func,ftor)           \
  using ::func;


  // *******************************************************************
  // Default mathematical operator overlaods
  // *******************************************************************

  VW_PIXEL_MATH_UNARY_FUNCTION(operator -, vw::ArgNegationFunctor)
  VW_PIXEL_MATH_BINARY_PP_FUNCTION(operator +, vw::ArgArgSumFunctor)
  VW_PIXEL_MATH_BINARY_PS_FUNCTION(operator +, vw::ArgValSumFunctor)
  VW_PIXEL_MATH_BINARY_SP_FUNCTION(operator +, vw::ValArgSumFunctor)
  VW_PIXEL_MATH_BINARY_IP_FUNCTION(operator +=, vw::ArgArgInPlaceSumFunctor)
  VW_PIXEL_MATH_BINARY_IS_FUNCTION(operator +=, vw::ArgValInPlaceSumFunctor)
  VW_PIXEL_MATH_BINARY_PP_FUNCTION(operator -, vw::ArgArgDifferenceFunctor)
  VW_PIXEL_MATH_BINARY_PS_FUNCTION(operator -, vw::ArgValDifferenceFunctor)
  VW_PIXEL_MATH_BINARY_SP_FUNCTION(operator -, vw::ValArgDifferenceFunctor)
  VW_PIXEL_MATH_BINARY_IP_FUNCTION(operator -=, vw::ArgArgInPlaceDifferenceFunctor)
  VW_PIXEL_MATH_BINARY_IS_FUNCTION(operator -=, vw::ArgValInPlaceDifferenceFunctor)
  VW_PIXEL_MATH_BINARY_PP_FUNCTION(operator *, vw::ArgArgProductFunctor)
  VW_PIXEL_MATH_BINARY_PS_FUNCTION(operator *, vw::ArgValProductFunctor)
  VW_PIXEL_MATH_BINARY_SP_FUNCTION(operator *, vw::ValArgProductFunctor)
  VW_PIXEL_MATH_BINARY_IP_FUNCTION(operator *=, vw::ArgArgInPlaceProductFunctor)
  VW_PIXEL_MATH_BINARY_IS_FUNCTION(operator *=, vw::ArgValInPlaceProductFunctor)
  VW_PIXEL_MATH_BINARY_PP_FUNCTION(operator /, vw::ArgArgQuotientFunctor)
  VW_PIXEL_MATH_BINARY_PS_FUNCTION(operator /, vw::ArgValQuotientFunctor)
  VW_PIXEL_MATH_BINARY_SP_FUNCTION(operator /, vw::ValArgQuotientFunctor)
  VW_PIXEL_MATH_BINARY_IP_FUNCTION(operator /=, vw::ArgArgInPlaceQuotientFunctor)
  VW_PIXEL_MATH_BINARY_IS_FUNCTION(operator /=, vw::ArgValInPlaceQuotientFunctor)

  template <class Pixel1T, class Pixel2T>
  typename boost::enable_if< CompoundIsCompatible<Pixel1T,Pixel2T>, bool >::type
  inline operator==( PixelMathBase<Pixel1T> const& pixel1, PixelMathBase<Pixel2T> const& pixel2 ) {
    for( size_t c=0; c<CompoundNumChannels<Pixel1T>::value; ++c )
      if( compound_select_channel<typename CompoundChannelType<Pixel1T>::type>(pixel1.impl(),c) !=
          compound_select_channel<typename CompoundChannelType<Pixel2T>::type>(pixel2.impl(),c) )
        return false;
    return true;
  }

  template <class Pixel1T, class Pixel2T>
  typename boost::enable_if< CompoundIsCompatible<Pixel1T,Pixel2T>, bool >::type
  inline operator!=( PixelMathBase<Pixel1T> const& pixel1, PixelMathBase<Pixel2T> const& pixel2 ) {
    return !(pixel1==pixel2);
  }

  // *******************************************************************
  // Default mathematical function overlaods
  // *******************************************************************

  VW_PIXEL_MATH_STD_UNARY_FUNCTION(acos, vw::math::ArgAcosFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(asin, vw::math::ArgAsinFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(atan, vw::math::ArgAtanFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(cos, vw::math::ArgCosFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(sin, vw::math::ArgSinFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(tan, vw::math::ArgTanFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(cosh, vw::math::ArgCoshFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(sinh, vw::math::ArgSinhFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(tanh, vw::math::ArgTanhFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(exp, vw::math::ArgExpFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(log, vw::math::ArgLogFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(log10, vw::math::ArgLog10Functor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(sqrt, vw::math::ArgSqrtFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(ceil, vw::math::ArgCeilFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(floor, vw::math::ArgFloorFunctor)

  VW_PIXEL_MATH_STD_BINARY_PP_FUNCTION(atan2, vw::math::ArgArgAtan2Functor)
  VW_PIXEL_MATH_BINARY_PS_FUNCTION(atan2, vw::math::ArgValAtan2Functor)
  VW_PIXEL_MATH_BINARY_SP_FUNCTION(atan2, vw::math::ValArgAtan2Functor)
  VW_PIXEL_MATH_STD_BINARY_PP_FUNCTION(pow, vw::math::ArgArgPowFunctor)
  VW_PIXEL_MATH_BINARY_PS_FUNCTION(pow, vw::math::ArgValPowFunctor)
  VW_PIXEL_MATH_BINARY_SP_FUNCTION(pow, vw::math::ValArgPowFunctor)

  VW_PIXEL_MATH_STD_BINARY_PP_FUNCTION(hypot, vw::math::ArgArgHypotFunctor)
  VW_PIXEL_MATH_BINARY_PS_FUNCTION(hypot, vw::math::ArgValHypotFunctor)
  VW_PIXEL_MATH_BINARY_SP_FUNCTION(hypot, vw::math::ValArgHypotFunctor)

#ifndef WIN32
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(acosh, vw::math::ArgAcoshFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(asinh, vw::math::ArgAsinhFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(atanh, vw::math::ArgAtanhFunctor)
#ifdef VW_HAVE_EXP2
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(exp2, vw::math::ArgExp2Functor)
#endif
#ifdef VW_HAVE_LOG2
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(log2, vw::math::ArgLog2Functor)
#endif
#ifdef VW_HAVE_TGAMMA
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(tgamma,vw::math::ArgTgammaFunctor)
#endif
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(expm1, vw::math::ArgExpm1Functor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(log1p, vw::math::ArgLog1pFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(cbrt, vw::math::ArgCbrtFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(erf, vw::math::ArgErfFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(erfc, vw::math::ArgErfcFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(lgamma,vw::math::ArgLgammaFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(round, vw::math::ArgRoundFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(trunc, vw::math::ArgTruncFunctor)

  VW_PIXEL_MATH_STD_BINARY_PP_FUNCTION(copysign, vw::math::ArgArgCopysignFunctor)
  VW_PIXEL_MATH_BINARY_PS_FUNCTION(copysign, vw::math::ArgValCopysignFunctor)
  VW_PIXEL_MATH_BINARY_SP_FUNCTION(copysign, vw::math::ValArgCopysignFunctor)
  VW_PIXEL_MATH_STD_BINARY_PP_FUNCTION(fdim, vw::math::ArgArgFdimFunctor)
  VW_PIXEL_MATH_BINARY_PS_FUNCTION(fdim, vw::math::ArgValFdimFunctor)
  VW_PIXEL_MATH_BINARY_SP_FUNCTION(fdim, vw::math::ValArgFdimFunctor)
#endif // WIN32

  VW_PIXEL_MATH_UNARY_FUNCTION(real, vw::math::ArgRealFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(imag, vw::math::ArgImagFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(abs,  vw::math::ArgAbsFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(conj, vw::math::ArgConjFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(square, vw::math::ArgSquareFunctor)

} // namespace vw

#endif // __VW_IMAGE_PIXELMATH_H__
