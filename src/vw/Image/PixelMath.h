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

/// \file PixelMath.h
/// 
/// Pixel type math base classes and supporting code.
///
#ifndef __VW_IMAGE_PIXELMATH_H__
#define __VW_IMAGE_PIXELMATH_H__

#include <boost/utility/enable_if.hpp>

#include <vw/Core/CompoundTypes.h>
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
    DerivedT& impl() { return static_cast<DerivedT&>(*this); }
    DerivedT const& impl() const { return static_cast<DerivedT const&>(*this); }
    template <class ArgT> inline DerivedT& operator+=( ArgT a ) { return *static_cast<DerivedT*>(this) = (static_cast<DerivedT const&>(*this) + a ); }
    template <class ArgT> inline DerivedT& operator-=( ArgT a ) { return *static_cast<DerivedT*>(this) = (static_cast<DerivedT const&>(*this) - a ); }
    template <class ArgT> inline DerivedT& operator*=( ArgT a ) { return *static_cast<DerivedT*>(this) = (static_cast<DerivedT const&>(*this) * a ); }
    template <class ArgT> inline DerivedT& operator/=( ArgT a ) { return *static_cast<DerivedT*>(this) = (static_cast<DerivedT const&>(*this) / a ); }
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
  typename boost::disable_if< IsCompound<ScalarT>, typename CompoundResult<ftor<ScalarT>,PixelT>::type >::type \
  inline func( PixelMathBase<PixelT> const& pixel, ScalarT scalar ) {   \
    return compound_apply(ftor<ScalarT>(scalar), pixel.impl() );        \
  }

#define VW_PIXEL_MATH_BINARY_SP_FUNCTION(func,ftor)                     \
  template <class PixelT, class ScalarT>                                \
  typename boost::disable_if< IsCompound<ScalarT>, typename CompoundResult<ftor<ScalarT>,PixelT>::type >::type \
  inline func( ScalarT scalar, PixelMathBase<PixelT> const& pixel ) {   \
    return compound_apply(ftor<ScalarT>(scalar), pixel.impl() );        \
  }


  // *******************************************************************
  // Default mathematical operator overlaods
  // *******************************************************************

  VW_PIXEL_MATH_UNARY_FUNCTION(operator -, vw::ArgNegationFunctor)
  VW_PIXEL_MATH_BINARY_PP_FUNCTION(operator +, vw::ArgArgSumFunctor)
  VW_PIXEL_MATH_BINARY_PS_FUNCTION(operator +, vw::ArgValSumFunctor)
  VW_PIXEL_MATH_BINARY_SP_FUNCTION(operator +, vw::ValArgSumFunctor)
  VW_PIXEL_MATH_BINARY_PP_FUNCTION(operator -, vw::ArgArgDifferenceFunctor)
  VW_PIXEL_MATH_BINARY_PS_FUNCTION(operator -, vw::ArgValDifferenceFunctor)
  VW_PIXEL_MATH_BINARY_SP_FUNCTION(operator -, vw::ValArgDifferenceFunctor)
  VW_PIXEL_MATH_BINARY_PP_FUNCTION(operator *, vw::ArgArgProductFunctor)
  VW_PIXEL_MATH_BINARY_PS_FUNCTION(operator *, vw::ArgValProductFunctor)
  VW_PIXEL_MATH_BINARY_SP_FUNCTION(operator *, vw::ValArgProductFunctor)
  VW_PIXEL_MATH_BINARY_PP_FUNCTION(operator /, vw::ArgArgQuotientFunctor)
  VW_PIXEL_MATH_BINARY_PS_FUNCTION(operator /, vw::ArgValQuotientFunctor)
  VW_PIXEL_MATH_BINARY_SP_FUNCTION(operator /, vw::ValArgQuotientFunctor)


  // *******************************************************************
  // Default mathematical function overlaods
  // *******************************************************************

  VW_PIXEL_MATH_UNARY_FUNCTION(acos, vw::math::ArgAcosFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(asin, vw::math::ArgAsinFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(atan, vw::math::ArgAtanFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(cos, vw::math::ArgCosFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(sin, vw::math::ArgSinFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(tan, vw::math::ArgTanFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(cosh, vw::math::ArgCoshFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(sinh, vw::math::ArgSinhFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(tanh, vw::math::ArgTanhFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(exp, vw::math::ArgExpFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(log, vw::math::ArgLogFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(log10, vw::math::ArgLog10Functor)
  VW_PIXEL_MATH_UNARY_FUNCTION(sqrt, vw::math::ArgSqrtFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(ceil, vw::math::ArgCeilFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(floor, vw::math::ArgFloorFunctor)

  VW_PIXEL_MATH_BINARY_PP_FUNCTION(atan2, vw::math::ArgArgAtan2Functor)
  VW_PIXEL_MATH_BINARY_PS_FUNCTION(atan2, vw::math::ArgValAtan2Functor)
  VW_PIXEL_MATH_BINARY_SP_FUNCTION(atan2, vw::math::ValArgAtan2Functor)
  VW_PIXEL_MATH_BINARY_PP_FUNCTION(pow, vw::math::ArgArgPowFunctor)
  VW_PIXEL_MATH_BINARY_PS_FUNCTION(pow, vw::math::ArgValPowFunctor)
  VW_PIXEL_MATH_BINARY_SP_FUNCTION(pow, vw::math::ValArgPowFunctor)

#ifndef WIN32
  VW_PIXEL_MATH_UNARY_FUNCTION(acosh, vw::math::ArgAcoshFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(asinh, vw::math::ArgAsinhFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(atanh, vw::math::ArgAtanhFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(exp2, vw::math::ArgExp2Functor)
  VW_PIXEL_MATH_UNARY_FUNCTION(expm1, vw::math::ArgExpm1Functor)
  VW_PIXEL_MATH_UNARY_FUNCTION(log2, vw::math::ArgLog2Functor)
  VW_PIXEL_MATH_UNARY_FUNCTION(log1p, vw::math::ArgLog1pFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(cbrt, vw::math::ArgCbrtFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(erf, vw::math::ArgErfFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(erfc, vw::math::ArgErfcFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(tgamma,vw::math::ArgTgammaFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(lgamma,vw::math::ArgLgammaFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(round, vw::math::ArgRoundFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(trunc, vw::math::ArgTruncFunctor)

  VW_PIXEL_MATH_BINARY_PP_FUNCTION(hypot, vw::math::ArgArgHypotFunctor)
  VW_PIXEL_MATH_BINARY_PS_FUNCTION(hypot, vw::math::ArgValHypotFunctor)
  VW_PIXEL_MATH_BINARY_SP_FUNCTION(hypot, vw::math::ValArgHypotFunctor)
  VW_PIXEL_MATH_BINARY_PP_FUNCTION(copysign, vw::math::ArgArgCopysignFunctor)
  VW_PIXEL_MATH_BINARY_PS_FUNCTION(copysign, vw::math::ArgValCopysignFunctor)
  VW_PIXEL_MATH_BINARY_SP_FUNCTION(copysign, vw::math::ValArgCopysignFunctor)
  VW_PIXEL_MATH_BINARY_PP_FUNCTION(fdim, vw::math::ArgArgFdimFunctor)
  VW_PIXEL_MATH_BINARY_PS_FUNCTION(fdim, vw::math::ArgValFdimFunctor)
  VW_PIXEL_MATH_BINARY_SP_FUNCTION(fdim, vw::math::ValArgFdimFunctor)
#endif // WIN32

  VW_PIXEL_MATH_UNARY_FUNCTION(real, vw::math::ArgRealFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(imag, vw::math::ArgImagFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(abs,  vw::math::ArgAbsFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(conj, vw::math::ArgConjFunctor)

} // namespace vw

#endif // __VW_IMAGE_PIXELMATH_H__
