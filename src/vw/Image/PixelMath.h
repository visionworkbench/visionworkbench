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
  inline func(PixelMathBase<PixelT> const& pixel) { \
    return compound_apply(ftor(), pixel.impl());    \
  }

#define VW_PIXEL_MATH_BINARY_PP_FUNCTION(func,ftor)                     \
  template <class Pixel1T, class Pixel2T>                               \
  typename boost::enable_if< CompoundIsCompatible<Pixel1T,Pixel2T>, typename CompoundResult<ftor,Pixel1T,Pixel2T>::type >::type \
  inline func(PixelMathBase<Pixel1T> const& pixel1, PixelMathBase<Pixel2T> const& pixel2) { \
    return compound_apply(ftor(), pixel1.impl(), pixel2.impl());       \
  }

#define VW_PIXEL_MATH_BINARY_PS_FUNCTION(func,ftor)                     \
  template <class PixelT, class ScalarT>                                \
  typename boost::enable_if< IsScalar<ScalarT>, typename CompoundResult<ftor<ScalarT>,PixelT>::type >::type \
  inline func(PixelMathBase<PixelT> const& pixel, ScalarT scalar) {   \
    return compound_apply(ftor<ScalarT>(scalar), pixel.impl());        \
  }

#define VW_PIXEL_MATH_BINARY_IP_FUNCTION(func,ftor)                     \
  template <class Pixel1T, class Pixel2T>                               \
  typename boost::enable_if< CompoundIsCompatible<Pixel1T,Pixel2T>, Pixel1T&>::type \
  inline func(PixelMathBase<Pixel1T>& pixel1, PixelMathBase<Pixel2T> const& pixel2) { \
    return compound_apply_in_place(ftor(), pixel1.impl(), pixel2.impl()); \
  }

#define VW_PIXEL_MATH_BINARY_IS_FUNCTION(func,ftor)                     \
  template <class PixelT, class ScalarT>                                \
  typename boost::enable_if< IsScalar<ScalarT>, PixelT&>::type          \
  inline func(PixelMathBase<PixelT>& pixel, ScalarT scalar) {         \
    return compound_apply_in_place(ftor<ScalarT>(scalar), pixel.impl()); \
  }

#define VW_PIXEL_MATH_BINARY_SP_FUNCTION(func,ftor)                     \
  template <class PixelT, class ScalarT>                                \
  typename boost::enable_if< IsScalar<ScalarT>, typename CompoundResult<ftor<ScalarT>,PixelT>::type >::type \
  inline func(ScalarT scalar, PixelMathBase<PixelT> const& pixel) {   \
    return compound_apply(ftor<ScalarT>(scalar), pixel.impl());        \
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
  inline operator==(PixelMathBase<Pixel1T> const& pixel1, PixelMathBase<Pixel2T> const& pixel2) {
    for(size_t c=0; c<CompoundNumChannels<Pixel1T>::value; ++c)
      if(compound_select_channel<typename CompoundChannelType<Pixel1T>::type>(pixel1.impl(),c) !=
         compound_select_channel<typename CompoundChannelType<Pixel2T>::type>(pixel2.impl(),c))
        return false;
    return true;
  }

  template <class Pixel1T, class Pixel2T>
  typename boost::enable_if< CompoundIsCompatible<Pixel1T,Pixel2T>, bool >::type
  inline operator!=(PixelMathBase<Pixel1T> const& pixel1, PixelMathBase<Pixel2T> const& pixel2) {
    return !(pixel1==pixel2);
  }

  // *******************************************************************
  // Default mathematical function overlaods
  // *******************************************************************

  VW_PIXEL_MATH_STD_UNARY_FUNCTION(acos, vw::ArgAcosFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(asin, vw::ArgAsinFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(atan, vw::ArgAtanFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(cos, vw::ArgCosFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(sin, vw::ArgSinFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(tan, vw::ArgTanFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(cosh, vw::ArgCoshFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(sinh, vw::ArgSinhFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(tanh, vw::ArgTanhFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(exp, vw::ArgExpFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(log, vw::ArgLogFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(log10, vw::ArgLog10Functor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(sqrt, vw::ArgSqrtFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(ceil, vw::ArgCeilFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(floor, vw::ArgFloorFunctor)

  VW_PIXEL_MATH_STD_BINARY_PP_FUNCTION(atan2, vw::ArgArgAtan2Functor)
  VW_PIXEL_MATH_BINARY_PS_FUNCTION(atan2, vw::ArgValAtan2Functor)
  VW_PIXEL_MATH_BINARY_SP_FUNCTION(atan2, vw::ValArgAtan2Functor)
  VW_PIXEL_MATH_STD_BINARY_PP_FUNCTION(pow, vw::ArgArgPowFunctor)
  VW_PIXEL_MATH_BINARY_PS_FUNCTION(pow, vw::ArgValPowFunctor)
  VW_PIXEL_MATH_BINARY_SP_FUNCTION(pow, vw::ValArgPowFunctor)

  VW_PIXEL_MATH_STD_BINARY_PP_FUNCTION(hypot, vw::ArgArgHypotFunctor)
  VW_PIXEL_MATH_BINARY_PS_FUNCTION(hypot, vw::ArgValHypotFunctor)
  VW_PIXEL_MATH_BINARY_SP_FUNCTION(hypot, vw::ValArgHypotFunctor)

#ifndef WIN32
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(acosh, vw::ArgAcoshFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(asinh, vw::ArgAsinhFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(atanh, vw::ArgAtanhFunctor)
#ifdef VW_HAVE_EXP2
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(exp2, vw::ArgExp2Functor)
#endif
#ifdef VW_HAVE_LOG2
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(log2, vw::ArgLog2Functor)
#endif
#ifdef VW_HAVE_TGAMMA
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(tgamma,vw::ArgTgammaFunctor)
#endif
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(expm1, vw::ArgExpm1Functor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(log1p, vw::ArgLog1pFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(cbrt, vw::ArgCbrtFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(erf, vw::ArgErfFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(erfc, vw::ArgErfcFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(lgamma,vw::ArgLgammaFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(round, vw::ArgRoundFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(trunc, vw::ArgTruncFunctor)

  VW_PIXEL_MATH_STD_BINARY_PP_FUNCTION(copysign, vw::ArgArgCopysignFunctor)
  VW_PIXEL_MATH_BINARY_PS_FUNCTION(copysign, vw::ArgValCopysignFunctor)
  VW_PIXEL_MATH_BINARY_SP_FUNCTION(copysign, vw::ValArgCopysignFunctor)
  VW_PIXEL_MATH_STD_BINARY_PP_FUNCTION(fdim, vw::ArgArgFdimFunctor)
  VW_PIXEL_MATH_BINARY_PS_FUNCTION(fdim, vw::ArgValFdimFunctor)
  VW_PIXEL_MATH_BINARY_SP_FUNCTION(fdim, vw::ValArgFdimFunctor)
#endif // WIN32

  VW_PIXEL_MATH_UNARY_FUNCTION(real, vw::ArgRealFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(imag, vw::ArgImagFunctor)
  VW_PIXEL_MATH_STD_UNARY_FUNCTION(abs,  vw::ArgAbsFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(conj, vw::ArgConjFunctor)
  VW_PIXEL_MATH_UNARY_FUNCTION(square, vw::ArgSquareFunctor)

} // namespace vw

#endif // __VW_IMAGE_PIXELMATH_H__
