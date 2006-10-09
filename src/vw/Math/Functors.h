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

/// \file Math/Functors.h
/// 
/// Mathematical functors.
/// 
/// This file provides polymorphic functor versions of the standard 
/// mathematical functions defined in e.g. math.h.
///
#ifndef __VW_MATH_FUNCTORS_H__
#define __VW_MATH_FUNCTORS_H__

#include <cstdlib>
#include <limits>
#include <complex>

#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Functors.h>
#include <vw/Math/Functions.h>

namespace vw {
namespace math {

  template <class T> struct StdMathType { typedef double type; };
  template <> struct StdMathType<float> { typedef float type; };
  template <> struct StdMathType<long double> { typedef long double type; };
  template <class T1, class T2> struct StdMathType2 {
    typedef typename StdMathType<typename PromoteType<T1,T2>::type>::type type;
  };

  // ********************************************************************
  // Standard Mathematical Functors
  // ********************************************************************

#define __VW_UNARY_IMAGE_MATH_FUNCTOR(name,func)      \
  struct Arg##name##Functor : UnaryReturnTemplateType<StdMathType> { \
    template <class ValT>                             \
    double operator()( ValT arg ) const {             \
      return func(double(arg));                       \
    }                                                 \
    float operator()( float arg ) const {             \
      return func##f(arg);                            \
    }                                                 \
    long double operator()( long double arg ) const { \
      return func##l(arg);                            \
    }                                                 \
  };                                                  \
  using ::func;

#define __VW_BINARY_IMAGE_MATH_FUNCTOR(name,func)                       \
  struct ArgArg##name##Functor : BinaryReturnTemplateType<StdMathType2> { \
    double apply( double arg1, double arg2 ) const {                    \
      return func(arg1,arg2);                                           \
    }                                                                   \
    float apply( float arg1, float arg2 ) const {                       \
      return func##f(arg1,arg2);                                        \
    }                                                                   \
    long double apply( long double arg1, long double arg2 ) const {     \
      return func##l(arg1,arg2);                                        \
    }                                                                   \
    template <class Arg1T, class Arg2T>                                 \
    typename result<ArgArg##name##Functor(Arg1T,Arg2T)>::type           \
    inline operator()( Arg1T const& arg1, Arg2T const& arg2 ) const {   \
      typedef typename result<ArgArg##name##Functor(Arg1T,Arg2T)>::type arg_type; \
      return apply( (arg_type)arg1, (arg_type)arg2 );                   \
    }                                                                   \
  };                                                                    \
  template <class ValT>                                                 \
  struct ArgVal##name##Functor : UnaryReturnBinaryTemplateBind2nd<StdMathType2,ValT> { \
    ValT m_val;                                                         \
    ArgVal##name##Functor(ValT val) : m_val(val) {}                     \
    template <class ArgT>                                               \
    typename StdMathType2<ArgT,ValT>::type                              \
    inline operator()( ArgT const& arg ) const {                        \
      return ArgArg##name##Functor()( arg, m_val );                     \
    }                                                                   \
  };                                                                    \
  template <class ValT>                                                 \
  struct ValArg##name##Functor : UnaryReturnBinaryTemplateBind1st<StdMathType2,ValT> { \
    ValT m_val;                                                         \
    ValArg##name##Functor(ValT val) : m_val(val) {}                     \
    template <class ArgT>                                               \
    typename StdMathType2<ValT,ArgT>::type                              \
    inline operator()( ArgT const& arg ) const {                        \
      return ArgArg##name##Functor()( m_val, arg );                     \
    }                                                                   \
  };                                                                    \
  using ::func;

  __VW_UNARY_IMAGE_MATH_FUNCTOR( Fabs, fabs )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Acos, acos )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Asin, asin )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Atan, atan )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Cos, cos )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Sin, sin )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Tan, tan )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Cosh, cosh )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Sinh, sinh )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Tanh, tanh )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Exp, exp )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Log, log )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Log10, log10 )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Sqrt, sqrt )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Ceil, ceil )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Floor, floor )

  __VW_BINARY_IMAGE_MATH_FUNCTOR( Atan2, atan2 )
  __VW_BINARY_IMAGE_MATH_FUNCTOR( Pow, pow )

#ifndef WIN32
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Acosh, acosh )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Asinh, asinh )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Atanh, atanh )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Exp2, exp2 )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Expm1, expm1 )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Log2, log2 )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Log1p, log1p )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Cbrt, cbrt )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Erf, erf )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Erfc, erfc )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Tgamma, tgamma )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Lgamma, lgamma )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Round, round )
  __VW_UNARY_IMAGE_MATH_FUNCTOR( Trunc, trunc )

  __VW_BINARY_IMAGE_MATH_FUNCTOR( Hypot, hypot )
  __VW_BINARY_IMAGE_MATH_FUNCTOR( Copysign, copysign )
  __VW_BINARY_IMAGE_MATH_FUNCTOR( Fdim, fdim )
#endif

#undef __VW_UNARY_IMAGE_MATH_FUNCTOR
#undef __VW_BINARY_IMAGE_MATH_FUNCTOR


  // Real part functor
  struct ArgRealFunctor : UnaryReturnTemplateType<MakeReal> {
    template <class ValT>
    ValT operator()( ValT const& val ) const {
      return val;
    }

    template <class ValT>
    ValT operator()( std::complex<ValT> const& val ) const {
      return std::real(val);
    }
  };


  // Imaginary part functor
  struct ArgImagFunctor : UnaryReturnTemplateType<MakeReal> {
    template <class ValT>
    ValT operator()( ValT const& val ) const {
      return ValT();
    }

    template <class ValT>
    ValT operator()( std::complex<ValT> const& val ) const {
      return std::imag(val);
    }
  };


  // Absolute value functor
  // This one's tricky because we have a bunch of distinct cases 
  // for integer types, floating-point types, and complex types.
  /// \cond INTERNAL
  // This is outside ArgAbsFunctor because explicit template
  // specialization doesn't work at class scope.
  template <bool IntegralN> struct DefaultAbsBehavior { template <class ValT> static inline int apply( ValT val ) { return std::abs(val); } };
  template <> struct DefaultAbsBehavior<false> { template <class ValT> static inline double apply( ValT val ) { return fabs(val); } };
  /// \endcond
  struct ArgAbsFunctor {
    template <class Args> struct result;
      
    template <class FuncT, class ValT>
    struct result<FuncT(ValT)> {
      typedef typename boost::mpl::if_c<std::numeric_limits<ValT>::is_integer, int, double>::type type;
    };

    template <class FuncT> struct result<FuncT(float)> { typedef float type; };
    template <class FuncT> struct result<FuncT(long double)> { typedef long double type; };
    template <class FuncT> struct result<FuncT(long)> { typedef long type; };
    template <class FuncT> struct result<FuncT(long long)> { typedef long long type; };
    template <class FuncT, class ValT> struct result<FuncT(std::complex<ValT>)> { typedef ValT type; };

    template <class ValT>
    typename result<ArgAbsFunctor(ValT)>::type
    inline operator()( ValT val ) const {
      return DefaultAbsBehavior<std::numeric_limits<ValT>::is_integer>::apply(val);
    }

    inline float operator()( float val ) const { return fabsf(val); }
    inline long double operator()( long double val ) const { return fabsl(val); }
    inline long operator()( long val ) const { return std::labs(val); }
#ifndef WIN32
    inline long long operator()( long long val ) const { return std::llabs(val); }
#endif

    template <class ValT>
    inline ValT operator()( std::complex<ValT> const& val ) const {
      return std::abs(val);
    }
  };


  // Complex conjugation functor
  struct ArgConjFunctor : UnaryReturnSameType {
    template <class ValT>
    ValT operator()( ValT const& val ) const {
      return val;
    }

    template <class ValT>
    std::complex<ValT> operator()( std::complex<ValT> const& val ) const {
      return std::conj(val);
    }
  };

} // namespace math
} // namespace vw

#endif  // __VW_MATH_FUNCTORS_H__
