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
#include <vector>
#include <algorithm>

#include <vw/vw_config.h>
#include <vw/Core/CompoundTypes.h>
#include <vw/Core/TypeDeduction.h>
#include <vw/Core/Functors.h>
#include <vw/Core/Exception.h>

#include <queue>

// The math.h header in FreeBSD (and possibly other platforms) does not
// include routines for manipulating long doubles.  We disable long
// double VW math routines here for certain platforms.
#if defined(__FreeBSD__) || defined(_WIN32)
#define __VW_MATH_DISABLE_LONG_DOUBLE_ARITHMETIC
#endif

namespace vw {
  namespace math {

    template <class T> struct StdMathType              { typedef double type; };
    template <>        struct StdMathType<float      > { typedef float type; };
    template <>        struct StdMathType<long double> { typedef long double type; };
    template <class T1, class T2> struct StdMathType2 {
      typedef typename StdMathType<typename PromoteType<T1,T2>::type>::type type;
    };

    template <class FuncT, class ArgT, bool ArgIsCompound>
    struct ArgUnaryFunctorTypeHelper : public StdMathType<ArgT> {};

    template <class FuncT, class ArgT>
    struct ArgUnaryFunctorTypeHelper<FuncT,ArgT,true> {
      typedef typename CompoundChannelType<ArgT>::type channel_type;
      typedef typename ArgUnaryFunctorTypeHelper<FuncT,channel_type,IsCompound<channel_type>::value>::type result_channel_type;
      typedef typename CompoundChannelCast<ArgT,result_channel_type>::type type;
    };

    template <class F, class T>
    struct ArgUnaryFunctorType : public ArgUnaryFunctorTypeHelper<F,T,IsCompound<T>::value> {};


    template <class FuncT, class Arg1T, class Arg2T, bool Arg1IsCompound, bool Arg2IsCompound>
    struct BinaryFunctorTypeHelper : public StdMathType2<Arg1T,Arg2T> {};

    template <class FuncT, class Arg1T, class Arg2T>
    struct BinaryFunctorTypeHelper<FuncT,Arg1T,Arg2T,true,false> {
      typedef typename CompoundChannelType<Arg1T>::type channel_type;
      typedef typename BinaryFunctorTypeHelper<FuncT,channel_type,Arg2T,IsCompound<channel_type>::value,false>::type result_channel_type;
      typedef typename CompoundChannelCast<Arg1T,result_channel_type>::type type;
    };

    template <class FuncT, class Arg1T, class Arg2T>
    struct BinaryFunctorTypeHelper<FuncT,Arg1T,Arg2T,false,true> {
      typedef typename CompoundChannelType<Arg2T>::type channel_type;
      typedef typename BinaryFunctorTypeHelper<FuncT,Arg1T,channel_type,false,IsCompound<channel_type>::value>::type result_channel_type;
      typedef typename CompoundChannelCast<Arg2T,result_channel_type>::type type;
    };

    template <class FuncT, class Arg1T, class Arg2T>
    struct BinaryFunctorTypeHelper<FuncT,Arg1T,Arg2T,true,true> {
      typedef typename CompoundChannelType<Arg1T>::type channel1_type;
      typedef typename CompoundChannelType<Arg2T>::type channel2_type;
      typedef typename BinaryFunctorTypeHelper<FuncT,channel1_type,channel2_type,IsCompound<channel1_type>::value,IsCompound<channel2_type>::value>::type result_channel_type;
      typedef typename boost::mpl::if_<CompoundIsCompatible<Arg1T,Arg2T>, typename CompoundChannelCast<Arg1T,result_channel_type>::type, TypeDeductionError<BinaryFunctorTypeHelper> >::type type;
    };

    template <class FuncT, class Arg1T, class Arg2T>
    struct BinaryFunctorType : public BinaryFunctorTypeHelper<FuncT,Arg1T,Arg2T,IsCompound<Arg1T>::value,IsCompound<Arg2T>::value> {};

    // ********************************************************************
    // Standard Mathematical Functors
    // ********************************************************************

#if !defined(__VW_MATH_DISABLE_LONG_DOUBLE_ARITHMETIC)
    // This will actually create a specialization for long doubles
#define __VW_MATH_UNARY_LONG_DOUBLE_IMPL(func)          \
    long double operator()( long double arg ) const {   \
      return func##l(arg);                              \
    }
#define __VW_MATH_BINARY_LONG_DOUBLE_IMPL(func)                         \
    long double operator()( long double arg1, long double arg2 ) const { \
      return func##l(arg1,arg2);                                        \
    }
#else
    // This is basically a no-op
#define __VW_MATH_UNARY_LONG_DOUBLE_IMPL(func)
#define __VW_MATH_BINARY_LONG_DOUBLE_IMPL(func)
#endif

// Macros to easily create functors with a certain interface
#define __VW_UNARY_MATH_FUNCTOR(name,func)                              \
    struct Arg##name##Functor                                           \
      : UnaryReturnBinaryTemplateBind1st<ArgUnaryFunctorType,Arg##name##Functor> { \
        template <class ValT>                                           \
          typename ArgUnaryFunctorType<Arg##name##Functor,ValT>::type   \
          operator()( ValT arg ) const {                                \
          return func(arg);                                             \
        }                                                               \
        float operator()( float arg ) const {                           \
          return func##f(arg);                                          \
        }                                                               \
        __VW_MATH_UNARY_LONG_DOUBLE_IMPL(func)                          \
          };                                                            \
    using ::func;
// END __VW_UNARY_MATH_FUNCTOR
#define __VW_BINARY_MATH_FUNCTOR(name,func)                             \
    struct ArgArg##name##Functor                                        \
      : BinaryReturnTernaryTemplateBind1st<BinaryFunctorType,ArgArg##name##Functor> { \
    float operator()( float arg1, float arg2 ) const {                  \
    return func##f(arg1,arg2);                                          \
  }                                                                     \
    __VW_MATH_BINARY_LONG_DOUBLE_IMPL(func)                             \
    template <class Arg1T, class Arg2T>                                 \
    typename BinaryFunctorType<ArgArg##name##Functor,Arg1T,Arg2T>::type \
    inline operator()( Arg1T const& arg1, Arg2T const& arg2 ) const {   \
    return func( arg1, arg2 );                                          \
  }                                                                     \
  };                                                                    \
    template <class ValT>                                               \
    struct ArgVal##name##Functor                                        \
      : UnaryReturnTernaryTemplateBind1st3rd<BinaryFunctorType,ArgArg##name##Functor,ValT> { \
    ValT m_val;                                                         \
    ArgVal##name##Functor(ValT val) : m_val(val) {}                     \
    template <class ArgT>                                               \
    typename BinaryFunctorType<ArgArg##name##Functor,ArgT,ValT>::type   \
    inline operator()( ArgT const& arg ) const {                        \
    return ArgArg##name##Functor()( arg, m_val );                       \
  }                                                                     \
  };                                                                    \
    template <class ValT>                                               \
    struct ValArg##name##Functor                                        \
      : UnaryReturnTernaryTemplateBind1st2nd<BinaryFunctorType,ArgArg##name##Functor,ValT> { \
    ValT m_val;                                                         \
    ValArg##name##Functor(ValT val) : m_val(val) {}                     \
    template <class ArgT>                                               \
    typename BinaryFunctorType<ArgArg##name##Functor,ValT,ArgT>::type   \
    inline operator()( ArgT const& arg ) const {                        \
    return ArgArg##name##Functor()( m_val, arg );                       \
  }                                                                     \
  };                                                                    \
    using ::func;
// END __VW_BINARY_MATH_FUNCTOR

    __VW_UNARY_MATH_FUNCTOR( Fabs, fabs )
    __VW_UNARY_MATH_FUNCTOR( Acos, acos )
    __VW_UNARY_MATH_FUNCTOR( Asin, asin )
    __VW_UNARY_MATH_FUNCTOR( Atan, atan )
    __VW_UNARY_MATH_FUNCTOR( Cos, cos )
    __VW_UNARY_MATH_FUNCTOR( Sin, sin )
    __VW_UNARY_MATH_FUNCTOR( Tan, tan )
    __VW_UNARY_MATH_FUNCTOR( Cosh, cosh )
    __VW_UNARY_MATH_FUNCTOR( Sinh, sinh )
    __VW_UNARY_MATH_FUNCTOR( Tanh, tanh )
    __VW_UNARY_MATH_FUNCTOR( Exp, exp )
    __VW_UNARY_MATH_FUNCTOR( Log, log )
    __VW_UNARY_MATH_FUNCTOR( Log10, log10 )
    __VW_UNARY_MATH_FUNCTOR( Sqrt, sqrt )
    __VW_UNARY_MATH_FUNCTOR( Ceil, ceil )
    __VW_UNARY_MATH_FUNCTOR( Floor, floor )

    __VW_BINARY_MATH_FUNCTOR( Atan2, atan2 )
    __VW_BINARY_MATH_FUNCTOR( Pow, pow )
    __VW_BINARY_MATH_FUNCTOR( Hypot, hypot )

#ifndef WIN32
    __VW_UNARY_MATH_FUNCTOR( Acosh, acosh )
    __VW_UNARY_MATH_FUNCTOR( Asinh, asinh )
    __VW_UNARY_MATH_FUNCTOR( Atanh, atanh )

#ifdef VW_HAVE_EXP2
    __VW_UNARY_MATH_FUNCTOR( Exp2, exp2 )
#endif
#ifdef VW_HAVE_LOG2
    __VW_UNARY_MATH_FUNCTOR( Log2, log2 )
#endif
#ifdef VW_HAVE_TGAMMA
    __VW_UNARY_MATH_FUNCTOR( Tgamma, tgamma )
#endif
    __VW_UNARY_MATH_FUNCTOR( Lgamma, lgamma )
    __VW_UNARY_MATH_FUNCTOR( Expm1, expm1 )
    __VW_UNARY_MATH_FUNCTOR( Log1p, log1p )
    __VW_UNARY_MATH_FUNCTOR( Cbrt, cbrt )
    __VW_UNARY_MATH_FUNCTOR( Erf, erf )
    __VW_UNARY_MATH_FUNCTOR( Erfc, erfc )
    __VW_UNARY_MATH_FUNCTOR( Round, round )
      __VW_UNARY_MATH_FUNCTOR( Trunc, trunc )

      __VW_BINARY_MATH_FUNCTOR( Copysign, copysign )
    __VW_BINARY_MATH_FUNCTOR( Fdim, fdim )
#endif

// Clean up the macros we are finished using
#undef __VW_UNARY_MATH_FUNCTOR
#undef __VW_BINARY_MATH_FUNCTOR


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
      ValT operator()( ValT const& /*val*/ ) const {
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
      template <class FuncT> struct result<FuncT(int32)> { typedef int32 type; };
      template <class FuncT> struct result<FuncT(int64)> { typedef int64 type; };
      template <class FuncT, class ValT> struct result<FuncT(std::complex<ValT>)> { typedef ValT type; };

      template <class ValT>
      typename result<ArgAbsFunctor(ValT)>::type
      inline operator()( ValT val ) const {
        return DefaultAbsBehavior<std::numeric_limits<ValT>::is_integer>::apply(val);
      }

      inline float operator()( float val ) const { return ::fabsf(val); }
      inline long operator()( long val ) const { return std::labs(val); }
#ifdef VW_HAVE_FABSL
      inline long double operator()( long double val ) const { return ::fabsl(val); }
#endif
#ifdef VW_HAVE_LLABS
      inline long long operator()( long long val ) const { return ::llabs(val); }
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

    // Square Functor (so we don't have to always invoke POW)
    struct ArgSquareFunctor : UnaryReturnSameType {
      template <class ValT>
      ValT operator()( ValT const& val ) const {
        return val*val;
      }
    };

    // General-purpose accumulation functor
    template <class AccumT, class FuncT = ArgArgInPlaceSumFunctor>
    struct Accumulator : ReturnFixedType<AccumT const&> {
    private:
      AccumT m_accum;
      FuncT m_func;
    public:
      typedef AccumT value_type;

      Accumulator() : m_accum(), m_func() {}
      Accumulator( FuncT const& func ) : m_accum(), m_func(func) {}
      Accumulator( AccumT const& accum ) : m_accum(accum), m_func() {}
      Accumulator( AccumT const& accum, FuncT const& func ) : m_accum(accum), m_func(func) {}

      template <class ArgT>
      inline AccumT const& operator()( ArgT const& arg ) {
        m_func(m_accum, arg);
        return m_accum;
      }

      inline AccumT const& value() const {
        return m_accum;
      }

      void reset( AccumT const& accum = AccumT() ) {
        m_accum = accum;
      }
    };


    /// Computes minimum and maximum values
    template <class ValT>
    class MinMaxAccumulator : public ReturnFixedType<void> {
      ValT m_minval, m_maxval;
      bool m_valid;
    public:
      typedef std::pair<ValT,ValT> value_type;

      MinMaxAccumulator() : m_minval(0), m_maxval(0), m_valid(false) {}

      void operator()( ValT const& arg ) {
        if ( ! m_valid ) {
          m_minval = m_maxval = arg;
          m_valid = true;
        }
        else {
          if( arg < m_minval ) m_minval = arg;
          if( m_maxval < arg ) m_maxval = arg;
        }
      }

      ValT minimum() const {
        VW_ASSERT(m_valid, ArgumentErr() << "MinMaxAccumulator: no valid samples");
        return m_minval;
      }

      ValT maximum() const {
        VW_ASSERT(m_valid, ArgumentErr() << "MinMaxAccumulator: no valid samples");
        return m_maxval;
      }

      std::pair<ValT,ValT> value() const {
        VW_ASSERT(m_valid, ArgumentErr() << "MinMaxAccumulator: no valid samples");
        return std::make_pair(m_minval,m_maxval);
      }
    };

    // Note: This function modifies the input!
    // Note that we always return double, for precision.
    template <class T>
    double destructive_median(std::vector<T> & vec) {
      int len = vec.size();
      VW_ASSERT(len, ArgumentErr() << "median: no valid samples.");
      std::sort(vec.begin(), vec.end());
      return len%2 ? vec[len/2] : (vec[len/2 - 1] + vec[len/2]) / 2.0;
    }

    // Compute the normalized median absolute deviation:
    // nmad = 1.4826 * median(abs(X - median(X)))  
    // Note: This function modifies the input!
    // This always return double, for precision.
    template <class T>
    double destructive_nmad(std::vector<T> & vec){
      int len = vec.size();
      VW_ASSERT(len, ArgumentErr() << "nmad: no valid samples.");
      
      // Find the median. This sorts the vector, but that is not a problem.
      double median = destructive_median(vec);
      
      // It is safer to make a copy of the vector specifically for the
      // mad calculation, so that not to mess up the sorted values in vec.
      std::vector<double> abs_diffs(len);
      for (size_t it = 0; it < vec.size(); it++)
        abs_diffs[it] = std::abs(vec[it] - median);
      
      return 1.4826 * destructive_median(abs_diffs);
    }
  
    // Compute the percentile using
    // https://en.wikipedia.org/wiki/Percentile#The_nearest-rank_method
    // Note: This function modifies the input!
    template <class T>
    T destructive_percentile(std::vector<T> & vec, double percentile){
      
      int len = vec.size();
      VW_ASSERT(len > 0, ArgumentErr() << "percentile: no valid samples.");
      VW_ASSERT(percentile >= 0 && percentile <= 100.0,
                ArgumentErr() << "Percentile must be between 0 and 100.");

      // Sorting is vital
      std::sort(vec.begin(), vec.end());

      int index = ceil((percentile/100.0) * double(len));

      // Account for the fact that in C++ indices start from 0 
      index--;

      if (index < 0) index = 0;
      if (index >= len) index = len-1;
        
      return vec[index];
    }

    // Computes the median of the values to which it is applied.
    template <class ValT>
    class MedianAccumulator : public ReturnFixedType<void> {
      std::vector<ValT> m_values;
    public:
      typedef ValT value_type;

      void operator()( ValT const& value ) {
        m_values.push_back( value );
      }

      // This is to check if there are any values
      size_t size() {
        return m_values.size();
      }
      
      double value() {
        return destructive_median(m_values);
      }
    };

    /// Computes the mean of the values to which it is applied.
    template <class ValT>
    class MeanAccumulator : public ReturnFixedType<void> {
      typedef typename CompoundChannelCast<ValT,double>::type accum_type;
      accum_type m_accum;
      double m_count;
    public:
      typedef accum_type value_type;

      MeanAccumulator() : m_accum(), m_count() {}

      void operator()( ValT const& value ) {
        m_accum += value;
        m_count += 1.0;
      }

      value_type value() const {
        VW_ASSERT(m_count, ArgumentErr() << "MeanAccumulator: no valid samples");
        return m_accum / m_count;
      }
    };


    /// Computes the standard deviation of the values to which it is applied.
    /// - This implementation normalizes by num_samples, not num_samples - 1.
    template <class ValT>
    class StdDevAccumulator : public ReturnFixedType<void> {
      typedef typename CompoundChannelCast<ValT,double>::type accum_type;
      accum_type mom1_accum, mom2_accum;
      double num_samples;
    public:
      typedef accum_type value_type;

      StdDevAccumulator() : mom1_accum(), mom2_accum(), num_samples() {}

      void operator()( ValT const& arg ) {
        mom1_accum += arg;
        mom2_accum += (accum_type) arg * arg;
        num_samples += 1.0;
      }

      /// Return the standard deviation
      value_type value() const {
        VW_ASSERT(num_samples, ArgumentErr() << "StdDevAccumulator(): no valid samples.");
        return sqrt(mom2_accum/num_samples - (mom1_accum/num_samples)*(mom1_accum/num_samples));
      }
      /// Return the mean
      value_type mean() const {
        VW_ASSERT(num_samples, ArgumentErr() << "StdDevAccumulator(): no valid samples.");
        return mom1_accum / num_samples;
      }
    };


    /// Compute the standard deviation of values in a fixed size list.
    /// - When a new value is added, the oldest value is removed from the statistics.
    /// - This implementation normalizes by num_samples, not num_samples - 1.
    class StdDevSlidingFunctor {
    public:
      /// Constructor set with the sliding window size.
      StdDevSlidingFunctor(const size_t max_size)
        : m_max_size(max_size), m_mean(0), m_squared(0) {}

      /// Implement push with the standard functor interface
      void operator()(double new_val) {
        push(new_val);
      }

      /// Add a new value and eject the oldest value.
      void push(double new_val) {

        // If we grew larger than the size limit remove the oldest value
        if (m_values.size() >= m_max_size)
          pop();

        // Now incorporate the newest value
        m_values.push(new_val); // Record the new value
        double count = static_cast<double>(m_values.size());

        // Update the statistics to add in the new value
        double delta = new_val - m_mean;
        m_mean += delta/count;

        double delta2 = new_val - m_mean;
        m_squared += delta*delta2;
      }

      /// Remove the oldest value
      void pop() {
        if (m_values.empty()) // Handle empty case
          return;

        double old_val = m_values.front();
        double count   = static_cast<double>(m_values.size());
        m_values.pop();

        // Update the statistics to account for removing the old value
        double shrunk_count = count - 1.0;
        double new_mean = (count*m_mean - old_val)/shrunk_count;

        m_squared -= (old_val - m_mean) * (old_val - new_mean);
        m_mean = new_mean;
      };

      /// Compute the standard deviation of all current values.
      double get_std_dev() {
        double count = static_cast<double>(m_values.size());
        if (count < 2.0)
          return 0;
        return sqrt(m_squared/count);
      }

    private:
        size_t m_max_size; ///< Max number of values to store.
        std::queue<double> m_values; ///< Store current values
        double m_mean;    ///< Store the mean
        double m_squared; ///< Store sum of differences squared
    }; // End class StdDevSlidingFunctor

  } // namespace math

  // I'm not even really sure why the math namespace exists anymore.... -MDH
  using math::Accumulator;
  using math::MinMaxAccumulator;
  using math::MedianAccumulator;
  using math::MeanAccumulator;
  using math::StdDevAccumulator;
  using math::StdDevSlidingFunctor;

} // namespace vw

#endif  // __VW_MATH_FUNCTORS_H__
