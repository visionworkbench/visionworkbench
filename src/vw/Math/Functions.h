// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Math/Functions.h
///
/// Mathematical functions.
///
/// This file provides implementations of standard mathematical
/// functions that may not be available on all platforms.  The
/// full set of mathematical operators is made available for
/// use in namespace vw.
///
#ifndef __VW_MATH_FUNCTIONS_H__
#define __VW_MATH_FUNCTIONS_H__

#include <cmath>
#include <cstdlib>
#include <limits>
#include <complex>

#include <vw/config.h>
#include <vw/Core/FundamentalTypes.h>

namespace vw {
namespace math {
  // These function implementations are defined in namespace impl
  // to avoid ambiguous collisions with system-defined implementations
  // on platforms that provide them.
namespace impl {

  /// A special inlinable implementation of round()
  /// for the common case of double->int32.
  inline int32 _round( double val ) {
    if( val < 0 ) return int32(val-0.5);
    else return int32(val+0.5);
  }

  /// A special inlinable implementation of round()
  /// for the common case of float->int32.
  inline int32 _round( float val ) {
    if( val < 0 ) return int32(val-0.5f);
    else return int32(val+0.5f);
  }

  /// A special inlinable implementation of floor()
  /// for the common case of double->int32.
  inline int32 _floor( double val ) {
    if( val < 0 ) {
      int32 ival = int32(val);
      if( double(ival)==val ) return ival;
      else return ival - 1;
    }
    else return int32(val);
  }

  /// A replacement for the C++ modulo operator % which sucks less.
  /// Always returns a value in the range [0,b), or (b,0] if b<0.
  template <class T>
  inline T mod( T a, T b ) {
    T m = a%b;
    return (b<0) ? ((m>0)?(m+b):(m)) : ((m<0)?(m+b):(m));
  }

  double erfc(double);

  /// The error function, i.e. the integral of the normal
  /// distribution from zero to x.  For small values of x
  /// we use a Taylor series expansion.  For large values
  /// of x we defer to erfc, below.
  inline double erf( double x ) {
    static const double two_sqrtpi= 1.128379167095513; // 2/sqrt(pi)
    if( fabs(x) > 2.0 ) return 1.0 - erfc(x);
    double xsqr=x*x, term=x, accum=x;
    int i=1;
    do {
      term *= -xsqr/i;
      accum += term / (2*i+1);
      ++i;
    } while( fabs(term/accum) > 1e-12 ); // About 40 bits of accuracy
    return two_sqrtpi * accum;
  }

  /// The complementary error function, i.e. one minus the
  /// error function.  For large values of x we use a
  /// continued fraction expansion.  For small values of x
  /// we defer to erf, above.
  inline double erfc( double x ) {
    static const double one_sqrtpi= 0.5641895835477563; // 1/sqrt(pi)
    if( fabs(x) < 2.0 ) return 1.0 - erf(x);
    if( x < 0 ) return 2.0 - erfc(-x);
    double n1=1, n2=x;                  // last two numerators
    double d1=x, d2=x*x+0.5;            // last two denominators
    double q1,   q2=n2/d2;              // last two approximations
    double coeff = 1.0;
    do {
      double tmp;
      tmp = n1*coeff + n2*x;
      n1 = n2;
      n2 = tmp;
      tmp = d1*coeff + d2*x;
      d1 = d2;
      d2 = tmp;
      coeff += 0.5;
      q1 = q2;
      q2 = n2/d2;
    } while ( fabs(q1-q2)/q2 > 1e-12 ); // About 40 bits of accuracy
    return one_sqrtpi * exp(-x*x) * q2;
  }

  // Notes for future function implementations:
  // * The asymptotic behavior of acosh(x) is ln(2x) for large x
  //   and sqrt(2(x-1)) for small x (i.e. x near 1).
  // * The asymptotic behavior of asinh(x) is also ln(2x).

  inline float hypotf (float x, float y) {
    return sqrtf(x*x + y*y);
  }
  inline double hypot (double x, double y) {
    return sqrt(x*x + y*y);
  }
  inline long double hypotl (long double x, long double y) {
    return sqrtl(x*x + y*y);
  }
} // namespace impl
} // namespace math

  using ::fabs;
  using ::acos;
  using ::asin;
  using ::atan;
  using ::cos;
  using ::sin;
  using ::tan;
  using ::cosh;
  using ::sinh;
  using ::tanh;
  using ::exp;
  using ::log;
  using ::log10;
  using ::sqrt;
  using ::ceil;
  using ::floor;
  using ::atan2;
  using ::pow;

  using math::impl::mod;

#ifdef WIN32
  using math::impl::erf;
  using math::impl::erfc;
  using math::impl::hypot;
  using math::impl::hypotl;
  using math::impl::hypotf;
#else
  using ::acosh;
  using ::asinh;
  using ::atanh;
#ifdef VW_HAVE_EXP2
  using ::exp2;
#endif
#ifdef VW_HAVE_LOG2
  using ::log2;
#endif
#ifdef VW_HAVE_TGAMMA
  using ::tgamma;
#endif
  using ::expm1;
  using ::log1p;
  using ::cbrt;
  using ::erf;
  using ::erfc;
  using ::lgamma;
  using ::round;
  using ::trunc;
  using ::hypot;
  using ::copysign;
  using ::fdim;
#endif // WIN32

#if defined(__APPLE__) && defined(__POWERPC__)
  inline long double erfl( long double arg ) {
    return erf( (double)arg );
  }

  inline long double erfcl( long double arg ) {
    return erfc( (double)arg );
  }
#endif

} // namespace vw

#endif  // __VW_MATH_FUNCTIONS_H__
