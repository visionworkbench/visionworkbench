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

#include <math.h>
#include <cstdlib>
#include <limits>
#include <complex>

namespace vw {
namespace math {
  // These function implementations are defined in namespace impl 
  // to avoid ambiguous collisions with system-defined implementations 
  // on platforms that provide them.
namespace impl {

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

#ifdef WIN32
  using math::impl::erf;
  using math::impl::erfc;
#else
  using ::acosh;
  using ::asinh;
  using ::atanh;
#ifndef __FreeBSD__
  using ::exp2;
  using ::log2;
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
#endif

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
