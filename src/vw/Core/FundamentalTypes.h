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

/// \file FundamentalTypes.h
/// 
/// Types and traits for fundamental floating point and integral types.
///
#ifndef __VW_CORE_FUNDAMENTAL_TYPES_H__
#define __VW_CORE_FUNDAMENTAL_TYPES_H__

#include <complex>
#include <boost/cstdint.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/not.hpp>

namespace vw {

  /// Basic signed integer types
  typedef boost::int8_t  int8;
  typedef boost::int16_t int16;
  typedef boost::int32_t int32;
  typedef boost::int64_t int64;

  /// Basic unsigned integer types
  typedef boost::uint8_t  uint8;
  typedef boost::uint16_t uint16;
  typedef boost::uint32_t uint32;
  typedef boost::uint64_t uint64;

  /// Basic floating-point types
  typedef float float32;
  typedef double float64;

  /// Given a type, these traits classes identify whether or not the
  /// type is a scalar (in the mathematical sense of the word.)  This
  /// includes the built-in arithmetic types as well as complex
  /// numbers.
  template <class T> struct IsScalar : public boost::is_arithmetic<T> {};
  template <class T> struct IsScalar<std::complex<T> > : public boost::true_type {};
  template <class T> struct IsScalar<const T> : public IsScalar<T> {};


  /// Given a type, these traits classes help to determine a suitable
  /// working type for accumulation operations or other intermediate 
  /// results that require computational headroom.
  template <class T> struct AccumulatorType {};
  template <> struct AccumulatorType<bool>        { typedef int         type; };
  template <> struct AccumulatorType<vw::uint8>   { typedef vw::uint32  type; };
  template <> struct AccumulatorType<vw::int8>    { typedef vw::int32   type; };
  template <> struct AccumulatorType<vw::uint16>  { typedef vw::uint32  type; };
  template <> struct AccumulatorType<vw::int16>   { typedef vw::int32   type; };
  template <> struct AccumulatorType<vw::uint32>  { typedef vw::uint64  type; };
  template <> struct AccumulatorType<vw::int32>   { typedef vw::int64   type; };
  template <> struct AccumulatorType<vw::uint64>  { typedef vw::uint64  type; };
  template <> struct AccumulatorType<vw::int64>   { typedef vw::int64   type; };
  template <> struct AccumulatorType<vw::float32> { typedef vw::float64 type; };
  template <> struct AccumulatorType<vw::float64> { typedef vw::float64 type; };
  template <class T> struct AccumulatorType<std::complex<T> > {
    typedef std::complex<typename AccumulatorType<T>::type> type;
  };

  /// This type computation class template computes a complex type 
  /// with the same storage type as the argument.  (This is the 
  /// identity operation for complex types.)
  template <class T> struct MakeComplex { typedef std::complex<T> type; };
  template <class T> struct MakeComplex<std::complex<T> > { typedef std::complex<T> type; };

  /// This type computation class template computes the real type 
  /// with the same storage type as the (possibly complex) argument.
  /// (This is the identity operation for non-complex types.)
  template <class T> struct MakeReal { typedef T type; };
  template <class T> struct MakeReal<std::complex<T> > { typedef T type; };

  /// This function is used to work around the fact that sending 
  /// a character type, the usual container for an 8-bit numeric 
  /// type, to a C++ output stream results in the character being 
  /// printed rather than the corresponding number.
  template <class T> inline T _numeric( T v ) { return v; }
  inline unsigned _numeric( uint8 v ) { return v; }
  inline int _numeric( int8 v ) { return v; }

  /// This traits class determines whether a dereferenceable or
  /// indexable type, such as an iterator or image view, returns
  /// elements whose address can be taken.  For example, a
  /// MemoryStridingPixelIterator is referenceable whereas a
  /// ProceduralPixelIterator is not: the latter's pixel values are
  /// determined on the fly via function evaluation and do not reside
  /// anywhere in memory.  This trait evaluates to false by default,
  /// which is always safe, but certain algorithms can make
  /// optimizations when it is true.
  template <class T>
  struct IsReferenceable : public boost::false_type {};

} // namespace vw

#endif // __VW_CORE_FUNDAMENTAL_TYPES_H__
