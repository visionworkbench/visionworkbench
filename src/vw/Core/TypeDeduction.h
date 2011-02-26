// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file TypeDeduction.h
///
/// Defines standard type deduction behavior.
///
/// Known issues:
///
/// The default type deduction behavior for compound types does not
/// use operation-specific specializations for the channel types.
/// Thus, for instance, SumType<A<B>,A<C>> by default reduces to
/// A<PromoteType<B,C>> rather than A<SumType<B,C>>.  This does not
/// currently impact any supported channel types.
///
#ifndef __VW_CORE_TYPE_DEDUCTION_H__
#define __VW_CORE_TYPE_DEDUCTION_H__

#include <boost/type_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/not.hpp>

#include <vw/Core/CompoundTypes.h>
#include <vw/Core/FundamentalTypes.h>

namespace vw {

  /// \cond INTERNAL

  // ********************************************************************
  // Operator Return Type Deduction Routines
  //
  // The following type-deduction metaprograms can be used to deduce the
  // return types of various operations on compound types.  This exists
  // and must be maintained until either C++ officially adopts the typeof
  // keyword (and compilers actually support it) or the BOOST folks get
  // their act together and agree on a single unified framework for this
  // sort of thing.  (In the latter case some code still needs to stay,
  // but it should get reworked to fit within said framework.)
  //
  // At the moment nothing fancy is supported.  For instance, compound
  // type sums are only supported when the two arguments have the same
  // type.  However, everything is in place to make it easy to extend
  // should the need arise.
  // ********************************************************************

  // Okay, I'll be the first to admit that this is a really brain-dead
  // way to do this.  Nevertheless, it works and will give us a standard
  // against which to compare better methods later.
#define __VW_STANDARD_ARITHMETIC_CONVERSIONS(Helper)                    \
    /* Combinations that return float128 */                          \
    template<> struct Helper<float128,float128> { typedef float128 type; }; \
    template<> struct Helper<float128,float64> { typedef float128 type; }; \
    template<> struct Helper<float64,float128> { typedef float128 type; }; \
    template<> struct Helper<float128,float> { typedef float128 type; }; \
    template<> struct Helper<float,float128> { typedef float128 type; }; \
    template<> struct Helper<float128,uint64> { typedef float128 type; }; \
    template<> struct Helper<uint64,float128> { typedef float128 type; }; \
    template<> struct Helper<float128,int64> { typedef float128 type; }; \
    template<> struct Helper<int64,float128> { typedef float128 type; }; \
    template<> struct Helper<float128,uint32> { typedef float128 type; }; \
    template<> struct Helper<uint32,float128> { typedef float128 type; }; \
    template<> struct Helper<float128,int32> { typedef float128 type; }; \
    template<> struct Helper<int32,float128> { typedef float128 type; }; \
    template<> struct Helper<float128,uint16> { typedef float128 type; }; \
    template<> struct Helper<uint16,float128> { typedef float128 type; }; \
    template<> struct Helper<float128,int16> { typedef float128 type; }; \
    template<> struct Helper<int16,float128> { typedef float128 type; }; \
    template<> struct Helper<float128,uint8> { typedef float128 type; }; \
    template<> struct Helper<uint8,float128> { typedef float128 type; }; \
    template<> struct Helper<float128,int8> { typedef float128 type; }; \
    template<> struct Helper<int8,float128> { typedef float128 type; }; \
    /* Combinations that return float64 */                               \
    template<> struct Helper<float64,float64> { typedef float64 type; }; \
    template<> struct Helper<float64,float> { typedef float64 type; }; \
    template<> struct Helper<float,float64> { typedef float64 type; }; \
    template<> struct Helper<float64,uint64> { typedef float64 type; }; \
    template<> struct Helper<uint64,float64> { typedef float64 type; }; \
    template<> struct Helper<float64,int64> { typedef float64 type; }; \
    template<> struct Helper<int64,float64> { typedef float64 type; }; \
    template<> struct Helper<float64,uint32> { typedef float64 type; }; \
    template<> struct Helper<uint32,float64> { typedef float64 type; }; \
    template<> struct Helper<float64,int32> { typedef float64 type; }; \
    template<> struct Helper<int32,float64> { typedef float64 type; }; \
    template<> struct Helper<float64,uint16> { typedef float64 type; }; \
    template<> struct Helper<uint16,float64> { typedef float64 type; }; \
    template<> struct Helper<float64,int16> { typedef float64 type; }; \
    template<> struct Helper<int16,float64> { typedef float64 type; }; \
    template<> struct Helper<float64,uint8> { typedef float64 type; }; \
    template<> struct Helper<uint8,float64> { typedef float64 type; }; \
    template<> struct Helper<float64,int8> { typedef float64 type; }; \
    template<> struct Helper<int8,float64> { typedef float64 type; }; \
    /* Combinations that return float */                                \
    template<> struct Helper<float,float> { typedef float type; };      \
    template<> struct Helper<uint64,float> { typedef float type; }; \
    template<> struct Helper<float,int64> { typedef float type; }; \
    template<> struct Helper<int64,float> { typedef float type; }; \
    template<> struct Helper<float,uint32> { typedef float type; }; \
    template<> struct Helper<uint32,float> { typedef float type; }; \
    template<> struct Helper<float,int32> { typedef float type; }; \
    template<> struct Helper<int32,float> { typedef float type; }; \
    template<> struct Helper<float,uint16> { typedef float type; }; \
    template<> struct Helper<uint16,float> { typedef float type; }; \
    template<> struct Helper<float,int16> { typedef float type; }; \
    template<> struct Helper<int16,float> { typedef float type; }; \
    template<> struct Helper<float,uint8> { typedef float type; }; \
    template<> struct Helper<uint8,float> { typedef float type; }; \
    template<> struct Helper<float,int8> { typedef float type; }; \
    template<> struct Helper<int8,float> { typedef float type; }; \
    /* Combinations that return uint64 */                        \
    template<> struct Helper<uint64,uint64> { typedef uint64 type; }; \
    template<> struct Helper<uint64,int64> { typedef uint64 type; }; \
    template<> struct Helper<int64,uint64> { typedef uint64 type; }; \
    template<> struct Helper<uint64,uint32> { typedef uint64 type; }; \
    template<> struct Helper<uint32,uint64> { typedef uint64 type; }; \
    template<> struct Helper<uint64,int32> { typedef uint64 type; }; \
    template<> struct Helper<int32,uint64> { typedef uint64 type; }; \
    template<> struct Helper<uint64,uint16> { typedef uint64 type; }; \
    template<> struct Helper<uint16,uint64> { typedef uint64 type; }; \
    template<> struct Helper<uint64,int16> { typedef uint64 type; }; \
    template<> struct Helper<int16,uint64> { typedef uint64 type; }; \
    template<> struct Helper<uint64,uint8> { typedef uint64 type; }; \
    template<> struct Helper<uint8,uint64> { typedef uint64 type; }; \
    template<> struct Helper<uint64,int8> { typedef uint64 type; }; \
    template<> struct Helper<int8,uint64> { typedef uint64 type; }; \
    /* Combinations that return int64 */                          \
    template<> struct Helper<int64,int64> { typedef int64 type; }; \
    template<> struct Helper<int64,uint32> { typedef int64 type; }; \
    template<> struct Helper<uint32,int64> { typedef int64 type; }; \
    template<> struct Helper<int64,int32> { typedef int64 type; }; \
    template<> struct Helper<int32,int64> { typedef int64 type; }; \
    template<> struct Helper<int64,uint16> { typedef int64 type; }; \
    template<> struct Helper<uint16,int64> { typedef int64 type; }; \
    template<> struct Helper<int64,int16> { typedef int64 type; }; \
    template<> struct Helper<int16,int64> { typedef int64 type; }; \
    template<> struct Helper<int64,uint8> { typedef int64 type; }; \
    template<> struct Helper<uint8,int64> { typedef int64 type; }; \
    template<> struct Helper<int64,int8> { typedef int64 type; }; \
    template<> struct Helper<int8,int64> { typedef int64 type; }; \
    /* Combinations that return uint32 */                         \
    template<> struct Helper<uint32,uint32> { typedef uint32 type; }; \
    template<> struct Helper<uint32,int32> { typedef uint32 type; }; \
    template<> struct Helper<int32,uint32> { typedef uint32 type; }; \
    template<> struct Helper<uint32,uint16> { typedef uint32 type; }; \
    template<> struct Helper<uint16,uint32> { typedef uint32 type; }; \
    template<> struct Helper<uint32,int16> { typedef uint32 type; }; \
    template<> struct Helper<int16,uint32> { typedef uint32 type; }; \
    template<> struct Helper<uint32,uint8> { typedef uint32 type; }; \
    template<> struct Helper<uint8,uint32> { typedef uint32 type; }; \
    template<> struct Helper<uint32,int8> { typedef uint32 type; }; \
    template<> struct Helper<int8,uint32> { typedef uint32 type; }; \
    /* Combinations that return int32 */                           \
    template<> struct Helper<int32,int32> { typedef int32 type; }; \
    template<> struct Helper<int32,uint16> { typedef int32 type; }; \
    template<> struct Helper<uint16,int32> { typedef int32 type; }; \
    template<> struct Helper<int32,int16> { typedef int32 type; }; \
    template<> struct Helper<int16,int32> { typedef int32 type; }; \
    template<> struct Helper<int32,uint8> { typedef int32 type; }; \
    template<> struct Helper<uint8,int32> { typedef int32 type; }; \
    template<> struct Helper<int32,int8> { typedef int32 type; }; \
    template<> struct Helper<int8,int32> { typedef int32 type; }; \
    /* Combinations that return uint16 */                       \
    template<> struct Helper<uint16,uint16> { typedef uint16 type; }; \
    template<> struct Helper<uint16,int16> { typedef uint16 type; }; \
    template<> struct Helper<int16,uint16> { typedef uint16 type; }; \
    template<> struct Helper<uint16,uint8> { typedef uint16 type; }; \
    template<> struct Helper<uint8,uint16> { typedef uint16 type; }; \
    template<> struct Helper<uint16,int8> { typedef uint16 type; }; \
    template<> struct Helper<int8,uint16> { typedef uint16 type; }; \
    /* Combinations that return int16 */                       \
    template<> struct Helper<int16,int16> { typedef int16 type; }; \
    template<> struct Helper<int16,uint8> { typedef int16 type; }; \
    template<> struct Helper<uint8,int16> { typedef int16 type; }; \
    template<> struct Helper<int16,int8> { typedef int16 type; }; \
    template<> struct Helper<int8,int16> { typedef int16 type; }; \
    /* Combinations that return uint8 */                        \
    template<> struct Helper<uint8,uint8> { typedef uint8 type; }; \
    template<> struct Helper<uint8,int8> { typedef uint8 type; }; \
    template<> struct Helper<int8,uint8> { typedef uint8 type; }; \
    /* Combinations that return int8 */                       \
    template<> struct Helper<int8,int8> { typedef int8 type; }; \
    /* Combinations with a user type return that type */                \
    template<class T> struct Helper<float128,T> { typedef T type; }; \
    template<class T> struct Helper<T,float128> { typedef T type; }; \
    template<class T> struct Helper<float64,T> { typedef T type; }; \
    template<class T> struct Helper<T,float64> { typedef T type; }; \
    template<class T> struct Helper<float,T> { typedef T type; }; \
    template<class T> struct Helper<T,float> { typedef T type; }; \
    template<class T> struct Helper<uint64,T> { typedef T type; }; \
    template<class T> struct Helper<T,uint64> { typedef T type; }; \
    template<class T> struct Helper<int64,T> { typedef T type; }; \
    template<class T> struct Helper<T,int64> { typedef T type; }; \
    template<class T> struct Helper<uint32,T> { typedef T type; }; \
    template<class T> struct Helper<T,uint32> { typedef T type; }; \
    template<class T> struct Helper<int32,T> { typedef T type; }; \
    template<class T> struct Helper<T,int32> { typedef T type; }; \
    template<class T> struct Helper<uint16,T> { typedef T type; }; \
    template<class T> struct Helper<T,uint16> { typedef T type; }; \
    template<class T> struct Helper<int16,T> { typedef T type; }; \
    template<class T> struct Helper<T,int16> { typedef T type; }; \
    template<class T> struct Helper<uint8,T> { typedef T type; }; \
    template<class T> struct Helper<T,uint8> { typedef T type; }; \
    template<class T> struct Helper<int8,T> { typedef T type; }; \
    template<class T> struct Helper<T,int8> { typedef T type; };

  template <class T>
  struct TypeDeductionError;

  template <class Arg1T, class Arg2T>
  struct TypeDeductionHelper {
    typedef TypeDeductionError<TypeDeductionHelper> type;
  };

  template <class ArgT>
  struct TypeDeductionHelper<ArgT, ArgT> {
    typedef ArgT type;
  };

  __VW_STANDARD_ARITHMETIC_CONVERSIONS(TypeDeductionHelper)


  // ********************************************************************
  // Now we set up the default promotion behavior for general types.
  // This includes promotions of compound types (e.g. pixel types)
  // when used in operations with fundamental types.  Thus for example
  // multiplying PixelRGB<int> by float should return PixelRGB<float>.
  // We address this by recursing on the template parameter.  Users
  // can explicitly specialize PromoteTypeSpecialization<> to support
  // custom behaviors when both arguments are custom types.
  // ********************************************************************

  // First we forward-declare the specialization class.
  template <class Arg1T, class Arg2T> class PromoteTypeSpecialization;

  // This helper class forwards the default case to TypeDeductionHelper.
  template <class Arg1T, class Arg2T, bool Arg1IsCompound, bool Arg2IsCompound>
  struct PromoteTypeSpecializationHelper : public TypeDeductionHelper<Arg1T,Arg2T> {};

  // Recursive behavior when Arg1T is a compound type.
  template <class Arg1T, class Arg2T>
  struct PromoteTypeSpecializationHelper<Arg1T,Arg2T,true,false> {
    typedef typename CompoundChannelCast<Arg1T,typename PromoteTypeSpecialization<typename CompoundChannelType<Arg1T>::type, Arg2T>::type>::type type;
  };

  // Recursive behavior when Arg2T is a compound type.
  template <class Arg1T, class Arg2T>
  struct PromoteTypeSpecializationHelper<Arg1T,Arg2T,false,true> {
    typedef typename CompoundChannelCast<Arg2T,typename PromoteTypeSpecialization<Arg1T, typename CompoundChannelType<Arg2T>::type>::type>::type type;
  };

  // Recursive behavior when both Arg1T and Arg2T are compound types.
  template <class Arg1T, class Arg2T>
  struct PromoteTypeSpecializationHelper<Arg1T,Arg2T,true,true> {
    typedef typename boost::is_same<Arg1T,typename CompoundChannelCast<Arg2T,typename CompoundChannelType<Arg1T>::type>::type>::type if_same_type;
    typedef typename CompoundChannelCast<Arg1T,typename PromoteTypeSpecialization<typename CompoundChannelType<Arg1T>::type,typename CompoundChannelType<Arg2T>::type>::type>::type is_same_type;
    typedef typename TypeDeductionHelper<Arg1T,Arg2T>::type not_same_type;
    typedef typename boost::mpl::if_<if_same_type,is_same_type,not_same_type>::type type;
  };

  // Dispatch to the appropriate helper based on which if any
  // argument is a compound type.
  template <class Arg1T, class Arg2T>
  class PromoteTypeSpecialization
    : public PromoteTypeSpecializationHelper<Arg1T, Arg2T, IsCompound<Arg1T>::value, IsCompound<Arg2T>::value>
  {};


  // ********************************************************************
  // Finally, we provide forwarding classes for the different operations
  // that the user can specifically specialize in special cases if
  // needed.  The classes the user *uses*, such as SumType, simply
  // strip off any CV-qualification and forward to the appropriate
  // specialization type.  If the user needs to override the type
  // deduction behavior, they should generally do so by overriding one
  // of the ...Specialization classes, so that they don't have to worry
  // about CV-qualification themselves.
  // ********************************************************************

  template <class Arg1T, class Arg2T>
  struct SumTypeSpecialization : public PromoteTypeSpecialization<Arg1T,Arg2T> {};

  template <class Arg1T, class Arg2T>
  struct DifferenceTypeSpecialization : public PromoteTypeSpecialization<Arg1T,Arg2T> {};

  template <class Arg1T, class Arg2T>
  struct ProductTypeSpecialization : public PromoteTypeSpecialization<Arg1T,Arg2T> {};

  template <class Arg1T, class Arg2T>
  struct QuotientTypeSpecialization : public PromoteTypeSpecialization<Arg1T,Arg2T> {};


  template <class Arg1T, class Arg2T>
  struct PromoteType {
    typedef typename PromoteTypeSpecialization<typename boost::remove_cv<Arg1T>::type,
                                               typename boost::remove_cv<Arg2T>::type>::type type;
  };

  template <class Arg1T, class Arg2T>
  struct SumType {
    typedef typename SumTypeSpecialization<typename boost::remove_cv<Arg1T>::type,
                                           typename boost::remove_cv<Arg2T>::type>::type type;
  };

  template <class Arg1T, class Arg2T>
  struct DifferenceType {
    typedef typename DifferenceTypeSpecialization<typename boost::remove_cv<Arg1T>::type,
                                                  typename boost::remove_cv<Arg2T>::type>::type type;
  };

  template <class Arg1T, class Arg2T>
  struct ProductType {
    typedef typename ProductTypeSpecialization<typename boost::remove_cv<Arg1T>::type,
                                               typename boost::remove_cv<Arg2T>::type>::type type;
  };

  template <class Arg1T, class Arg2T>
  struct QuotientType {
    typedef typename QuotientTypeSpecialization<typename boost::remove_cv<Arg1T>::type,
                                                typename boost::remove_cv<Arg2T>::type>::type type;
  };

  /// \endcond

} // namespace vw

#endif // __VW_TYPE_DEDUCTION_H__
