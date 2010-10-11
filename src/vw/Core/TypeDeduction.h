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
    /* Combinations that return long double */                          \
    template<> struct Helper<long double,long double> { typedef long double type; }; \
    template<> struct Helper<long double,double> { typedef long double type; }; \
    template<> struct Helper<double,long double> { typedef long double type; }; \
    template<> struct Helper<long double,float> { typedef long double type; }; \
    template<> struct Helper<float,long double> { typedef long double type; }; \
    template<> struct Helper<long double,unsigned long> { typedef long double type; }; \
    template<> struct Helper<unsigned long,long double> { typedef long double type; }; \
    template<> struct Helper<long double,signed long> { typedef long double type; }; \
    template<> struct Helper<signed long,long double> { typedef long double type; }; \
    template<> struct Helper<long double,unsigned int> { typedef long double type; }; \
    template<> struct Helper<unsigned int,long double> { typedef long double type; }; \
    template<> struct Helper<long double,signed int> { typedef long double type; }; \
    template<> struct Helper<signed int,long double> { typedef long double type; }; \
    template<> struct Helper<long double,unsigned short> { typedef long double type; }; \
    template<> struct Helper<unsigned short,long double> { typedef long double type; }; \
    template<> struct Helper<long double,signed short> { typedef long double type; }; \
    template<> struct Helper<signed short,long double> { typedef long double type; }; \
    template<> struct Helper<long double,unsigned char> { typedef long double type; }; \
    template<> struct Helper<unsigned char,long double> { typedef long double type; }; \
    template<> struct Helper<long double,signed char> { typedef long double type; }; \
    template<> struct Helper<signed char,long double> { typedef long double type; }; \
    /* Combinations that return double */                               \
    template<> struct Helper<double,double> { typedef double type; };   \
    template<> struct Helper<double,float> { typedef double type; };    \
    template<> struct Helper<float,double> { typedef double type; };    \
    template<> struct Helper<double,unsigned long> { typedef double type; }; \
    template<> struct Helper<unsigned long,double> { typedef double type; }; \
    template<> struct Helper<double,signed long> { typedef double type; }; \
    template<> struct Helper<signed long,double> { typedef double type; }; \
    template<> struct Helper<double,unsigned int> { typedef double type; }; \
    template<> struct Helper<unsigned int,double> { typedef double type; }; \
    template<> struct Helper<double,signed int> { typedef double type; }; \
    template<> struct Helper<signed int,double> { typedef double type; }; \
    template<> struct Helper<double,unsigned short> { typedef double type; }; \
    template<> struct Helper<unsigned short,double> { typedef double type; }; \
    template<> struct Helper<double,signed short> { typedef double type; }; \
    template<> struct Helper<signed short,double> { typedef double type; }; \
    template<> struct Helper<double,unsigned char> { typedef double type; }; \
    template<> struct Helper<unsigned char,double> { typedef double type; }; \
    template<> struct Helper<double,signed char> { typedef double type; }; \
    template<> struct Helper<signed char,double> { typedef double type; }; \
    /* Combinations that return float */                                \
    template<> struct Helper<float,float> { typedef float type; };      \
    template<> struct Helper<unsigned long,float> { typedef float type; }; \
    template<> struct Helper<float,signed long> { typedef float type; }; \
    template<> struct Helper<signed long,float> { typedef float type; }; \
    template<> struct Helper<float,unsigned int> { typedef float type; }; \
    template<> struct Helper<unsigned int,float> { typedef float type; }; \
    template<> struct Helper<float,signed int> { typedef float type; }; \
    template<> struct Helper<signed int,float> { typedef float type; }; \
    template<> struct Helper<float,unsigned short> { typedef float type; }; \
    template<> struct Helper<unsigned short,float> { typedef float type; }; \
    template<> struct Helper<float,signed short> { typedef float type; }; \
    template<> struct Helper<signed short,float> { typedef float type; }; \
    template<> struct Helper<float,unsigned char> { typedef float type; }; \
    template<> struct Helper<unsigned char,float> { typedef float type; }; \
    template<> struct Helper<float,signed char> { typedef float type; }; \
    template<> struct Helper<signed char,float> { typedef float type; }; \
    /* Combinations that return unsigned long */                        \
    template<> struct Helper<unsigned long,unsigned long> { typedef unsigned long type; }; \
    template<> struct Helper<unsigned long,signed long> { typedef unsigned long type; }; \
    template<> struct Helper<signed long,unsigned long> { typedef unsigned long type; }; \
    template<> struct Helper<unsigned long,unsigned int> { typedef unsigned long type; }; \
    template<> struct Helper<unsigned int,unsigned long> { typedef unsigned long type; }; \
    template<> struct Helper<unsigned long,signed int> { typedef unsigned long type; }; \
    template<> struct Helper<signed int,unsigned long> { typedef unsigned long type; }; \
    template<> struct Helper<unsigned long,unsigned short> { typedef unsigned long type; }; \
    template<> struct Helper<unsigned short,unsigned long> { typedef unsigned long type; }; \
    template<> struct Helper<unsigned long,signed short> { typedef unsigned long type; }; \
    template<> struct Helper<signed short,unsigned long> { typedef unsigned long type; }; \
    template<> struct Helper<unsigned long,unsigned char> { typedef unsigned long type; }; \
    template<> struct Helper<unsigned char,unsigned long> { typedef unsigned long type; }; \
    template<> struct Helper<unsigned long,signed char> { typedef unsigned long type; }; \
    template<> struct Helper<signed char,unsigned long> { typedef unsigned long type; }; \
    /* Combinations that return signed long */                          \
    template<> struct Helper<signed long,signed long> { typedef signed long type; }; \
    template<> struct Helper<signed long,unsigned int> { typedef signed long type; }; \
    template<> struct Helper<unsigned int,signed long> { typedef signed long type; }; \
    template<> struct Helper<signed long,signed int> { typedef signed long type; }; \
    template<> struct Helper<signed int,signed long> { typedef signed long type; }; \
    template<> struct Helper<signed long,unsigned short> { typedef signed long type; }; \
    template<> struct Helper<unsigned short,signed long> { typedef signed long type; }; \
    template<> struct Helper<signed long,signed short> { typedef signed long type; }; \
    template<> struct Helper<signed short,signed long> { typedef signed long type; }; \
    template<> struct Helper<signed long,unsigned char> { typedef signed long type; }; \
    template<> struct Helper<unsigned char,signed long> { typedef signed long type; }; \
    template<> struct Helper<signed long,signed char> { typedef signed long type; }; \
    template<> struct Helper<signed char,signed long> { typedef signed long type; }; \
    /* Combinations that return unsigned int */                         \
    template<> struct Helper<unsigned int,unsigned int> { typedef unsigned int type; }; \
    template<> struct Helper<unsigned int,signed int> { typedef unsigned int type; }; \
    template<> struct Helper<signed int,unsigned int> { typedef unsigned int type; }; \
    template<> struct Helper<unsigned int,unsigned short> { typedef unsigned int type; }; \
    template<> struct Helper<unsigned short,unsigned int> { typedef unsigned int type; }; \
    template<> struct Helper<unsigned int,signed short> { typedef unsigned int type; }; \
    template<> struct Helper<signed short,unsigned int> { typedef unsigned int type; }; \
    template<> struct Helper<unsigned int,unsigned char> { typedef unsigned int type; }; \
    template<> struct Helper<unsigned char,unsigned int> { typedef unsigned int type; }; \
    template<> struct Helper<unsigned int,signed char> { typedef unsigned int type; }; \
    template<> struct Helper<signed char,unsigned int> { typedef unsigned int type; }; \
    /* Combinations that return signed int */                           \
    template<> struct Helper<signed int,signed int> { typedef signed int type; }; \
    template<> struct Helper<signed int,unsigned short> { typedef signed int type; }; \
    template<> struct Helper<unsigned short,signed int> { typedef signed int type; }; \
    template<> struct Helper<signed int,signed short> { typedef signed int type; }; \
    template<> struct Helper<signed short,signed int> { typedef signed int type; }; \
    template<> struct Helper<signed int,unsigned char> { typedef signed int type; }; \
    template<> struct Helper<unsigned char,signed int> { typedef signed int type; }; \
    template<> struct Helper<signed int,signed char> { typedef signed int type; }; \
    template<> struct Helper<signed char,signed int> { typedef signed int type; }; \
    /* Combinations that return unsigned short */                       \
    template<> struct Helper<unsigned short,unsigned short> { typedef unsigned short type; }; \
    template<> struct Helper<unsigned short,signed short> { typedef unsigned short type; }; \
    template<> struct Helper<signed short,unsigned short> { typedef unsigned short type; }; \
    template<> struct Helper<unsigned short,unsigned char> { typedef unsigned short type; }; \
    template<> struct Helper<unsigned char,unsigned short> { typedef unsigned short type; }; \
    template<> struct Helper<unsigned short,signed char> { typedef unsigned short type; }; \
    template<> struct Helper<signed char,unsigned short> { typedef unsigned short type; }; \
    /* Combinations that return signed short */                       \
    template<> struct Helper<signed short,signed short> { typedef signed short type; }; \
    template<> struct Helper<signed short,unsigned char> { typedef signed short type; }; \
    template<> struct Helper<unsigned char,signed short> { typedef signed short type; }; \
    template<> struct Helper<signed short,signed char> { typedef signed short type; }; \
    template<> struct Helper<signed char,signed short> { typedef signed short type; }; \
    /* Combinations that return unsigned char */                        \
    template<> struct Helper<unsigned char,unsigned char> { typedef unsigned char type; }; \
    template<> struct Helper<unsigned char,signed char> { typedef unsigned char type; }; \
    template<> struct Helper<signed char,unsigned char> { typedef unsigned char type; }; \
    /* Combinations that return signed char */                       \
    template<> struct Helper<signed char,signed char> { typedef signed char type; }; \
    /* Combinations with a user type return that type */                \
    template<class T> struct Helper<long double,T> { typedef T type; }; \
    template<class T> struct Helper<T,long double> { typedef T type; }; \
    template<class T> struct Helper<double,T> { typedef T type; };      \
    template<class T> struct Helper<T,double> { typedef T type; };      \
    template<class T> struct Helper<float,T> { typedef T type; };       \
    template<class T> struct Helper<T,float> { typedef T type; };       \
    template<class T> struct Helper<unsigned long,T> { typedef T type; }; \
    template<class T> struct Helper<T,unsigned long> { typedef T type; }; \
    template<class T> struct Helper<signed long,T> { typedef T type; }; \
    template<class T> struct Helper<T,signed long> { typedef T type; }; \
    template<class T> struct Helper<unsigned int,T> { typedef T type; }; \
    template<class T> struct Helper<T,unsigned int> { typedef T type; };  \
    template<class T> struct Helper<signed int,T> { typedef T type; }; \
    template<class T> struct Helper<T,signed int> { typedef T type; };  \
    template<class T> struct Helper<unsigned short,T> { typedef T type; }; \
    template<class T> struct Helper<T,unsigned short> { typedef T type; }; \
    template<class T> struct Helper<signed short,T> { typedef T type; }; \
    template<class T> struct Helper<T,signed short> { typedef T type; }; \
    template<class T> struct Helper<unsigned char,T> { typedef T type; }; \
    template<class T> struct Helper<T,unsigned char> { typedef T type; }; \
    template<class T> struct Helper<signed char,T> { typedef T type; }; \
    template<class T> struct Helper<T,signed char> { typedef T type; };

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
