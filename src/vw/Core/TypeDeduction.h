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

#include <boost/config.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/greater.hpp>
#include <boost/mpl/identity.hpp>

#include <vw/Core/CompoundTypes.h>

namespace vw {
  // If you have an error message that leads you here, you tried to apply a
  // functor to two user-defined (non-arithmetic) types, and we didn't know
  // what type to return.
  template <class T>
  struct TypeDeductionError;

  namespace core { namespace detail {
    struct unknown_type_tag : public boost::mpl::int_<10000> {};

    // This should be considered an implementation detail
#define _VW_INTERNAL_TYPE_DEDUCTION(Type, Number)            \
    BOOST_STATIC_ASSERT(boost::is_arithmetic<Type>::value);  \
    template <>                                              \
    struct TypeDeductionIndex<Type> {                        \
      typedef boost::mpl::int_<Number> type;                 \
      BOOST_STATIC_CONSTANT(unsigned, value = type::value);  \
    };                                                       \

    template <class T>
    struct TypeDeductionIndex {
      typedef unknown_type_tag type;
      BOOST_STATIC_CONSTANT(unsigned, value = type::value);
    };

    _VW_INTERNAL_TYPE_DEDUCTION(              char,  100);
    _VW_INTERNAL_TYPE_DEDUCTION(       signed char,  200);
    _VW_INTERNAL_TYPE_DEDUCTION(     unsigned char,  300);
    _VW_INTERNAL_TYPE_DEDUCTION(      signed short,  400);
    _VW_INTERNAL_TYPE_DEDUCTION(    unsigned short,  500);
    _VW_INTERNAL_TYPE_DEDUCTION(        signed int,  600);
    _VW_INTERNAL_TYPE_DEDUCTION(      unsigned int,  700);
    _VW_INTERNAL_TYPE_DEDUCTION(       signed long,  800);
    _VW_INTERNAL_TYPE_DEDUCTION(     unsigned long,  900);
#if defined(BOOST_HAS_LONG_LONG)
    _VW_INTERNAL_TYPE_DEDUCTION(  signed long long, 1000);
    _VW_INTERNAL_TYPE_DEDUCTION(unsigned long long, 1100);
#endif
#if defined(BOOST_HAS_MS_INT64)
    _VW_INTERNAL_TYPE_DEDUCTION(          __int64;, 1200);
    _VW_INTERNAL_TYPE_DEDUCTION( unsigned __int64;, 1300);
#endif
    _VW_INTERNAL_TYPE_DEDUCTION(             float, 1400);
    _VW_INTERNAL_TYPE_DEDUCTION(            double, 1500);
    _VW_INTERNAL_TYPE_DEDUCTION(       long double, 1600);
  }}

  template <class T1, class T2>
  struct TypeDeductionHelper {
    typedef typename core::detail::TypeDeductionIndex<T1>::type I1;
    typedef typename core::detail::TypeDeductionIndex<T2>::type I2;

    // If an error message leads you here, you tried to apply a functor to two
    // user-defined (non-arithmetic) types, and we didn't know what type to
    // return. (We catch the case where the types are actually the same in the
    // next specialization.)
    BOOST_STATIC_ASSERT((!boost::mpl::equal_to<I1, I2>::value));

    // If an error message leads you here, you have an arithmetic type that we
    // haven't identified. Please report this as a bug.
    //
    // (Implementation note: _VW_INTERNAL_TYPE_DEDUCTION contains an
    // is_arithmetic check which makes this logic correct)
    BOOST_STATIC_ASSERT(!(boost::mpl::and_<boost::is_arithmetic<T1>, boost::mpl::equal_to<I1, core::detail::unknown_type_tag> >::value));
    BOOST_STATIC_ASSERT(!(boost::mpl::and_<boost::is_arithmetic<T2>, boost::mpl::equal_to<I2, core::detail::unknown_type_tag> >::value));

    typedef typename
      boost::mpl::eval_if< boost::mpl::greater<I1, I2>
      , boost::mpl::identity<T1>
      , boost::mpl::identity<T2>
      >::type type;
  };

  // Both types are the same? Just use that.
  template <class ArgT>
  struct TypeDeductionHelper<ArgT, ArgT> {
    typedef ArgT type;
  };

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

  template <class Arg1T, class Arg2T, bool SameWrapper>
  struct PromoteTypeSpecializationHelper2;

  template <class Arg1T, class Arg2T>
  struct PromoteTypeSpecializationHelper2<Arg1T, Arg2T, true> {
    typedef typename CompoundChannelCast<Arg1T, typename PromoteTypeSpecialization<typename CompoundChannelType<Arg1T>::type, typename CompoundChannelType<Arg2T>::type>::type>::type type;
  };

  template <class Arg1T, class Arg2T>
  struct PromoteTypeSpecializationHelper2<Arg1T, Arg2T, false> {
    typedef typename TypeDeductionHelper<Arg1T,Arg2T>::type type;
  };

  // Recursive behavior when both Arg1T and Arg2T are compound types.
  template <class Arg1T, class Arg2T>
  struct PromoteTypeSpecializationHelper<Arg1T,Arg2T,true,true> {
    typedef typename PromoteTypeSpecializationHelper2<Arg1T, Arg2T,
      boost::is_same<Arg1T, typename CompoundChannelCast<Arg2T,typename CompoundChannelType<Arg1T>::type>::type>::value>::type type;
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

#endif
