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
#include <boost/mpl/greater.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/sizeof.hpp>
#include <boost/mpl/identity.hpp>

#include <vw/Core/CompoundTypes.h>

namespace vw {

  // If you have an error message that leads you here, you tried to apply a
  // functor to two user-defined (non-arithmetic) types, and we didn't know
  // what type to return. If you know, you can put a template specialization
  // for it, like:
  //    VW_NEW_TYPE_DEDUCTION1(user1_t, result_t);
  //    VW_NEW_TYPE_DEDUCTION2(user1_t, user2_t, result_t);
  template <class T>
  struct TypeDeductionError;

  #define VW_NEW_TYPE_DEDUCTION1(T1, Result)                             \
    namespace vw { namespace core { namespace detail {                   \
    template <class T2>                                                  \
    struct TypeDeductionHelperImpl<T1, T2, false, false, false, false> { \
      typedef Result type;                                               \
    };                                                                   \
    template <class T2>                                                  \
    struct TypeDeductionHelperImpl<T2, T1, false, false, false, false> { \
      typedef Result type;                                               \
    };                                                                   \
    }}}

  #define VW_NEW_TYPE_DEDUCTION2(T1, T2, Result)                         \
    namespace vw { namespace core { namespace detail {                   \
    template <>                                                          \
    struct TypeDeductionHelperImpl<T1, T2, false, false, false, false> { \
      typedef Result type;                                               \
    };                                                                   \
    template <>                                                          \
    struct TypeDeductionHelperImpl<T2, T1, false, false, false, false> { \
      typedef Result type;                                               \
    };                                                                   \
    }}}

  namespace core { namespace detail {
    template <class T1, class T2, bool T1Float, bool T2Float, bool T1Arith, bool T2Arith>
    struct TypeDeductionHelperImpl;

    // Both types are different user-defined types [since they're not arithmetic].
    // No solution. Error.
    template <class T1, class T2, bool T1Float, bool T2Float>
    struct TypeDeductionHelperImpl<T1, T2, T1Float, T2Float, false, false> {
      typedef TypeDeductionError<TypeDeductionHelperImpl> type;
    };
    // Only one type is arithmetic, pick the user-defined type.
    template <class T1, class T2, bool T1Float, bool T2Float>
    struct TypeDeductionHelperImpl<T1, T2, T1Float, T2Float, true, false> {
      typedef T2 type;
    };
    // Only one type is arithmetic, pick the user-defined type.
    template <class T1, class T2, bool T1Float, bool T2Float>
    struct TypeDeductionHelperImpl<T1, T2, T1Float, T2Float, false, true> {
      typedef T1 type;
    };
    // both arithmetic, only one is floating. pick the float.
    template <class T1, class T2>
    struct TypeDeductionHelperImpl<T1, T2, false, true, true, true> {
      typedef T2 type;
    };
    // both arithmetic, only one is floating. pick the float.
    template <class T1, class T2>
    struct TypeDeductionHelperImpl<T1, T2, true, false, true, true> {
      typedef T1 type;
    };
    // both arithmetic, both float. pick the bigger float.
    template <class T1, class T2>
    struct TypeDeductionHelperImpl<T1, T2, true, true, true, true> {
      typedef typename boost::mpl::if_<
        typename boost::mpl::greater<typename boost::mpl::sizeof_<T1>, typename boost::mpl::sizeof_<T2> >::type, T1, T2
      >::type type;
    };
    // both arithmetic, neither float. We use a modification of the C++ arithmetic expression rules here.
    // if sizeof(T1) == sizeof(T2):
    //   if is_unsigned(T1) or is_unsigned(T2): return make_unsigned(T1)
    //   else: return make_signed(T1)
    // else:
    //  return typeof(bigger_type)
    template <class T1, class T2>
    struct TypeDeductionHelperImpl<T1, T2, false, false, true, true> {
      BOOST_STATIC_ASSERT(boost::mpl::sizeof_<T1>::value > 0);
      BOOST_STATIC_ASSERT(boost::mpl::sizeof_<T2>::value > 0);

      typedef typename
        boost::mpl::eval_if< boost::mpl::equal_to<boost::mpl::sizeof_<T1>, boost::mpl::sizeof_<T2> >
          , boost::mpl::eval_if< boost::mpl::or_<boost::is_unsigned<T1>, boost::is_unsigned<T2> >
            , boost::make_unsigned<T2>
            , boost::make_signed<T2> >
          , boost::mpl::eval_if< boost::mpl::greater<boost::mpl::sizeof_<T1>, boost::mpl::sizeof_<T2> >
            , boost::mpl::identity<T1>
            , boost::mpl::identity<T2> >
      >::type type;
    };
  }}

  template <class T1, class T2>
  struct TypeDeductionHelper;

  // Both types are the same? Just use that.
  template <class ArgT>
  struct TypeDeductionHelper<ArgT, ArgT> {
    typedef ArgT type;
  };

  template <class T1, class T2>
  struct TypeDeductionHelper {
    typedef typename core::detail::TypeDeductionHelperImpl<T1, T2
            , boost::is_floating_point<T1>::value
            , boost::is_floating_point<T2>::value
            , boost::is_arithmetic<T1>::value
            , boost::is_arithmetic<T2>::value
            >::type type;
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
