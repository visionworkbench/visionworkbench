// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file CompoundTypes.h
///
/// Types and traits for compound (i.e. multi-channel) types.
///
#ifndef __VW_CORE_COMPOUNDTYPES_H__
#define __VW_CORE_COMPOUNDTYPES_H__

#include <boost/utility/result_of.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/and.hpp>

#include <vw/Core/FundamentalTypes.h>

namespace vw {

  // *******************************************************************
  // Compount type traits classes.
  // *******************************************************************

  // Default compound type traits templates.  Compound types are mainly used
  // as pixel types, but the type traits machinery is defined here in more
  // general terms to avoid undesirable dependencies.
  template <class T> struct CompoundChannelType { typedef T type; };
  template <class T> struct CompoundNumChannels { static const size_t value = 1; };
  template <class T, class ChannelT> struct CompoundChannelCast { typedef ChannelT type; };
  template <class T> struct IsCompound
    : public boost::mpl::not_< boost::is_same< typename CompoundChannelType<T>::type, T > >::type {};
  template <class T> struct IsScalarOrCompound
    : public boost::mpl::or_< typename IsScalar<T>::type, typename IsCompound<T>::type >::type {};
  template <class T1, class T2> struct CompoundIsCompatible
    : public boost::is_same< typename CompoundChannelCast<T1, typename CompoundChannelType<T2>::type>::type, T2 >::type {};
  template <class T> struct CompoundAccumulatorType {
    typedef typename CompoundChannelCast<T, typename AccumulatorType<typename CompoundChannelType<T>::type>::type>::type type;
  };

  // Default specializations of the compound type traits for const types.
  template <class T> struct CompoundChannelType<const T> : public CompoundChannelType<T> {};
  template <class T> struct CompoundNumChannels<const T> : public CompoundNumChannels<T> {};
  template <class T, class ChannelT> struct CompoundChannelCast<const T, ChannelT> : public CompoundChannelCast<T, ChannelT> {};
  template <class T> struct IsCompound<const T> : public IsCompound<T> {};
  template <class T1, class T2> struct CompoundIsCompatible<T1, const T2> : public CompoundIsCompatible<T1,T2> {};

  // These functions take a required ResultT template parameter so the caller
  // can specify whether to return by value or by reference, and likewise
  // whether to accept the first argument by value or reference.  This is all
  // rather annoying, and results in four totally incomprehensible versions of
  // a function to do one simple thing.  Is there a better way?  This would all
  // be somewhat cleaner if we could do the function enabling in the return type
  // rather than the second argument, but that still breaks on some not-so-old
  // compilers.
  template <class ResultT, class PixelT>
  inline ResultT compound_select_channel( PixelT& pixel, typename boost::enable_if<typename boost::mpl::and_<typename boost::mpl::not_< IsCompound<PixelT> >::type, typename boost::is_reference<ResultT>::type>, size_t>::type /*channel*/ ) {
    return pixel;
  }

  template <class ResultT, class PixelT>
  inline ResultT compound_select_channel( PixelT const& pixel, typename boost::enable_if<typename boost::mpl::and_<typename boost::mpl::not_< IsCompound<PixelT> >::type, typename boost::mpl::not_<typename boost::is_reference<ResultT>::type>::type >, size_t>::type /*channel*/ ) {
    return pixel;
  }

  template <class ResultT, class PixelT>
  inline ResultT compound_select_channel( PixelT& pixel, typename boost::enable_if<typename boost::mpl::and_<IsCompound<PixelT>, typename boost::is_reference<ResultT>::type>, size_t>::type channel ) {
    return pixel[channel];
  }

  template <class ResultT, class PixelT>
  inline ResultT compound_select_channel( PixelT const& pixel, typename boost::enable_if<typename boost::mpl::and_<IsCompound<PixelT>, typename boost::mpl::not_<typename boost::is_reference<ResultT>::type>::type >, size_t>::type channel ) {
    return pixel[channel];
  }

  // Computes the mean value of a compound type.  Not especially efficient.
  // I'm not really sure this is the best place for this to live....
  template <class T>
  typename boost::enable_if< IsScalarOrCompound<T>, double >::type
  inline mean_channel_value( T const& arg ) {
    typedef typename CompoundChannelType<T>::type channel_type;
    size_t num_channels = CompoundNumChannels<T>::value;
    double accum = 0;
    for( size_t i=0; i<num_channels; ++i )
      accum += compound_select_channel<channel_type const&>( arg, i );
    return accum / num_channels;
  }

  // *******************************************************************
  // Binary elementwise compound type functor.
  // *******************************************************************

  template <class FuncT, class FuncArg1T, class FuncArg2T>
  class BinaryCompoundFunctor {
    FuncT func;

    // The general multi-channel case
    template <bool CompoundB, size_t ChannelsN, class ResultT, class Arg1T, class Arg2T>
    struct Helper {
      static inline ResultT construct( FuncT const& func, Arg1T const& arg1, Arg2T const& arg2 ) {
        ResultT result;
        for( size_t i=0; i<ChannelsN; ++i ) result[i] = func(arg1[i],arg2[i]);
        return result;
      }
    };

    // Specialization for non-compound types
    template <class ResultT, class Arg1T, class Arg2T>
    struct Helper<false,1,ResultT,Arg1T,Arg2T> {
      static inline ResultT construct( FuncT const& func, Arg1T const& arg1, Arg2T const& arg2 ) {
        return func(arg1,arg2);
      }
    };

    // Specialization for one-channel types
    template <class ResultT, class Arg1T, class Arg2T>
    struct Helper<true,1,ResultT,Arg1T,Arg2T> {
      static inline ResultT construct( FuncT const& func, Arg1T const& arg1, Arg2T const& arg2 ) {
        return ResultT( func(arg1[0],arg2[0]) );
      }
    };

    // Specialization for two-channel types
    template <class ResultT, class Arg1T, class Arg2T>
    struct Helper<true,2,ResultT,Arg1T,Arg2T> {
      static inline ResultT construct( FuncT const& func, Arg1T const& arg1, Arg2T const& arg2 ) {
        return ResultT( func(arg1[0],arg2[0]), func(arg1[1],arg2[1]) );
      }
    };

    // Specialization for three-channel types
    template <class ResultT, class Arg1T, class Arg2T>
    struct Helper<true,3,ResultT,Arg1T,Arg2T> {
      static inline ResultT construct( FuncT const& func, Arg1T const& arg1, Arg2T const& arg2 ) {
        return ResultT( func(arg1[0],arg2[0]), func(arg1[1],arg2[1]), func(arg1[2],arg2[2]) );
      }
    };

    // Specialization for four-channel types
    template <class ResultT, class Arg1T, class Arg2T>
    struct Helper<true,4,ResultT,Arg1T,Arg2T> {
      static inline ResultT construct( FuncT const& func, Arg1T const& arg1, Arg2T const& arg2 ) {
        return ResultT( func(arg1[0],arg2[0]), func(arg1[1],arg2[1]), func(arg1[2],arg2[2]), func(arg1[3],arg2[3]) );
      }
    };

  public:
    BinaryCompoundFunctor() : func() {}
    BinaryCompoundFunctor( FuncT const& func ) : func(func) {}

    template <class ArgsT> struct result {};

    template <class F, class Arg1T, class Arg2T>
    struct result<F(Arg1T,Arg2T)> {
      typedef typename CompoundChannelType<Arg1T>::type arg1_type;
      typedef typename CompoundChannelType<Arg2T>::type arg2_type;
      typedef typename boost::result_of<FuncT(arg1_type,arg2_type)>::type result_type;
      typedef typename CompoundChannelCast<Arg1T,result_type>::type type;
    };

    template <class Arg1T, class Arg2T>
    typename result<BinaryCompoundFunctor(Arg1T,Arg2T)>::type
    inline operator()( Arg1T const& arg1, Arg2T const& arg2 ) const {
      typedef typename result<BinaryCompoundFunctor(Arg1T,Arg2T)>::type result_type;
      return Helper<IsCompound<result_type>::value,CompoundNumChannels<result_type>::value,result_type,Arg1T,Arg2T>::construct(func,arg1,arg2);
    }
  };

  template <class FuncT, class Arg1T, class Arg2T=void>
  struct CompoundResult {
    typedef typename boost::result_of<BinaryCompoundFunctor<FuncT, Arg1T, Arg2T>(Arg1T,Arg2T)>::type type;
  };

  template <class FuncT, class Arg1T, class Arg2T>
  typename CompoundResult<FuncT,Arg1T,Arg2T>::type
  inline compound_apply( FuncT const& func, Arg1T const& arg1, Arg2T const& arg2 ) {
    return BinaryCompoundFunctor<FuncT,Arg1T,Arg2T>(func)(arg1,arg2);
  }


  // *******************************************************************
  // Binary in-place elementwise compound type functor.
  // *******************************************************************

  template <class FuncT, class FuncArg1T, class FuncArg2T>
  class BinaryInPlaceCompoundFunctor {
    FuncT func;

    // The general multi-channel case
    template <bool CompoundB, size_t ChannelsN, class Arg1T, class Arg2T>
    struct Helper {
      static inline Arg1T& apply( FuncT const& func, Arg1T& arg1, Arg2T const& arg2 ) {
        for( size_t i=0; i<ChannelsN; ++i ) func(arg1[i],arg2[i]);
        return arg1;
      }
    };

    // Specialization for non-compound types
    template <class Arg1T, class Arg2T>
    struct Helper<false,1,Arg1T,Arg2T> {
      static inline Arg1T& apply( FuncT const& func, Arg1T& arg1, Arg2T const& arg2 ) {
        func(arg1,arg2);
        return arg1;
      }
    };

    // Specialization for one-channel types
    template <class Arg1T, class Arg2T>
    struct Helper<true,1,Arg1T,Arg2T> {
      static inline Arg1T& apply( FuncT const& func, Arg1T& arg1, Arg2T const& arg2 ) {
        func(arg1[0],arg2[0]);
        return arg1;
      }
    };

    // Specialization for two-channel types
    template <class Arg1T, class Arg2T>
    struct Helper<true,2,Arg1T,Arg2T> {
      static inline Arg1T& apply( FuncT const& func, Arg1T& arg1, Arg2T const& arg2 ) {
        func(arg1[0],arg2[0]);
        func(arg1[1],arg2[1]);
        return arg1;
      }
    };

    // Specialization for three-channel types
    template <class Arg1T, class Arg2T>
    struct Helper<true,3,Arg1T,Arg2T> {
      static inline Arg1T& apply( FuncT const& func, Arg1T& arg1, Arg2T const& arg2 ) {
        func(arg1[0],arg2[0]);
        func(arg1[1],arg2[1]);
        func(arg1[2],arg2[2]);
        return arg1;
      }
    };

    // Specialization for four-channel types
    template <class Arg1T, class Arg2T>
    struct Helper<true,4,Arg1T,Arg2T> {
      static inline Arg1T& apply( FuncT const& func, Arg1T& arg1, Arg2T const& arg2 ) {
        func(arg1[0],arg2[0]);
        func(arg1[1],arg2[1]);
        func(arg1[2],arg2[2]);
        func(arg1[3],arg2[3]);
        return arg1;
      }
    };

  public:
    BinaryInPlaceCompoundFunctor() : func() {}
    BinaryInPlaceCompoundFunctor( FuncT const& func ) : func(func) {}

    template <class ArgsT> struct result {};

    template <class F, class Arg1T, class Arg2T>
    struct result<F(Arg1T,Arg2T)> {
      typedef Arg1T& type;
    };

    template <class Arg1T, class Arg2T>
    typename result<BinaryInPlaceCompoundFunctor(Arg1T,Arg2T)>::type
    inline operator()( Arg1T& arg1, Arg2T const& arg2 ) const {
      return Helper<IsCompound<Arg1T>::value,CompoundNumChannels<Arg1T>::value,Arg1T,Arg2T>::apply(func,arg1,arg2);
    }
  };

  template <class FuncT, class Arg1T, class Arg2T>
  inline Arg1T& compound_apply_in_place( FuncT const& func, Arg1T& arg1, Arg2T const& arg2 ) {
    return BinaryInPlaceCompoundFunctor<FuncT,Arg1T,Arg2T>(func)(arg1,arg2);
  }


  // *******************************************************************
  // Unary elementwise compound type functor.
  // *******************************************************************

  template <class FuncT, class FuncArgT>
  class UnaryCompoundFunctor {
    FuncT func;

    // The general multi-channel case
    template <bool CompoundB, size_t ChannelsN, class ResultT, class ArgT>
    struct Helper {
      static inline ResultT construct( FuncT const& func, ArgT const& arg ) {
        ResultT result;
        for( size_t i=0; i<ChannelsN; ++i ) result[i] = func(arg[i]);
        return result;
      }
    };

    // Specialization for non-compound types
    template <class ResultT, class ArgT>
    struct Helper<false,1,ResultT,ArgT> {
      static inline ResultT construct( FuncT const& func, ArgT const& arg ) {
        return func(arg);
      }
    };

    // Specialization for single-channel types
    template <class ResultT, class ArgT>
    struct Helper<true,1,ResultT,ArgT> {
      static inline ResultT construct( FuncT const& func, ArgT const& arg ) {
        return ResultT( func(arg[0]) );
      }
    };

    // Specialization for two-channel types
    template <class ResultT, class ArgT>
    struct Helper<true,2,ResultT,ArgT> {
      static inline ResultT construct( FuncT const& func, ArgT const& arg ) {
        return ResultT( func(arg[0]), func(arg[1]) );
      }
    };

    // Specialization for three-channel types
    template <class ResultT, class ArgT>
    struct Helper<true,3,ResultT,ArgT> {
      static inline ResultT construct( FuncT const& func, ArgT const& arg ) {
        return ResultT( func(arg[0]), func(arg[1]), func(arg[2]) );
      }
    };

    // Specialization for four-channel types
    template <class ResultT, class ArgT>
    struct Helper<true,4,ResultT,ArgT> {
      static inline ResultT construct( FuncT const& func, ArgT const& arg ) {
        return ResultT( func(arg[0]), func(arg[1]), func(arg[2]), func(arg[3]) );
      }
    };

  public:
    UnaryCompoundFunctor() : func() {}
    UnaryCompoundFunctor( FuncT const& func ) : func(func) {}

    template <class ArgsT> struct result {};

    template <class F, class ArgT>
    struct result<F(ArgT)> {
      typedef typename CompoundChannelType<ArgT>::type arg_type;
      typedef typename boost::result_of<FuncT(arg_type)>::type result_type;
      typedef typename CompoundChannelCast<ArgT,result_type>::type type;
    };

    template <class ArgT>
    typename result<UnaryCompoundFunctor(ArgT)>::type
    inline operator()( ArgT const& arg ) const {
      typedef typename result<UnaryCompoundFunctor(ArgT)>::type result_type;
      return Helper<IsCompound<result_type>::value,CompoundNumChannels<result_type>::value,result_type,ArgT>::construct(func,arg);
    }
  };

  template <class FuncT, class ArgT>
  struct CompoundResult<FuncT,ArgT,void> {
    typedef typename boost::result_of<UnaryCompoundFunctor<FuncT,ArgT>(ArgT)>::type type;
  };

  template <class FuncT, class ArgT>
  typename CompoundResult<FuncT,ArgT>::type
  inline compound_apply( FuncT const& func, ArgT const& arg ) {
    return UnaryCompoundFunctor<FuncT,ArgT>(func)(arg);
  }


  // *******************************************************************
  // Unary in-place elementwise compound type functor.
  // *******************************************************************

  template <class FuncT, class FuncArgT>
  class UnaryInPlaceCompoundFunctor {
    FuncT func;
    typedef typename boost::add_reference<FuncT>::type func_ref;

    // The general multi-channel case
    template <bool CompoundB, size_t ChannelsN, class ArgT>
    struct Helper {
      static inline ArgT& apply( func_ref func, ArgT& arg ) {
        for( size_t i=0; i<ChannelsN; ++i ) func(arg[i]);
        return arg;
      }
    };

    // Specialization for non-compound types
    template <class ArgT>
    struct Helper<false,1,ArgT> {
      static inline ArgT& apply( func_ref func, ArgT& arg ) {
        func(arg);
        return arg;
      }
    };

    // Specialization for single-channel types
    template <class ArgT>
    struct Helper<true,1,ArgT> {
      static inline ArgT& apply( func_ref func, ArgT& arg ) {
        func(arg[0]);
        return arg;
      }
    };

    // Specialization for two-channel types
    template <class ArgT>
    struct Helper<true,2,ArgT> {
      static inline ArgT& apply( func_ref func, ArgT& arg ) {
        func(arg[0]);
        func(arg[1]);
        return arg;
      }
    };

    // Specialization for three-channel types
    template <class ArgT>
    struct Helper<true,3,ArgT> {
      static inline ArgT& apply( func_ref func, ArgT& arg ) {
        func(arg[0]);
        func(arg[1]);
        func(arg[2]);
        return arg;
      }
    };

    // Specialization for four-channel types
    template <class ArgT>
    struct Helper<true,4,ArgT> {
      static inline ArgT& apply( func_ref func, ArgT& arg ) {
        func(arg[0]);
        func(arg[1]);
        func(arg[2]);
        func(arg[3]);
        return arg;
      }
    };

  public:
    UnaryInPlaceCompoundFunctor() : func() {}
    UnaryInPlaceCompoundFunctor( func_ref func ) : func(func) {}

    template <class ArgsT> struct result {};

    /// FIXME: This seems not to respect the constness of ArgT?  Weird?
    template <class F, class ArgT>
    struct result<F(ArgT)> {
      typedef ArgT& type;
    };

    template <class ArgT>
    inline ArgT& operator()( ArgT& arg ) const {
      return Helper<IsCompound<ArgT>::value,CompoundNumChannels<ArgT>::value,ArgT>::apply(func,arg);
    }

    template <class ArgT>
    inline const ArgT& operator()( const ArgT& arg ) const {
      return Helper<IsCompound<ArgT>::value,CompoundNumChannels<ArgT>::value,const ArgT>::apply(func,arg);
    }
  };

  template <class FuncT, class ArgT>
  inline ArgT& compound_apply_in_place( FuncT& func, ArgT& arg ) {
    return UnaryInPlaceCompoundFunctor<FuncT&,ArgT>(func)(arg);
  }

  template <class FuncT, class ArgT>
  inline const ArgT& compound_apply_in_place( FuncT& func, ArgT const& arg ) {
    return UnaryInPlaceCompoundFunctor<FuncT&,ArgT>(func)(arg);
  }

  template <class FuncT, class ArgT>
  inline ArgT& compound_apply_in_place( FuncT const& func, ArgT& arg ) {
    return UnaryInPlaceCompoundFunctor<FuncT const&,ArgT>(func)(arg);
  }

  template <class FuncT, class ArgT>
  inline const ArgT& compound_apply_in_place( FuncT const& func, ArgT const& arg ) {
    return UnaryInPlaceCompoundFunctor<FuncT const&,ArgT>(func)(arg);
  }

} // namespace vw

#endif // __VW_CORE_COMPOUNDTYPES_H__
