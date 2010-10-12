// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file PixelMask.h
///
/// Defines the useful pixel utility type that can wrap any existing
/// pixel type and add mask semantics.  Any operations with an
/// "invalid" pixel returns an invalid pixel as a result.
///
#ifndef __VW_IMAGE_PIXELMASK_H__
#define __VW_IMAGE_PIXELMASK_H__

#include <ostream>

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/static_assert.hpp>

#include <vw/Image/PixelMath.h>
#include <vw/Image/PixelTypeInfo.h>

namespace vw {

  // *******************************************************************
  // The PixelMask wrapper pixel type.
  // *******************************************************************

  /// A generic wrapper for any of the above pixel types that adds an
  /// additional "valid" bit to the pixel.  Math operations that
  /// include invalide pixels will produce resulting pixels that are
  /// themselves invalid.
  template <class ChildT>
  struct PixelMask : public PixelMathBase< PixelMask<ChildT> >
  {
    typedef typename CompoundChannelType<ChildT>::type channel_type;

  private:
    ChildT m_child;
    channel_type m_valid;

  public:

    // Default constructor (zero value).  Pixel is not valid by
    // default.
    PixelMask() {
      m_child = ChildT();
      m_valid = ChannelRange<channel_type>::min();
    }

    /// Conversion from other PixelMask<> types.
    template <class T>
    PixelMask( PixelMask<T> other ) {
      // We let the child's built-in conversions do their work here.
      // This will fail if there is no conversion defined from
      // OtherT to ChildT.
      m_child = ChildT(other.child());

      if (other.valid())
        m_valid = ChannelRange<channel_type>::max();
      else
        m_valid = ChannelRange<channel_type>::min();
    }

    /// implicit construction from the raw channel value, or from the
    /// child pixel type value.  Values constructed in this manner
    /// are considered valid.
    template <class T>
    PixelMask( T const& pix) {
      m_child = ChildT(pix);
      m_valid = ChannelRange<channel_type>::max();
    }

    /// Constructs a pixel with the given channel values (use only when child has 2 channels)
    PixelMask( channel_type const& a0, channel_type const& a1 ) {
      m_child[0]=a0; m_child[1]=a1; m_valid = ChannelRange<channel_type>::max();
    }

    /// Constructs a pixel with the given channel values (use only when child has 3 channels)
    PixelMask( channel_type const& a0, channel_type const& a1, channel_type const& a2 ) {
      m_child[0]=a0; m_child[1]=a1; m_child[2]=a2; m_valid = ChannelRange<channel_type>::max();
    }

    /// Constructs a pixel with the given channel values (use only when child has 4 channels)
    PixelMask( channel_type const& a0, channel_type const& a1, channel_type const& a2, channel_type const& a3 ) {
      m_child[0]=a0; m_child[1]=a1; m_child[2]=a2; m_child[3]=a3; m_valid = ChannelRange<channel_type>::max();
    }

    /// Returns the value from the valid channel
    channel_type valid() const { return m_valid; }

    /// Invalidates this pixel, setting its valid bit to zero.
    void invalidate() { m_valid = ChannelRange<channel_type>::min(); }

    /// Invalidates this pixel, setting its valid bit to 1;
    void validate() { m_valid = ChannelRange<channel_type>::max(); }

    /// Switching valid to what it was not previously
    void toggle() {
      if ( m_valid )
        m_valid = ChannelRange<channel_type>::min();
      else
        m_valid = ChannelRange<channel_type>::max();
    }

    /// Returns the child pixel type
    ChildT & child() { return m_child; }
    ChildT const& child() const { return m_child; }

    /// Automatic down-cast to the raw channel value in numeric
    /// contexts.  This should only work for pixels that contain one
    /// data channel (plus the mask channel).  We add a
    /// BOOST_STATIC_ASSERT here to make sure that this is the case.
    operator channel_type() const {
       BOOST_STATIC_ASSERT(CompoundNumChannels<ChildT>::value == 1);
       return compound_select_channel<channel_type const&>(m_child,0);
    }

    /// Channel indexing operator.
    inline channel_type& operator[](size_t i) {
      if (i == CompoundNumChannels<ChildT>::value)
        return m_valid;
      else
        return compound_select_channel<channel_type&>(m_child,i);
     }
    /// Channel indexing operator (const overload).
    inline channel_type const& operator[](size_t i) const {
      if (i == CompoundNumChannels<ChildT>::value)
        return m_valid;
      else
        return compound_select_channel<channel_type const&>(m_child,i);
    }
    /// Channel indexing operator.
    inline channel_type& operator()(size_t i) {
      if (i == CompoundNumChannels<ChildT>::value)
        return valid;
      else
        return compound_select_channel<channel_type&>(m_child,i);
    }
    /// Channel indexing operator (const overload).
    inline channel_type const& operator()(size_t i) const {
      if (i == CompoundNumChannels<ChildT>::value)
        return valid;
      else
        return compound_select_channel<channel_type const&>(m_child,i);
    }
  };

  // *******************************************************************
  // Generic Mask Manipulation Methods
  // *******************************************************************

  // Overload for the pixel transparency traits class.
  template <class ChildT>
  bool is_transparent(PixelMask<ChildT> const& pixel) { return !pixel.valid(); }

  // Overload for the pixel transparency traits class.
  template <class ChildT>
  bool is_opaque(PixelMask<ChildT> const& pixel) { return pixel.valid(); }

  // A pixel validity checker.
  template <class PixelT>
  bool is_valid(PixelT const& /* pixel */) { return true; }

  // Overload for the pixel transparency traits class.
  template <class ChildT>
  bool is_valid(PixelMask<ChildT> const& pixel) { return pixel.valid(); }

  // Generic method for "validating" pixel (setting the mask bit).
  // This is a no-op by default, but it actually calls px.validate()
  // for PixelMask<> types.
  template <class PixelT>
  inline void validate(PixelT & /* pixel */ ) { return; }

  template <class ChildPixelT>
  inline void validate(PixelMask<ChildPixelT> &pixel) { pixel.validate(); }

  // Generic method for "invalidating" pixel (removing the mask bit).
  // This is a no-op by default, but it actually calls px.validate()
  // for PixelMask<> types.
  template <class PixelT>
  inline void invalidate(PixelT& /*pixel*/) { return; }

  template <class ChildPixelT>
  inline void invalidate(PixelMask<ChildPixelT> &pixel) { pixel.invalidate(); }

  // Generic method for "toggling" a pixel (not'ing the mask bit).
  // This is a no-op by default, but it actually calls px.validate()
  // for PixelMask<> types.
  template <class PixelT>
  inline void toggle(PixelT &/*pixel*/) { return; }

  template <class ChildPixelT>
  inline void toggle(PixelMask<ChildPixelT> &pixel) { pixel.toggle(); }

  // Remove mask
  template <class PixelT>
  inline PixelT& remove_mask(PixelT& pixel) { return pixel; }

  template <class PixelT>
  inline const PixelT& remove_mask(const PixelT& pixel) { return pixel; }

  template <class PixelT>
  inline PixelT& remove_mask(PixelMask<PixelT>& pixel) { return pixel.child(); }

  template <class PixelT>
  inline const PixelT& remove_mask(const PixelMask<PixelT>& pixel) { return pixel.child(); }

  /// Print a PixelMask to a debugging stream.
  template <class ChildT>
  std::ostream& operator<<( std::ostream& os, PixelMask<ChildT> const& pix ) {
    return os << "PixelMask( " << pix.child() << " : " << _numeric(pix.valid()) << " )";
  }

  // *******************************************************************
  // Type Traits
  // *******************************************************************

  template <class T>
  struct CompoundChannelType<PixelMask<T> > {
    typedef typename CompoundChannelType<T>::type type;
  };
  template <class T>
  struct CompoundNumChannels<PixelMask<T> > {
    static const size_t value = CompoundNumChannels<T>::value + 1;
  };
  template <class OldT, class NewChT>
  struct CompoundChannelCast<PixelMask<OldT>, NewChT> {
    typedef PixelMask<typename CompoundChannelCast<OldT,NewChT>::type> type;
  };
  template <class OldT, class NewChT>
  struct CompoundChannelCast<PixelMask<OldT>, const NewChT> {
    typedef const PixelMask<typename CompoundChannelCast<OldT,NewChT>::type> type;
  };
  template <class T> struct IsMasked : public boost::false_type::type {};
  template <class T> struct IsMasked<PixelMask<T> > : public boost::true_type::type {};


  // Computes the mean value of a compound PixelMask<> type.  Not
  // especially efficient.
  template <class T>
  typename boost::enable_if< IsScalarOrCompound<T>, double >::type
  inline mean_channel_value( PixelMask<T> const& arg ) {
    typedef typename CompoundChannelType<T>::type channel_type;
    if (arg.valid()) {
      size_t num_channels = CompoundNumChannels<T>::value;
      double accum = 0;
      for( size_t i=0; i<num_channels-1; ++i )
        accum += compound_select_channel<channel_type const&>( arg, i );
      return accum / num_channels;
    } else {
      return 0;
    }
  }

  // These are handy for determining what the masked type is for a
  // given pixel type.
  template <class PixelT>
  struct MaskedPixelType { typedef PixelMask<PixelT> type; };

  // Most importantly, we want to make sure that we only wrap each
  // pixel in one layer of PixelMask<>.
  template <class ChildT>
  struct MaskedPixelType<PixelMask<ChildT> > { typedef PixelMask<ChildT> type; };

  // These are handy for determining what the masked type is for a
  // given pixel type.
  template <class PixelT>
  struct UnmaskedPixelType { typedef PixelT type; };

  // Most importantly, we want to make sure that we only wrap each
  // pixel in one layer of PixelMask<>.
  template <class ChildT>
  struct UnmaskedPixelType<PixelMask<ChildT> > { typedef ChildT type; };

  // *******************************************************************
  // Binary elementwise compound type functor.
  // *******************************************************************

  template <class FuncT, class ChildPixel1T, class ChildPixel2T>
  class BinaryCompoundFunctor<FuncT, PixelMask<ChildPixel1T>, PixelMask<ChildPixel2T> > {
    FuncT func;

    // The general multi-channel case
    template <bool CompoundB, size_t ChannelsN, class ResultT, class Arg1T, class Arg2T>
    struct Helper {
      static inline ResultT construct( FuncT const& func, Arg1T const& arg1, Arg2T const& arg2 ) {
        ResultT result;
        for( size_t i=0; i<ChannelsN-1; ++i ) result[i] = func(arg1[i],arg2[i]);
        if (arg1.valid() && arg2.valid())
          result.validate();
        return result;
      }
    };

    // Specialization for one-channel + 1 "valid pixel" channel
    template <class ResultT, class Arg1T, class Arg2T>
    struct Helper<true,2,ResultT,Arg1T,Arg2T> {
      static inline ResultT construct( FuncT const& func, Arg1T const& arg1, Arg2T const& arg2 ) {
        ResultT result( func(arg1[0],arg2[0]) );
        if (!arg1.valid() || !arg2.valid())
          result.invalidate();
        return result;
      }
    };

    // Specialization for two-channel + 1 valid pixel channel types
    template <class ResultT, class Arg1T, class Arg2T>
    struct Helper<true,3,ResultT,Arg1T,Arg2T> {
      static inline ResultT construct( FuncT const& func, Arg1T const& arg1, Arg2T const& arg2 ) {
        ResultT result( func(arg1[0],arg2[0]), func(arg1[1],arg2[1]) );
        if (!arg1.valid() || !arg2.valid())
          result.invalidate();
        return result;
      }
    };

    // Specialization for three-channel + 1 valid pixel channel types
    template <class ResultT, class Arg1T, class Arg2T>
    struct Helper<true,4,ResultT,Arg1T,Arg2T> {
      static inline ResultT construct( FuncT const& func, Arg1T const& arg1, Arg2T const& arg2 ) {
        ResultT result( func(arg1[0],arg2[0]), func(arg1[1],arg2[1]), func(arg1[2],arg2[2]) );
        if (!arg1.valid() || !arg2.valid())
          result.invalidate();
        return result;
      }
    };

    // Specialization for four-channel + 1 valid pixel channel types
    template <class ResultT, class Arg1T, class Arg2T>
    struct Helper<true,5,ResultT,Arg1T,Arg2T> {
      static inline ResultT construct( FuncT const& func, Arg1T const& arg1, Arg2T const& arg2 ) {
        ResultT result( func(arg1[0],arg2[0]), func(arg1[1],arg2[1]), func(arg1[2],arg2[2]), func(arg1[3],arg2[3]) );
        if (!arg1.valid() || !arg2.valid())
          result.invalidate();
        return result;
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

  // *******************************************************************
  // Binary in-place elementwise compound type functor.
  // *******************************************************************

  template <class FuncT, class ChildPixel1T, class ChildPixel2T>
  class BinaryInPlaceCompoundFunctor<FuncT, PixelMask<ChildPixel1T>, PixelMask<ChildPixel2T> > {
    FuncT func;

    // The general multi-channel case
    template <bool CompoundB, size_t ChannelsN, class Arg1T, class Arg2T>
    struct Helper {
      static inline Arg1T& apply( FuncT const& func, Arg1T& arg1, Arg2T const& arg2 ) {
        for( size_t i=0; i<ChannelsN-1; ++i ) func(arg1[i],arg2[i]);
        if (!arg2.valid())
          arg1.invalidate();
        return arg1;
      }
    };

    // Specialization for one-channel types + 1 "valid pixel" channel
    template <class Arg1T, class Arg2T>
    struct Helper<true,2,Arg1T,Arg2T> {
      static inline Arg1T& apply( FuncT const& func, Arg1T& arg1, Arg2T const& arg2 ) {
        func(arg1[0],arg2[0]);
        if (!arg2.valid())
          arg1.invalidate();
        return arg1;
      }
    };

    // Specialization for two-channel types + 1 "valid pixel" channel
    template <class Arg1T, class Arg2T>
    struct Helper<true,3,Arg1T,Arg2T> {
      static inline Arg1T& apply( FuncT const& func, Arg1T& arg1, Arg2T const& arg2 ) {
        func(arg1[0],arg2[0]);
        func(arg1[1],arg2[1]);
        if (!arg2.valid())
          arg1.invalidate();
        return arg1;
      }
    };

    // Specialization for three-channel types + 1 "valid pixel" channel
    template <class Arg1T, class Arg2T>
    struct Helper<true,4,Arg1T,Arg2T> {
      static inline Arg1T& apply( FuncT const& func, Arg1T& arg1, Arg2T const& arg2 ) {
        func(arg1[0],arg2[0]);
        func(arg1[1],arg2[1]);
        func(arg1[2],arg2[2]);
        if (!arg2.valid())
          arg1.invalidate();
        return arg1;
      }
    };

    // Specialization for four-channel types + 1 "valid pixel" channel
    template <class Arg1T, class Arg2T>
    struct Helper<true,5,Arg1T,Arg2T> {
      static inline Arg1T& apply( FuncT const& func, Arg1T& arg1, Arg2T const& arg2 ) {
        func(arg1[0],arg2[0]);
        func(arg1[1],arg2[1]);
        func(arg1[2],arg2[2]);
        func(arg1[3],arg2[3]);
        if (!arg2.valid())
          arg1.invalidate();
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

  // *******************************************************************
  // Unary elementwise compound type functor.
  // *******************************************************************

  template <class FuncT, class ChildPixelT>
  class UnaryCompoundFunctor<FuncT, PixelMask<ChildPixelT> > {
    FuncT func;

    // The general multi-channel case
    template <bool CompoundB, size_t ChannelsN, class ResultT, class ArgT>
    struct Helper {
      static inline ResultT construct( FuncT const& func, ArgT const& arg ) {
        ResultT result;
        for( size_t i=0; i<ChannelsN-1; ++i ) result[i] = func(arg[i]);
        if (arg.valid())
          result.validate();
        return result;
      }
    };

    // Specialization for one-channel + 1 "valid pixel" channel
    template <class ResultT, class ArgT>
    struct Helper<true,2,ResultT,ArgT> {
      static inline ResultT construct( FuncT const& func, ArgT const& arg ) {
        ResultT result( func(arg[0]) );
        if (!arg.valid())
          result.invalidate();
        return result;
      }
    };

    // Specialization for two-channel + 1 valid pixel channel types
    template <class ResultT, class ArgT>
    struct Helper<true,3,ResultT,ArgT> {
      static inline ResultT construct( FuncT const& func, ArgT const& arg ) {
        ResultT result( func(arg[0]), func(arg[1]) );
        if (!arg.valid())
          result.invalidate();
        return result;
      }
    };

    // Specialization for three-channel + 1 valid pixel channel types
    template <class ResultT, class ArgT>
    struct Helper<true,4,ResultT,ArgT> {
      static inline ResultT construct( FuncT const& func, ArgT const& arg ) {
        ResultT result( func(arg[0]), func(arg[1]), func(arg[2]) );
        if (!arg.valid())
          result.invalidate();
        return result;
      }
    };

    // Specialization for four-channel + 1 valid pixel channel types
    template <class ResultT, class ArgT>
    struct Helper<true,5,ResultT,ArgT> {
      static inline ResultT construct( FuncT const& func, ArgT const& arg ) {
        ResultT result( func(arg[0]), func(arg[1]), func(arg[2]), func(arg[3]) );
        if (!arg.valid())
          result.invalidate();
        return result;
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

  // *******************************************************************
  // Unary in-place elementwise compound type functor.
  // *******************************************************************

  template <class FuncT, class ChildPixelT>
  class UnaryInPlaceCompoundFunctor<FuncT, PixelMask<ChildPixelT> > {
    FuncT func;
    typedef typename boost::add_reference<FuncT>::type func_ref;

    // The general multi-channel case
    template <bool CompoundB, size_t ChannelsN, class ArgT>
    struct Helper {
      static inline ArgT& apply( func_ref func, ArgT& arg ) {
        for( size_t i=0; i<ChannelsN-1; ++i ) func(arg[i]);
        return arg;
      }
    };

    // Specialization for one-channel types + 1 "valid pixel" channel
    template <class ArgT>
    struct Helper<true,2,ArgT> {
      static inline ArgT& apply( func_ref func, ArgT& arg ) {
        func(arg[0]);
        return arg;
      }
    };

    // Specialization for two-channel types + 1 "valid pixel" channel
    template <class ArgT>
    struct Helper<true,3,ArgT> {
      static inline ArgT& apply( func_ref func, ArgT& arg ) {
        func(arg[0]);
        func(arg[1]);
        return arg;
      }
    };

    // Specialization for three-channel types + 1 "valid pixel" channel
    template <class ArgT>
    struct Helper<true,4,ArgT> {
      static inline ArgT& apply( func_ref func, ArgT& arg ) {
        func(arg[0]);
        func(arg[1]);
        func(arg[2]);
        return arg;
      }
    };

    // Specialization for four-channel types + 1 "valid pixel" channel
    template <class ArgT>
    struct Helper<true,5,ArgT> {
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

  // ******************************************************************
  // Special Quotient Safe Functors for Pixel Mask types
  struct ArgArgMaskedSafeQuotientFunctor : BinaryReturnTemplateType<QuotientType> {
    template <class Arg1T, class Arg2T>
    inline typename QuotientType<PixelMask<Arg1T>, PixelMask<Arg2T> >::type
    operator()( PixelMask<Arg1T> const& arg1, PixelMask<Arg2T> const& arg2 ) const {
      if ( arg2.child() == Arg2T() ) {
        if ( is_valid(arg1) && is_valid(arg2) )
          return typename QuotientType<PixelMask<Arg1T>,PixelMask<Arg2T> >::type(0); // Valid
        else
          return typename QuotientType<PixelMask<Arg1T>,PixelMask<Arg2T> >::type();
      }
      else return ( arg1 / arg2 );
    }
  };

  template <class ValT>
  struct ValArgMaskedSafeQuotientFunctor : UnaryReturnBinaryTemplateBind1st<QuotientType,ValT> {
  private:
    const ValT m_val;
  public:
    ValArgMaskedSafeQuotientFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline typename QuotientType<ValT, PixelMask<ArgT> >::type
    operator()( PixelMask<ArgT> const& arg ) const {
      if ( arg.child()==ArgT() ) {
        if ( is_valid(arg) && is_valid(m_val) )
          return typename QuotientType<ValT,PixelMask<ArgT> >::type(0);
        else
          return typename QuotientType<ValT,PixelMask<ArgT> >::type();
      }
      else return (m_val / arg);
    }
  };

  struct ArgArgInPlaceMaskedSafeQuotientFunctor : BinaryReturn1stType {
    template <class Arg1T, class Arg2T>
    inline PixelMask<Arg1T>& operator()( PixelMask<Arg1T>& arg1,
                                         PixelMask<Arg2T> const& arg2 ) const {
      if ( arg2.child()==Arg2T() ) {
        if ( is_valid(arg1) && is_valid(arg2) )
          return arg1=PixelMask<Arg1T>(0);
        else
          return arg1=PixelMask<Arg1T>();
      } else
        return arg1=(PixelMask<Arg1T>)(arg1/arg2);
    }
  };

  // *******************************************************************
  // The Pixel Math specialization for math against PixelMask<Scalar>
  //
  // Yet this should specifically not do PixelMask<scalar> against
  // PixelMask<scalar> as that becomes an ambiguous call against the
  // normal math operators.
  // *******************************************************************
#define VW_PIXEL_MASK_MATH_BINARY_PS_FUNCTION(func,ftor)                                \
  template <class PixelT, class ScalarT>                                                \
    typename boost::enable_if< boost::mpl::and_<boost::mpl::and_<IsScalar<ScalarT>, IsMasked<PixelT> >, \
                               boost::mpl::not_<boost::is_same<PixelT,PixelMask<ScalarT> > > >, \
                               typename CompoundResult<ftor<ScalarT>,PixelT>::type >::type \
  inline func( PixelMathBase<PixelT> const& pixel, PixelMask<ScalarT> masked_scalar ) { \
    if (!masked_scalar.valid()) {                                                       \
      PixelT px = pixel.impl();                                                         \
      px.invalidate();                                                                  \
      return compound_apply(ftor<ScalarT>(masked_scalar.child()), px );                 \
    } else {                                                                            \
      return compound_apply(ftor<ScalarT>(masked_scalar.child()), pixel.impl() );       \
    }                                                                                   \
  }

#define VW_PIXEL_MASK_MATH_BINARY_SP_FUNCTION(func,ftor)                                \
  template <class PixelT, class ScalarT>                                                \
    typename boost::enable_if< boost::mpl::and_<boost::mpl::and_<IsScalar<ScalarT>, IsMasked<PixelT> >, \
                               boost::mpl::not_<boost::is_same<PixelT,PixelMask<ScalarT> > > >, \
                               typename CompoundResult<ftor<ScalarT>,PixelT>::type >::type \
  inline func( PixelMask<ScalarT> masked_scalar, PixelMathBase<PixelT> const& pixel ) { \
    if (!masked_scalar.valid()) {                                                       \
      PixelT px = pixel.impl();                                                         \
      px.invalidate();                                                                  \
      return compound_apply(ftor<ScalarT>(masked_scalar.child()), px );                 \
    } else {                                                                            \
      return compound_apply(ftor<ScalarT>(masked_scalar.child()), pixel.impl() );       \
    }                                                                                   \
  }

#define VW_PIXEL_MASK_MATH_BINARY_IS_FUNCTION(func,ftor)                                 \
  template <class PixelT, class ScalarT>                                                 \
    typename boost::enable_if< boost::mpl::and_<IsScalar<ScalarT>, IsMasked<PixelT> >, PixelT&>::type \
  inline func( PixelMathBase<PixelT>& pixel, PixelMask<ScalarT> masked_scalar ) {        \
    if (!masked_scalar.valid())                                                          \
      pixel.impl().invalidate();                                                         \
    return compound_apply_in_place(ftor<ScalarT>(masked_scalar.child()), pixel.impl() ); \
  }

  VW_PIXEL_MASK_MATH_BINARY_PS_FUNCTION(operator +, vw::ArgValSumFunctor)
  VW_PIXEL_MASK_MATH_BINARY_SP_FUNCTION(operator +, vw::ValArgSumFunctor)
  VW_PIXEL_MASK_MATH_BINARY_IS_FUNCTION(operator +=, vw::ArgValInPlaceSumFunctor)
  VW_PIXEL_MASK_MATH_BINARY_PS_FUNCTION(operator -, vw::ArgValDifferenceFunctor)
  VW_PIXEL_MASK_MATH_BINARY_SP_FUNCTION(operator -, vw::ValArgDifferenceFunctor)
  VW_PIXEL_MASK_MATH_BINARY_IS_FUNCTION(operator -=, vw::ArgValInPlaceDifferenceFunctor)
  VW_PIXEL_MASK_MATH_BINARY_PS_FUNCTION(operator *, vw::ArgValProductFunctor)
  VW_PIXEL_MASK_MATH_BINARY_SP_FUNCTION(operator *, vw::ValArgProductFunctor)
  VW_PIXEL_MASK_MATH_BINARY_IS_FUNCTION(operator *=, vw::ArgValInPlaceProductFunctor)
  VW_PIXEL_MASK_MATH_BINARY_PS_FUNCTION(operator /, vw::ArgValQuotientFunctor)
  VW_PIXEL_MASK_MATH_BINARY_SP_FUNCTION(operator /, vw::ValArgQuotientFunctor)
  VW_PIXEL_MASK_MATH_BINARY_IS_FUNCTION(operator /=, vw::ArgValInPlaceQuotientFunctor)



} // namespace vw

#endif // __VW_IMAGE_PIXELMASK_H__
