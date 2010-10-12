// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file PixelTypeInfo.h
///
/// Base classes, support classes, and other infrastructure
/// to support the pixel types.
///
#ifndef __VW_IMAGE_PIXELTYPEINFO_H__
#define __VW_IMAGE_PIXELTYPEINFO_H__

#include <complex>
#include <boost/integer_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/if.hpp>

#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/CompoundTypes.h>
#include <vw/Core/Functors.h>
#include <vw/Math/Functions.h>

namespace vw {

  // Forward declaration, used to disable certain pixel functions when
  // the arguments are actually images.
  template <class ImageT> struct ImageViewBase;

  // *******************************************************************
  // Basic pixel type manipulation logic
  //
  // Here we provide basic pixel type manipulation classes, most of
  // which are just thin wrappers provided for backwards
  // compatability.  We also provide tools to convert between real
  // and complex pixel types.
  // *******************************************************************

  template <class T> struct PixelChannelType : CompoundChannelType<T> {};
  template <class T> struct PixelNumChannels : CompoundNumChannels<T> {};
  template <class T, class ChannelT> struct PixelChannelCast : CompoundChannelCast<T,ChannelT> {};
  template <class T> struct PixelHasAlpha : false_type {};
  template <class T> struct PixelWithAlpha {};
  template <class T> struct PixelWithoutAlpha { typedef T type; };

  template <class PixelT>
  struct PixelMakeComplex {
    typedef typename PixelChannelCast<PixelT,typename MakeComplex<typename PixelChannelType<PixelT>::type>::type>::type type;
  };

  template <class PixelT>
  struct PixelMakeReal {
    typedef typename PixelChannelCast<PixelT,typename MakeReal<typename PixelChannelType<PixelT>::type>::type>::type type;
  };


  // *******************************************************************
  // Pixel channel standard range computation logic
  //
  // In some functions it is necessary to assume that the values in a
  // pixel span some meaningful range.  For example, when compositing
  // with alpha you have to know what value of alpha corresponds to
  // full opacity.  For floating-point types the range is generally
  // from zero to one.  These classes allow you to compute the range
  // for an arbitrary pixel type.
  // *******************************************************************

  // Channel range helper for integer types
  template <class T, bool>
  struct ChannelRangeHelper {
    static inline T max() {
      return boost::integer_traits<T>::max();
    }
    static inline T min() {
      return T();
    }
  };

  // Channel range helper for floating-point types
  template <class T>
  struct ChannelRangeHelper<T,false> {
    static inline T max() {
      return (T)(1.0);
    }
    static inline T min() {
      return (T)(0.0);
    }
  };

  // Channel range helper for complex types
  template <class T> struct ChannelRangeHelper<std::complex<T>,false> : public ChannelRangeHelper<T,boost::integer_traits<T>::is_integer> {};

  /// A channel range computation class.  This class is templatized on
  /// the pixel type and provides to static functions, max() and min(),
  /// that return the top and bottom of the standard range of the
  /// underlying channel type.  That is, max() corresponds to 1.0 in
  /// the floating-point world and min() corresponds to 0.0.  The value
  /// of min() is generally zero, *not* some negative number.  If you
  /// really just want to know the range of values you can store in a
  /// given type, use std::numeric_limits instead.
  template <class T> struct ChannelRange : public ChannelRangeHelper<typename CompoundChannelType<T>::type, boost::integer_traits<typename CompoundChannelType<T>::type>::is_integer> {};


  // *******************************************************************
  // Basic pixel alpha channel logic
  //
  // Here we provide basic pixel type manipulation classes for
  // determining whether a pixel is transparent (i.e. its alpha is
  // zero).  For pixel types without alpha, the pixel is always
  // considered to be non-transparent.
  // *******************************************************************
  template <class PixelT>
  inline typename boost::enable_if< typename IsScalarOrCompound<PixelT>::type, typename boost::enable_if<typename PixelHasAlpha<PixelT>::type, bool>::type>::type
  is_transparent(PixelT const& pixel) { return !(pixel.a()); }

  template <class PixelT>
  inline typename boost::enable_if< typename IsScalarOrCompound<PixelT>::type, typename boost::enable_if<typename PixelHasAlpha<PixelT>::type, bool>::type>::type
  is_opaque(PixelT const& pixel) { return pixel.a() == ChannelRange<PixelT>::max(); }

  template <class PixelT>
  inline typename boost::enable_if< typename IsScalarOrCompound<PixelT>::type, typename boost::disable_if<typename PixelHasAlpha<PixelT>::type, bool>::type>::type
  is_transparent(PixelT const& /*pixel*/) { return false; }

  template <class PixelT>
  inline typename boost::enable_if< typename IsScalarOrCompound<PixelT>::type, typename boost::disable_if<typename PixelHasAlpha<PixelT>::type, bool>::type>::type
  is_opaque(PixelT const& /*pixel*/) { return true; }

  template <class PixelT>
  inline typename boost::enable_if< typename IsScalarOrCompound<PixelT>::type, typename PixelChannelType<PixelT>::type >::type
  alpha_channel( PixelT& /*pix*/ ) {
    return ChannelRange<PixelT>::max();
  }

  template <class PixelT>
  inline typename boost::enable_if< typename IsScalarOrCompound<PixelT>::type, PixelT& >::type
  non_alpha_channels( PixelT& pix ) {
    return pix;
  }

  template <class PixelT>
  inline typename boost::enable_if< typename IsScalarOrCompound<PixelT>::type, const PixelT& >::type
  non_alpha_channels( const PixelT& pix ) {
    return pix;
  }


  // *******************************************************************
  // Pixel channel casting and rescaling logic.
  //
  // Here we defined channel_cast and channel_cast_rescale to operate
  // on pixels.  We've forward-declared ImageViewBase so that we can
  // disable these functions for images, allowing other whole-image
  // overloads to work in that case.  The simple channel_cast could be
  // put into Core/CompoundTypes.h instead, but it's only used in the
  // context of pixel casting and mirrors channel_cast_rescale nicely,
  // so we leave it here.
  //
  // FIXME The _clamp version clamps to the min/max integer values of
  // the destination.  This function does not work at all for floating-
  // point destination types, and the user might expect it to clamp to
  // the ChannelRange instead.  Probably the _clamp behavior should be
  // promoted to the default behavior for integer destination types.
  // Anyone who is relying on channel_cast's overvlow behavior probably
  // deserves to lose anyway.  At the moment the only thing using this
  // code is BicubicInterpolation.
  // *******************************************************************

  template <class DestT>
  class ChannelCastFunctor : public ReturnFixedType<DestT> {
  public:
    template <class SourceT>
    inline DestT operator()( SourceT source ) const {
      return (DestT)source;
    }
  };

  template <class ChannelT, class PixelT>
  typename boost::enable_if< typename IsScalarOrCompound<PixelT>::type, typename CompoundChannelCast<PixelT, ChannelT>::type >::type
  inline channel_cast( PixelT pixel ) {
    return compound_apply( ChannelCastFunctor<ChannelT>(), pixel );
  }

  template <class DestT>
  class ChannelCastRescaleFunctor : public ReturnFixedType<DestT> {
  public:
    template <class SourceT>
    inline DestT operator()( SourceT source ) const {
      // Clamping semantics are more reasonable for float->int rescaling.
      if( boost::is_floating_point<SourceT>::value && ! boost::is_floating_point<DestT>::value) {
        if( source > ChannelRange<SourceT>::max() ) source = ChannelRange<SourceT>::max();
        else if( source < ChannelRange<SourceT>::min() ) source = ChannelRange<SourceT>::min();
      }
      return (DestT)(source*(((double)(ChannelRange<DestT>::max()))/(ChannelRange<SourceT>::max())));
    }
    inline DestT operator()( DestT source ) const {
      return source;
    }
  };

  template <class ChannelT, class PixelT>
  typename boost::enable_if< typename IsScalarOrCompound<PixelT>::type, typename CompoundChannelCast<PixelT, ChannelT>::type >::type
  inline channel_cast_rescale( PixelT pixel ) {
    return compound_apply( ChannelCastRescaleFunctor<ChannelT>(), pixel );
  }

  template <class DestT>
  class ChannelCastClampFunctor : public ReturnFixedType<DestT> {
  public:
    template <class SourceT>
    inline DestT operator()( SourceT source ) const {
      if( source > boost::integer_traits<DestT>::max() ) return boost::integer_traits<DestT>::max();
      else if( source < boost::integer_traits<DestT>::min() ) return boost::integer_traits<DestT>::min();
      else return DestT( source );
    }
    inline DestT operator()( DestT source ) const {
      return source;
    }
  };

  template <class ChannelT, class PixelT>
  typename boost::enable_if< typename IsScalarOrCompound<PixelT>::type, typename CompoundChannelCast<PixelT, ChannelT>::type >::type
  inline channel_cast_clamp( PixelT pixel ) {
    return compound_apply( ChannelCastClampFunctor<ChannelT>(), pixel );
  }

  template <class ChannelT, class PixelT>
  typename boost::enable_if< typename IsScalarOrCompound<PixelT>::type, typename CompoundChannelCast<PixelT, ChannelT>::type >::type
  inline channel_cast_clamp_if_int( PixelT pixel ) {
    typedef typename boost::is_floating_point<ChannelT>::type is_float_type;
    typedef typename boost::mpl::if_<is_float_type, ChannelCastFunctor<ChannelT>, ChannelCastClampFunctor<ChannelT> >::type functor_type;
    return compound_apply( functor_type(), pixel );
  }

  template <class DestT>
  class ChannelCastRoundFunctor : public ReturnFixedType<DestT> {
  public:
    template <class SourceT>
    inline DestT operator()( SourceT source ) const {
      return DestT( math::impl::_round( source ) );
    }
    inline DestT operator()( DestT source ) const {
      return source;
    }
  };

  template <class ChannelT, class PixelT>
  typename boost::enable_if< typename IsScalarOrCompound<PixelT>::type, typename CompoundChannelCast<PixelT, ChannelT>::type >::type
  inline channel_cast_round( PixelT pixel ) {
    return compound_apply( ChannelCastRoundFunctor<ChannelT>(), pixel );
  }

  template <class ChannelT, class PixelT>
  typename boost::enable_if< typename IsScalarOrCompound<PixelT>::type, typename CompoundChannelCast<PixelT, ChannelT>::type >::type
  inline channel_cast_round_if_int( PixelT pixel ) {
    typedef typename boost::is_floating_point<ChannelT>::type is_float_type;
    typedef typename boost::mpl::if_<is_float_type, ChannelCastFunctor<ChannelT>, ChannelCastRoundFunctor<ChannelT> >::type functor_type;
    return compound_apply( functor_type(), pixel );
  }

  template <class DestT>
  class ChannelCastRoundClampFunctor : public ReturnFixedType<DestT> {
  public:
    template <class SourceT>
    inline DestT operator()( SourceT source ) const {
      if( source > boost::integer_traits<DestT>::max() ) return boost::integer_traits<DestT>::max();
      else if( source < boost::integer_traits<DestT>::min() ) return boost::integer_traits<DestT>::min();
      else return DestT( math::impl::_round( source ) );
    }
    inline DestT operator()( DestT source ) const {
      return source;
    }
  };

  template <class ChannelT, class PixelT>
  typename boost::enable_if< typename IsScalarOrCompound<PixelT>::type, typename CompoundChannelCast<PixelT, ChannelT>::type >::type
  inline channel_cast_round_and_clamp( PixelT pixel ) {
    return compound_apply( ChannelCastRoundClampFunctor<ChannelT>(), pixel );
  }

  template <class ChannelT, class PixelT>
  typename boost::enable_if< typename IsScalarOrCompound<PixelT>::type, typename CompoundChannelCast<PixelT, ChannelT>::type >::type
  inline channel_cast_round_and_clamp_if_int( PixelT pixel ) {
    typedef typename boost::is_floating_point<ChannelT>::type is_float_type;
    typedef typename boost::mpl::if_<is_float_type, ChannelCastFunctor<ChannelT>, ChannelCastRoundClampFunctor<ChannelT> >::type functor_type;
    return compound_apply( functor_type(), pixel );
  }

  // *******************************************************************
  // Pixel casting and rescaling logic.
  //
  // We implement pixel casting indirectly through this pixel_cast<>
  // template so that users can easily control the behavior of casts
  // between user-defined and standard pixel types.  When the source
  // and destination channel types differ, the channel conversion is
  // handled separately from the pixel format conversion.  The pixel
  // conversion itself takes place in whichever channel space has
  // greater precision, and occurs via a call to the corresponding
  // pixel_cast function.  Thus, if the user wishes to overload the
  // conversion behavior between two pixel formats, they need only
  // overload it for the case when the channel types are the same.
  // The rescaling variant, pixel_cast_rescale<>, performs a scale
  // conversion when converting the channel type.
  // *******************************************************************

  // Default behavior, invoked when the source and destination channel
  // type are the same.
  template <bool SameN, bool SrcN, bool RescaleN>
  struct PixelCastHelper {
    template <class DestT, class SrcT>
    static inline DestT convert( SrcT src ) {
      return DestT( src );
    }
  };

  // Non-rescaling pixel cast free function
  template <class DestT, class SrcT>
  typename boost::enable_if< typename IsScalarOrCompound<SrcT>::type, DestT >::type
  inline pixel_cast( SrcT src ) {
    typedef typename CompoundChannelType<SrcT>::type src_ch;
    typedef typename CompoundChannelType<DestT>::type dest_ch;
    typedef typename boost::is_same<src_ch,dest_ch>::type is_same;
    typedef PixelCastHelper< is_same::value, (sizeof(src_ch)>sizeof(dest_ch)), false > helper;
    return helper::template convert<DestT>( src );
  }

  // Rescaling pixel cast free function
  template <class DestT, class SrcT>
  typename boost::enable_if< typename IsScalarOrCompound<SrcT>::type, DestT >::type
  inline pixel_cast_rescale( SrcT src ) {
    typedef typename CompoundChannelType<SrcT>::type src_ch;
    typedef typename CompoundChannelType<DestT>::type dest_ch;
    typedef typename boost::is_same<src_ch,dest_ch>::type is_same;
    typedef PixelCastHelper< is_same::value, (sizeof(src_ch)>sizeof(dest_ch)), true > helper;
    return helper::template convert<DestT>( src );
  }

  // Non-rescaling behavior when SrcT has less or same precision
  template <>
  struct PixelCastHelper<false, false, false> {
    template <class DestT, class SrcT>
    static inline DestT convert( SrcT src ) {
      typedef typename CompoundChannelType<DestT>::type dest_ch;
      return pixel_cast<DestT>( channel_cast<dest_ch>( src ) );
    }
  };

  // Non-rescaling behavior when SrcT has greater precision
  template <>
  struct PixelCastHelper<false, true, false> {
    template <class DestT, class SrcT>
    static inline DestT convert( SrcT src ) {
      typedef typename CompoundChannelType<SrcT>::type src_ch;
      typedef typename CompoundChannelType<DestT>::type dest_ch;
      typedef typename CompoundChannelCast<DestT,src_ch>::type dest_px;
      return channel_cast<dest_ch>( pixel_cast<dest_px>( src ) );
    }
  };

  // Rescaling behavior when SrcT has less or same precision
  template <>
  struct PixelCastHelper<false, false, true> {
    template <class DestT, class SrcT>
    static inline DestT convert( SrcT src ) {
      typedef typename CompoundChannelType<DestT>::type dest_ch;
      return pixel_cast<DestT>( channel_cast_rescale<dest_ch>( src ) );
    }
  };

  // Rescaling behavior when SrcT has greater precision
  template <>
  struct PixelCastHelper<false, true, true> {
    template <class DestT, class SrcT>
    static inline DestT convert( SrcT src ) {
      typedef typename CompoundChannelType<SrcT>::type src_ch;
      typedef typename CompoundChannelType<DestT>::type dest_ch;
      typedef typename CompoundChannelCast<DestT,src_ch>::type dest_px;
      return channel_cast_rescale<dest_ch>( pixel_cast<dest_px>( src ) );
    }
  };


  // *******************************************************************
  // A pixel type convenience macro and forward declrations.
  // *******************************************************************

  /// This macro provides the appropriate specializations of
  /// the compound type traits classes for a new pixel type
  /// with a fixed number of channels (the common case).
#define VW_DECLARE_PIXEL_TYPE(PIXELT,NCHANNELS)              \
  template <class ChannelT>                                  \
  struct CompoundChannelType<PIXELT<ChannelT> > {            \
    typedef ChannelT type;                                   \
  };                                                         \
  template <class ChannelT>                                  \
  struct CompoundNumChannels<PIXELT<ChannelT> > {            \
    static const size_t value = NCHANNELS;                   \
  };                                                         \
  template <class OldChT, class NewChT>                      \
  struct CompoundChannelCast<PIXELT<OldChT>, NewChT> {       \
    typedef PIXELT<NewChT> type;                             \
  };                                                         \
  template <class OldChT, class NewChT>                      \
  struct CompoundChannelCast<PIXELT<OldChT>, const NewChT> { \
    typedef const PIXELT<NewChT> type;                       \
  }                                                          \

  /// This macro provides the appropriate specializations of
  /// the compound type traits classes for a new pixel type
  /// with a variable number of channels.
#define VW_DECLARE_PIXEL_TYPE_NCHANNELS(PIXELT)              \
  template <class ChannelT, size_t SizeN>                    \
  struct CompoundChannelType<PIXELT<ChannelT,SizeN> > {      \
    typedef ChannelT type;                                   \
  };                                                         \
  template <class ChannelT, size_t SizeN>                    \
  struct CompoundNumChannels<PIXELT<ChannelT,SizeN> > {      \
    static const size_t value = SizeN;                       \
  };                                                         \
  template <class OldChT, class NewChT, size_t SizeN>        \
  struct CompoundChannelCast<PIXELT<OldChT,SizeN>, NewChT> { \
    typedef PIXELT<NewChT,SizeN> type;                       \
  };                                                         \
  template <class OldChT, class NewChT, size_t SizeN>        \
  struct CompoundChannelCast<PIXELT<OldChT,SizeN>, const NewChT> { \
    typedef const PIXELT<NewChT,SizeN> type;                 \
  }

#define VW_ALPHA_PIXEL_TRAITS(opaque, alpha) \
  template <class T> struct PixelHasAlpha<alpha<T> > : true_type {};                  \
  template <class T> struct PixelWithAlpha<opaque<T> >    { typedef alpha<T> type; }; \
  template <class T> struct PixelWithAlpha<alpha<T> >     { typedef alpha<T> type; }; \
  template <class T> struct PixelWithoutAlpha<alpha<T> >  { typedef opaque<T> type; } \


  // Forward pixel type declarations for complex pixel types
  template <class ChannelT> class PixelGray;
  template <class ChannelT> class PixelGrayA;
  template <class ChannelT> class PixelRGB;
  template <class ChannelT> class PixelRGBA;
  template <class ChannelT> class PixelHSV;
  template <class ChannelT> class PixelXYZ;
  template <class ChannelT> class PixelLuv;
  template <class ChannelT> class PixelLab;

  // Forward pixel type declarations for the pixel mask wrapper
  template <class ChildT> struct PixelMask;

  // *******************************************************************
  // Run-time pixel type manipulation routines.
  // *******************************************************************

  enum PixelFormatEnum {
    VW_PIXEL_UNKNOWN = 0,
    VW_PIXEL_SCALAR = 1,
    VW_PIXEL_GRAY = 2,
    VW_PIXEL_GRAYA = 3,
    VW_PIXEL_RGB = 4,
    VW_PIXEL_RGBA = 5,
    VW_PIXEL_HSV = 6,
    VW_PIXEL_XYZ = 7,
    VW_PIXEL_LUV = 8,
    VW_PIXEL_LAB = 9,
    VW_PIXEL_UNKNOWN_MASKED = 40,
    VW_PIXEL_SCALAR_MASKED = 41,
    VW_PIXEL_GRAY_MASKED = 42,
    VW_PIXEL_GRAYA_MASKED = 43,
    VW_PIXEL_RGB_MASKED = 44,
    VW_PIXEL_RGBA_MASKED = 45,
    VW_PIXEL_HSV_MASKED = 46,
    VW_PIXEL_XYZ_MASKED = 47,
    VW_PIXEL_LUV_MASKED = 48,
    VW_PIXEL_LAB_MASKED = 49,
    VW_PIXEL_GENERIC_1_CHANNEL = 90,
    VW_PIXEL_GENERIC_2_CHANNEL = 91,
    VW_PIXEL_GENERIC_3_CHANNEL = 92,
    VW_PIXEL_GENERIC_4_CHANNEL = 93,
    VW_PIXEL_GENERIC_5_CHANNEL = 94,
    VW_PIXEL_GENERIC_6_CHANNEL = 95,
    VW_PIXEL_USER = 100
  };

  enum ChannelTypeEnum {
    VW_CHANNEL_UNKNOWN = 0,
    VW_CHANNEL_INT8 = 1,
    VW_CHANNEL_UINT8 = 2,
    VW_CHANNEL_INT16 = 3,
    VW_CHANNEL_UINT16 = 4,
    VW_CHANNEL_INT32 = 5,
    VW_CHANNEL_UINT32 = 6,
    VW_CHANNEL_INT64 = 7,
    VW_CHANNEL_UINT64 = 8,
    VW_CHANNEL_FLOAT16 = 9,
    VW_CHANNEL_FLOAT32 = 10,
    VW_CHANNEL_FLOAT64 = 11,
    VW_CHANNEL_BOOL = 12,
    VW_CHANNEL_CHAR = 13, // A tag to identify e.g. string medadata
    VW_CHANNEL_GENERIC_1_BYTE = 90,
    VW_CHANNEL_GENERIC_2_BYTE = 91,
    VW_CHANNEL_GENERIC_4_BYTE = 92,
    VW_CHANNEL_GENERIC_8_BYTE = 93,
    VW_CHANNEL_USER = 100
  };

  // PixelFormatID<>
  template <class PixelT> struct PixelFormatID { static const PixelFormatEnum value = VW_PIXEL_UNKNOWN; };
  template<> struct PixelFormatID<vw::int8>    { static const PixelFormatEnum value = VW_PIXEL_SCALAR; };
  template<> struct PixelFormatID<vw::uint8>   { static const PixelFormatEnum value = VW_PIXEL_SCALAR; };
  template<> struct PixelFormatID<vw::int16>   { static const PixelFormatEnum value = VW_PIXEL_SCALAR; };
  template<> struct PixelFormatID<vw::uint16>  { static const PixelFormatEnum value = VW_PIXEL_SCALAR; };
  template<> struct PixelFormatID<vw::int32>   { static const PixelFormatEnum value = VW_PIXEL_SCALAR; };
  template<> struct PixelFormatID<vw::uint32>  { static const PixelFormatEnum value = VW_PIXEL_SCALAR; };
  template<> struct PixelFormatID<vw::int64>   { static const PixelFormatEnum value = VW_PIXEL_SCALAR; };
  template<> struct PixelFormatID<vw::uint64>  { static const PixelFormatEnum value = VW_PIXEL_SCALAR; };
  template<> struct PixelFormatID<vw::float32> { static const PixelFormatEnum value = VW_PIXEL_SCALAR; };
  template<> struct PixelFormatID<vw::float64> { static const PixelFormatEnum value = VW_PIXEL_SCALAR; };
  template<> struct PixelFormatID<bool>        { static const PixelFormatEnum value = VW_PIXEL_SCALAR; };
  template <class ChT> struct PixelFormatID<PixelGray<ChT> >  { static const PixelFormatEnum value = VW_PIXEL_GRAY; };
  template <class ChT> struct PixelFormatID<PixelGrayA<ChT> > { static const PixelFormatEnum value = VW_PIXEL_GRAYA; };
  template <class ChT> struct PixelFormatID<PixelRGB<ChT> >   { static const PixelFormatEnum value = VW_PIXEL_RGB; };
  template <class ChT> struct PixelFormatID<PixelRGBA<ChT> >  { static const PixelFormatEnum value = VW_PIXEL_RGBA; };
  template <class ChT> struct PixelFormatID<PixelHSV<ChT> >   { static const PixelFormatEnum value = VW_PIXEL_HSV; };
  template <class ChT> struct PixelFormatID<PixelXYZ<ChT> >   { static const PixelFormatEnum value = VW_PIXEL_XYZ; };
  template <class ChT> struct PixelFormatID<PixelLuv<ChT> >   { static const PixelFormatEnum value = VW_PIXEL_LUV; };
  template <class ChT> struct PixelFormatID<PixelLab<ChT> >   { static const PixelFormatEnum value = VW_PIXEL_LAB; };

  // PixelFormatID<> specialized for masked pixel types
  template<> struct PixelFormatID<PixelMask<vw::int8> >    { static const PixelFormatEnum value = VW_PIXEL_SCALAR_MASKED; };
  template<> struct PixelFormatID<PixelMask<vw::uint8> >   { static const PixelFormatEnum value = VW_PIXEL_SCALAR_MASKED; };
  template<> struct PixelFormatID<PixelMask<vw::int16> >   { static const PixelFormatEnum value = VW_PIXEL_SCALAR_MASKED; };
  template<> struct PixelFormatID<PixelMask<vw::uint16> >  { static const PixelFormatEnum value = VW_PIXEL_SCALAR_MASKED; };
  template<> struct PixelFormatID<PixelMask<vw::int32> >   { static const PixelFormatEnum value = VW_PIXEL_SCALAR_MASKED; };
  template<> struct PixelFormatID<PixelMask<vw::uint32> >  { static const PixelFormatEnum value = VW_PIXEL_SCALAR_MASKED; };
  template<> struct PixelFormatID<PixelMask<vw::int64> >   { static const PixelFormatEnum value = VW_PIXEL_SCALAR_MASKED; };
  template<> struct PixelFormatID<PixelMask<vw::uint64> >  { static const PixelFormatEnum value = VW_PIXEL_SCALAR_MASKED; };
  template<> struct PixelFormatID<PixelMask<vw::float32> > { static const PixelFormatEnum value = VW_PIXEL_SCALAR_MASKED; };
  template<> struct PixelFormatID<PixelMask<vw::float64> > { static const PixelFormatEnum value = VW_PIXEL_SCALAR_MASKED; };
  template<> struct PixelFormatID<PixelMask<bool> >        { static const PixelFormatEnum value = VW_PIXEL_SCALAR_MASKED; };
  template <class ChT> struct PixelFormatID<PixelMask<PixelGray<ChT> > >  { static const PixelFormatEnum value = VW_PIXEL_GRAY_MASKED; };
  template <class ChT> struct PixelFormatID<PixelMask<PixelGrayA<ChT> > > { static const PixelFormatEnum value = VW_PIXEL_GRAYA_MASKED; };
  template <class ChT> struct PixelFormatID<PixelMask<PixelRGB<ChT> > >   { static const PixelFormatEnum value = VW_PIXEL_RGB_MASKED; };
  template <class ChT> struct PixelFormatID<PixelMask<PixelRGBA<ChT> > >  { static const PixelFormatEnum value = VW_PIXEL_RGBA_MASKED; };
  template <class ChT> struct PixelFormatID<PixelMask<PixelHSV<ChT> > >   { static const PixelFormatEnum value = VW_PIXEL_HSV_MASKED; };
  template <class ChT> struct PixelFormatID<PixelMask<PixelXYZ<ChT> > >   { static const PixelFormatEnum value = VW_PIXEL_XYZ_MASKED; };
  template <class ChT> struct PixelFormatID<PixelMask<PixelLuv<ChT> > >   { static const PixelFormatEnum value = VW_PIXEL_LUV_MASKED; };
  template <class ChT> struct PixelFormatID<PixelMask<PixelLab<ChT> > >   { static const PixelFormatEnum value = VW_PIXEL_LAB_MASKED; };

  // ChannelTypeID<>
  template <class ChannelT> struct ChannelTypeID { static const ChannelTypeEnum value = VW_CHANNEL_UNKNOWN; };
  template<> struct ChannelTypeID<vw::int8>      { static const ChannelTypeEnum value = VW_CHANNEL_INT8; };
  template<> struct ChannelTypeID<vw::uint8>     { static const ChannelTypeEnum value = VW_CHANNEL_UINT8; };
  template<> struct ChannelTypeID<vw::int16>     { static const ChannelTypeEnum value = VW_CHANNEL_INT16; };
  template<> struct ChannelTypeID<vw::uint16>    { static const ChannelTypeEnum value = VW_CHANNEL_UINT16; };
  template<> struct ChannelTypeID<vw::int32>     { static const ChannelTypeEnum value = VW_CHANNEL_INT32; };
  template<> struct ChannelTypeID<vw::uint32>    { static const ChannelTypeEnum value = VW_CHANNEL_UINT32; };
  template<> struct ChannelTypeID<vw::int64>     { static const ChannelTypeEnum value = VW_CHANNEL_INT64; };
  template<> struct ChannelTypeID<vw::uint64>    { static const ChannelTypeEnum value = VW_CHANNEL_UINT64; };
  template<> struct ChannelTypeID<vw::float32>   { static const ChannelTypeEnum value = VW_CHANNEL_FLOAT32; };
  template<> struct ChannelTypeID<vw::float64>   { static const ChannelTypeEnum value = VW_CHANNEL_FLOAT64; };
  template<> struct ChannelTypeID<bool>          { static const ChannelTypeEnum value = VW_CHANNEL_BOOL; };

  int32 channel_size( ChannelTypeEnum type );
  const char *channel_type_name( ChannelTypeEnum type );
  int32 num_channels( PixelFormatEnum format );
  const char *pixel_format_name( PixelFormatEnum format );
  ChannelTypeEnum channel_name_to_enum( const std::string& name );

} // namespace vw

#endif // __VW_IMAGE_PIXELTYPEINFO_H__
