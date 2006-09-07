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

/// \file PixelTypeInfo.h
/// 
/// Base classes, support classes, and other infrastructure 
/// to support the pixel types.
///
#ifndef __VW_IMAGE__PIXEL_TYPE_INFO_H__
#define __VW_IMAGE__PIXEL_TYPE_INFO_H__

#include <complex>
#include <boost/integer_traits.hpp>
#include <boost/utility/enable_if.hpp>

#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/CompoundTypes.h>
#include <vw/Core/Functors.h>
#include <vw/Image/ImageViewBase.h> // For IsImageView<>

namespace vw {

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
  template <class T> struct ChannelRange : public ChannelRangeHelper<typename CompoundChannelType<T>::type, boost::integer_traits<T>::is_integer> {};


  // *******************************************************************
  // Pixel channel casting and rescaling logic.
  //
  // Here we defined channel_cast and channel_cast_rescale to operate
  // on pixels.  We've included ImageViewBase.h so that we can disable
  // these functions for images, allowing other whole-image overloads
  // to work in that case.  The simple channel_cast could be put into
  // Core/CompoundTypes.h instead, but it's only used in the context
  // of pixel casting and mirrors channel_cast_rescale nicely, so we
  // leave it here.
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
  typename boost::disable_if< IsImageView<PixelT>, typename CompoundChannelCast<PixelT, ChannelT>::type >::type
  inline channel_cast( PixelT pixel ) {
    return compound_apply( ChannelCastFunctor<ChannelT>(), pixel );
  }

  template <class DestT>
  class ChannelCastRescaleFunctor : public ReturnFixedType<DestT> {
  public:
    template <class SourceT>
    inline DestT operator()( SourceT source ) const {
      return (DestT)(source*(((double)(ChannelRange<DestT>::max()))/(ChannelRange<SourceT>::max())));
    }
    inline DestT operator()( DestT source ) const {
      return source;
    }
  };

  template <class ChannelT, class PixelT>
  typename boost::disable_if< IsImageView<PixelT>, typename CompoundChannelCast<PixelT, ChannelT>::type >::type
  inline channel_cast_rescale( PixelT pixel ) {
    return compound_apply( ChannelCastRescaleFunctor<ChannelT>(), pixel );
  }


  // *******************************************************************
  // Pixel math base class and free functions.
  //
  // Here we define the PixelMathBase CRTP base class, which can be 
  // used as a mix-in to give any pixel-like type the usual pixel math 
  // behavior.
  // *******************************************************************

  /// This CRTP base class, and the related free operator overloads,
  /// provide some default pixel math behavior used by pixel types
  /// such as PixelRGB.  It makes extensive use of type computation
  /// mechanisms to handle things like type promotion, so that for
  /// example a PixelRGB<int> multiplied by a float returns a
  /// PixelRGB<float>.  The following operators are currently
  /// supported:
  /// 
  ///  - <TT>-</TT> pixel
  ///  - pixel <TT>+</TT> pixel
  ///  - pixel <TT>+=</TT> pixel
  ///  - pixel <TT>-</TT> pixel
  ///  - pixel <TT>-=</TT> pixel
  ///  - pixel <TT>*</TT> scalar
  ///  - scalar <TT>*</TT> pixel
  ///  - pixel <TT>*=</TT> scalar
  ///  - pixel <TT>/</TT> scalar
  ///  - pixel <TT>/=</TT> scalar
  ///
  /// To use this class, simply derive your pixel type from it in 
  /// the usual CRTP manner.  Your pixel type must provide operator[]
  /// and you must specialize the usual compound type traits:
  /// CompoundChannelType, CompoundNumChannels, and CompoundChannelCast.
  template <class DerivedT>
  class PixelMathBase {
  public:
    template <class ArgT> inline DerivedT& operator+=( ArgT a ) { return *static_cast<DerivedT*>(this) = (static_cast<DerivedT const&>(*this) + a ); }
    template <class ArgT> inline DerivedT& operator-=( ArgT a ) { return *static_cast<DerivedT*>(this) = (static_cast<DerivedT const&>(*this) - a ); }
    template <class ArgT> inline DerivedT& operator*=( ArgT a ) { return *static_cast<DerivedT*>(this) = (static_cast<DerivedT const&>(*this) * a ); }
    template <class ArgT> inline DerivedT& operator/=( ArgT a ) { return *static_cast<DerivedT*>(this) = (static_cast<DerivedT const&>(*this) / a ); }
  };

  template <class PixelT>
  inline PixelT operator-( PixelMathBase<PixelT> const& pixel ) {
    return compound_apply(ArgNegationFunctor(), static_cast<PixelT const&>(pixel) );
  }

  template <class Pixel1T, class Pixel2T>
  typename boost::enable_if< CompoundIsCompatible<Pixel1T,Pixel2T>, typename CompoundResult<ArgArgSumFunctor, Pixel1T, Pixel2T>::type >::type
  inline operator+( PixelMathBase<Pixel1T> const& pixel1, PixelMathBase<Pixel2T> const& pixel2 ) {
    return compound_apply(ArgArgSumFunctor(), static_cast<Pixel1T const&>(pixel1), static_cast<Pixel2T const&>(pixel2) );
  }

  template <class Pixel1T, class Pixel2T>
  typename boost::enable_if< CompoundIsCompatible<Pixel1T,Pixel2T>, typename CompoundResult<ArgArgDifferenceFunctor, Pixel1T, Pixel2T>::type >::type
  inline operator-( PixelMathBase<Pixel1T> const& pixel1, PixelMathBase<Pixel2T> const& pixel2 ) {
    return compound_apply(ArgArgDifferenceFunctor(), static_cast<Pixel1T const&>(pixel1), static_cast<Pixel2T const&>(pixel2) );
  }

  template <class PixelT, class ScalarT>
  typename boost::enable_if< IsScalar<ScalarT>, typename CompoundResult<ArgValProductFunctor<ScalarT>, PixelT>::type >::type
  inline operator*( PixelMathBase<PixelT> const& pixel, ScalarT scalar ) {
    return compound_apply(ArgValProductFunctor<ScalarT>(scalar), static_cast<PixelT const&>(pixel) );
  }

  template <class ScalarT, class PixelT>
  typename boost::enable_if< IsScalar<ScalarT>, typename CompoundResult<ValArgProductFunctor<ScalarT>, PixelT>::type >::type
  inline operator*( ScalarT scalar, PixelMathBase<PixelT> const& pixel ) {
    return compound_apply(ValArgProductFunctor<ScalarT>(scalar), static_cast<PixelT const&>(pixel) );
  }

  template <class PixelT, class ScalarT>
  typename boost::enable_if< IsScalar<ScalarT>, typename CompoundResult<ArgValQuotientFunctor<ScalarT>, PixelT>::type >::type
  inline operator/( PixelMathBase<PixelT> const& pixel, ScalarT scalar ) {
    return compound_apply(ArgValQuotientFunctor<ScalarT>(scalar), static_cast<PixelT const&>(pixel) );
  }


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
    static const unsigned value = NCHANNELS;                 \
  };                                                         \
  template <class OldChT, class NewChT>                      \
  struct CompoundChannelCast<PIXELT<OldChT>, NewChT> {       \
    typedef PIXELT<NewChT> type;                             \
  };                                                         \
  template <class OldChT, class NewChT>                      \
  struct CompoundChannelCast<PIXELT<OldChT>, const NewChT> { \
    typedef const PIXELT<NewChT> type;                       \
  }

  // Forward pixel type declarations for complex pixel types
  template <class ChannelT> class PixelGray;
  template <class ChannelT> class PixelGrayA;
  template <class ChannelT> class PixelRGB;
  template <class ChannelT> class PixelRGBA;
  template <class ChannelT> class PixelHSV;

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
    VW_CHANNEL_USER = 100
  };
  
  template <class PixelT> struct PixelFormatID { static const PixelFormatEnum value = VW_PIXEL_UNKNOWN; };
  template<> struct PixelFormatID<vw::int8>    { static const PixelFormatEnum value = VW_PIXEL_SCALAR; };
  template<> struct PixelFormatID<vw::uint8>   { static const PixelFormatEnum value = VW_PIXEL_SCALAR; };
  template<> struct PixelFormatID<vw::int16>   { static const PixelFormatEnum value = VW_PIXEL_SCALAR; };
  template<> struct PixelFormatID<vw::uint16>  { static const PixelFormatEnum value = VW_PIXEL_SCALAR; };
  template<> struct PixelFormatID<vw::int32>   { static const PixelFormatEnum value = VW_PIXEL_SCALAR; };
  template<> struct PixelFormatID<vw::uint32>  { static const PixelFormatEnum value = VW_PIXEL_SCALAR; };
  template<> struct PixelFormatID<vw::int64>   { static const PixelFormatEnum value = VW_PIXEL_SCALAR; };
  template<> struct PixelFormatID<vw::uint64>  { static const PixelFormatEnum value = VW_PIXEL_SCALAR; };
  template<> struct PixelFormatID<float>       { static const PixelFormatEnum value = VW_PIXEL_SCALAR; };
  template<> struct PixelFormatID<double>      { static const PixelFormatEnum value = VW_PIXEL_SCALAR; };
  template<> struct PixelFormatID<bool>        { static const PixelFormatEnum value = VW_PIXEL_SCALAR; };
  template <class ChT> struct PixelFormatID<PixelGray<ChT> >  { static const PixelFormatEnum value = VW_PIXEL_GRAY; };
  template <class ChT> struct PixelFormatID<PixelGrayA<ChT> > { static const PixelFormatEnum value = VW_PIXEL_GRAYA; };
  template <class ChT> struct PixelFormatID<PixelRGB<ChT> >   { static const PixelFormatEnum value = VW_PIXEL_RGB; };
  template <class ChT> struct PixelFormatID<PixelRGBA<ChT> >  { static const PixelFormatEnum value = VW_PIXEL_RGBA; };
  template <class ChT> struct PixelFormatID<PixelHSV<ChT> >   { static const PixelFormatEnum value = VW_PIXEL_HSV; };

  template <class ChannelT> struct ChannelTypeID { static const ChannelTypeEnum value = VW_CHANNEL_UNKNOWN; };
  template<> struct ChannelTypeID<vw::int8>      { static const ChannelTypeEnum value = VW_CHANNEL_INT8; };
  template<> struct ChannelTypeID<vw::uint8>     { static const ChannelTypeEnum value = VW_CHANNEL_UINT8; };
  template<> struct ChannelTypeID<vw::int16>     { static const ChannelTypeEnum value = VW_CHANNEL_INT16; };
  template<> struct ChannelTypeID<vw::uint16>    { static const ChannelTypeEnum value = VW_CHANNEL_UINT16; };
  template<> struct ChannelTypeID<vw::int32>     { static const ChannelTypeEnum value = VW_CHANNEL_INT32; };
  template<> struct ChannelTypeID<vw::uint32>    { static const ChannelTypeEnum value = VW_CHANNEL_UINT32; };
  template<> struct ChannelTypeID<vw::int64>     { static const ChannelTypeEnum value = VW_CHANNEL_INT64; };
  template<> struct ChannelTypeID<vw::uint64>    { static const ChannelTypeEnum value = VW_CHANNEL_UINT64; };
  template<> struct ChannelTypeID<float>         { static const ChannelTypeEnum value = VW_CHANNEL_FLOAT32; };
  template<> struct ChannelTypeID<double>        { static const ChannelTypeEnum value = VW_CHANNEL_FLOAT64; };
  template<> struct ChannelTypeID<bool>          { static const ChannelTypeEnum value = VW_CHANNEL_BOOL; };

  inline int num_channels( PixelFormatEnum format ) {
    switch( format ) {
    case VW_PIXEL_SCALAR:
    case VW_PIXEL_GRAY:
      return 1;
    case VW_PIXEL_GRAYA:
      return 2;
    case VW_PIXEL_RGB:
    case VW_PIXEL_HSV:
      return 3;
    case VW_PIXEL_RGBA:
      return 4;
    default:
      throw ArgumentErr() << "Unrecognized or unsupported pixel format (" << format << ").";
    }
  }

  inline int channel_size( ChannelTypeEnum type ) {
    switch( type ) {
    case VW_CHANNEL_BOOL:
      return sizeof(bool);
    case VW_CHANNEL_INT8:
    case VW_CHANNEL_UINT8:
      return 1;
    case VW_CHANNEL_INT16:
    case VW_CHANNEL_UINT16:
      return 2;
    case VW_CHANNEL_INT32:
    case VW_CHANNEL_UINT32:
    case VW_CHANNEL_FLOAT32:
      return 4;
    case VW_CHANNEL_INT64:
    case VW_CHANNEL_UINT64:
    case VW_CHANNEL_FLOAT64:
      return 8;
    default:
      throw ArgumentErr() << "Unrecognized or unsupported channel type (" << type << ").";
    }
  }

} // namespace vw

#endif // __VW_IMAGE__PIXEL_TYPE_INFO_H__
