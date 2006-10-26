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

/// \file PixelTypes.h
/// 
/// Defines the standard pixel types.  Currently this includes
/// PixelGray, PixelGrayA, PixelRGB, PixelRGBA, PixelHSV, and
/// PixelXYZ.
///
/// The alpha pixel types (PixelGrayA and PixelRGBA) have 
/// arithmetic operators defined that assume the pixels 
/// are represented in pre-multiplied form.
///
#ifndef __VW_IMAGE_PIXELTYPES_H__
#define __VW_IMAGE_PIXELTYPES_H__

#include <ostream>

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

#include <vw/Image/PixelTypeInfo.h>
#include <vw/Image/PixelMath.h>
#include <vw/Math/Vector.h>

namespace vw {

  /// \cond INTERNAL

  // Default weights used in the weighted_rgb_to_gray() functions
#define VW_RGB_TO_GRAY_R_WEIGHT (0.30)
#define VW_RGB_TO_GRAY_G_WEIGHT (0.59)
#define VW_RGB_TO_GRAY_B_WEIGHT (0.11)

  // Default weights used for RGB to XYZ conversion (D65 white point)
#define VW_RGB_TO_XYZ_RX_WEIGHT (0.412453)
#define VW_RGB_TO_XYZ_RY_WEIGHT (0.212671)
#define VW_RGB_TO_XYZ_RZ_WEIGHT (0.019334)
#define VW_RGB_TO_XYZ_GX_WEIGHT (0.357580)
#define VW_RGB_TO_XYZ_GY_WEIGHT (0.715160)
#define VW_RGB_TO_XYZ_GZ_WEIGHT (0.119193)
#define VW_RGB_TO_XYZ_BX_WEIGHT (0.180423)
#define VW_RGB_TO_XYZ_BY_WEIGHT (0.072169)
#define VW_RGB_TO_XYZ_BZ_WEIGHT (0.950227)

  // Default weights used for XYZ to RGB conversion (D65 white point)
#define VW_XYZ_TO_RGB_XR_WEIGHT (3.240479)
#define VW_XYZ_TO_RGB_XG_WEIGHT (-0.969256)
#define VW_XYZ_TO_RGB_XB_WEIGHT (0.055648)
#define VW_XYZ_TO_RGB_YR_WEIGHT (-1.537150)
#define VW_XYZ_TO_RGB_YG_WEIGHT (1.875992)
#define VW_XYZ_TO_RGB_YB_WEIGHT (-0.204043)
#define VW_XYZ_TO_RGB_ZR_WEIGHT (-0.498535)
#define VW_XYZ_TO_RGB_ZG_WEIGHT (0.041556)
#define VW_XYZ_TO_RGB_ZB_WEIGHT (1.057311)

  // This function is used in messy color space conversions 
  // that operate internally using floating-point numbers 
  // even when the destination channel type is integral.
  template <class T> inline T _round_if_needed( typename boost::enable_if<boost::is_floating_point<T>,double>::type v) { return v; }
  template <class T> inline T _round_if_needed( typename boost::disable_if<boost::is_floating_point<T>,double>::type v) {
    double vr = round(v);
    if( vr < ChannelRange<T>::min() ) return ChannelRange<T>::min();
    if( vr > ChannelRange<T>::max() ) return ChannelRange<T>::max();
    return vr;
  }

  /// \endcond

  // *******************************************************************
  // The PixelGray grayscale pixel type.
  // *******************************************************************

  /// A grayscale pixel type with only one channel, the pixel's value 
  /// or luminance.  This pixel type can be mixed reasonably freely 
  /// with pure scalars in mathematical expressions.
  template <class ChannelT>
  class PixelGray : public PixelMathBase< PixelGray<ChannelT> >
  {
    ChannelT m_ch[1];

  public:
    // Default constructor (zero value).
    PixelGray() { m_ch[0]=0; }

    /// Implicit construction from the raw channel value.
    PixelGray( ChannelT v ) { m_ch[0]=v; }

    /// Explicit conversion from other PixelGray types.
    template <class OtherT> explicit PixelGray( PixelGray<OtherT> other ) {
      m_ch[0] = ChannelT(other[0]);
    }
    
    /// Explicit conversion from PixelGrayA types.  Discards the alpha
    /// channel, which is equivalent to compositing over a black
    /// background.
    template <class OtherT> explicit PixelGray( PixelGrayA<OtherT> other ) {
      m_ch[0]=ChannelT(other[0]);
    }

    /// Explicit conversion from PixelRGB types.  Computes the
    /// luminance value by averaging the R, G, and B channels.
    template <class OtherT> explicit PixelGray( PixelRGB<OtherT> other ) {
      typedef typename AccumulatorType<ChannelT>::type a_t;
      m_ch[0] = ChannelT( ( a_t(other[0]) + a_t(other[1]) + a_t(other[2]) ) / 3 );
    }

    /// Explicit conversion from PixelRGBA types.  Computes the
    /// luminance value by averaging the R, G, and B channels.
    /// Discards the alpha channel, which is equivalent to compositing
    /// over a black background.
    template <class OtherT> explicit PixelGray( PixelRGBA<OtherT> other ) {
      typedef typename AccumulatorType<ChannelT>::type a_t;
      m_ch[0] = ChannelT( ( a_t(other[0]) + a_t(other[1]) + a_t(other[2]) ) / 3 );
    }

    /// Automatic down-cast to the raw channel value in numeric contexts.
    operator ChannelT() const { return m_ch[0]; }

    /// Channel indexing operator.
    inline ChannelT& operator[](int i) { return m_ch[i]; }
    /// Channel indexing operator (const overload).
    inline ChannelT const& operator[](int i) const { return m_ch[i]; }
    /// Channel indexing operator.
    inline ChannelT& operator()(int i) { return m_ch[i]; }
    /// Channel indexing operator (const overload).
    inline ChannelT const& operator()(int i) const { return m_ch[i]; }

    /// Value channel accessor.
    inline ChannelT& v() { return m_ch[0]; }
    /// Value channel accessor (const overload).
    inline ChannelT const& v() const { return m_ch[0]; }
  };

  /// Declare the standard 1-channel pixel traits for PixelGray.
  VW_DECLARE_PIXEL_TYPE(PixelGray,1);

  /// Print a PixelGrayA to a debugging stream.
  template <class ChannelT>
  std::ostream& operator<<( std::ostream& os, PixelGray<ChannelT> const& pix ) {
    return os << "PixelGray(" << _numeric(pix.v()) << ")";
  }


  // *******************************************************************
  // The PixelGrayA grayscale/alpha pixel type.
  // *******************************************************************

  /// A grayscale pixel type with an alpha channel.
  template <class ChannelT>
  class PixelGrayA : public PixelMathBase< PixelGrayA<ChannelT> >
  {
    ChannelT m_ch[2];

  public:
    // Default constructor (zero value, i.e. transparent).
    PixelGrayA() { m_ch[0]=m_ch[1]=0; }

    /// Constructs an opaque pixel from the raw luminance value.  This
    /// is marked as explicit to prevent you from accidentially
    /// initializing a pixel to zero when what you really want is the
    /// default-constructed value.
    explicit PixelGrayA( ChannelT v ) { 
      m_ch[0]=v;
      m_ch[1]=ChannelRange<ChannelT>::max();
    }

    /// Constructs a pixel with the given channel values.
    PixelGrayA( ChannelT v, ChannelT a ) {
      m_ch[0]=v;
      m_ch[1]=a;
    }

    /// Explicit conversion from PixelGray types.
    template <class OtherT> explicit PixelGrayA( PixelGray<OtherT> other ) {
      m_ch[0] = ChannelT(other[0]);
      m_ch[1] = ChannelRange<ChannelT>::max();
    }

    /// Explicit conversion from other PixelGrayA types.
    template <class OtherT> explicit PixelGrayA( PixelGrayA<OtherT> other ) {
      m_ch[0] = ChannelT(other[0]);
      m_ch[1] = ChannelT(other[1]);
    }

    /// Explicit conversion from PixelRGB types.  Computes the
    /// luminance value by averaging the R, G, and B channels.
    template <class OtherT> explicit PixelGrayA( PixelRGB<OtherT> other ) {
      typedef typename AccumulatorType<ChannelT>::type a_t;
      m_ch[0] = ChannelT( ( a_t(other[0]) + a_t(other[1]) + a_t(other[2]) ) / 3 );
      m_ch[1] = ChannelRange<ChannelT>::max();
    }

    /// Explicit conversion from PixelRGBA types.  Computes the
    /// luminance value by averaging the R, G, and B channels.
    template <class OtherT> explicit PixelGrayA( PixelRGBA<OtherT> other ) {
      typedef typename AccumulatorType<ChannelT>::type a_t;
      m_ch[0] = ChannelT( ( a_t(other[0]) + a_t(other[1]) + a_t(other[2]) ) / 3 );
      m_ch[1] = ChannelT(other[3]);
    }

    /// Channel indexing operator.
    inline ChannelT& operator[](int i) { return m_ch[i]; }
    /// Channel indexing operator (const overload).
    inline ChannelT const& operator[](int i) const { return m_ch[i]; }
    /// Channel indexing operator.
    inline ChannelT& operator()(int i) { return m_ch[i]; }
    /// Channel indexing operator (const overload).
    inline ChannelT const& operator()(int i) const { return m_ch[i]; }

    /// Value channel accessor.
    inline ChannelT& v() { return m_ch[0]; }
    /// Value channel accessor (const overload).
    inline ChannelT const& v() const { return m_ch[0]; }
    /// Alpha channel accessor.
    inline ChannelT& a() { return m_ch[1]; }
    /// Alpha channel accessor (const overload).
    inline ChannelT const& a() const { return m_ch[1]; }
  };

  /// Declare the standard 2-channel pixel traits for PixelGrayA.
  VW_DECLARE_PIXEL_TYPE(PixelGrayA,2);
  template <class T> struct PixelHasAlpha<PixelGrayA<T> > : true_type {};

  /// Print a PixelGrayA to a debugging stream.
  template <class ChannelT>
  std::ostream& operator<<( std::ostream& os, PixelGrayA<ChannelT> const& pix ) {
    return os << "PixelGrayA(" << _numeric(pix.v()) << "," << _numeric(pix.a()) << ")";
  }


  // *******************************************************************
  // The PixelRGB red/green/blue pixel type.
  // *******************************************************************

  /// A full-color RGB pixel type.
  template <class ChannelT>
  class PixelRGB : public PixelMathBase< PixelRGB<ChannelT> >
  {
    ChannelT m_ch[3];

  public:
    // Default constructor (zero value).
    PixelRGB() { m_ch[0]=m_ch[1]=m_ch[2]=0; }

    /// Explicitly constructs of a gray pixel with the given luminance value.
    PixelRGB( ChannelT const& v ) { m_ch[0]=m_ch[1]=m_ch[2]=v; };

    /// Constructs a pixel with the given channel values.
    PixelRGB( ChannelT const& r, ChannelT const& g, ChannelT const& b ) {
      m_ch[0]=r; m_ch[1]=g; m_ch[2]=b;
    }

    /// Explicit conversion from PixelGray types.
    template <class OtherT> explicit PixelRGB( PixelGray<OtherT> const& other ) {
      m_ch[0] = m_ch[1] = m_ch[2] = ChannelT(other[0]);
    }
    
    /// Explicit conversion from PixelGrayA types.
    template <class OtherT> explicit PixelRGB( PixelGrayA<OtherT> const& other ) {
      m_ch[0] = m_ch[1] = m_ch[2] = ChannelT(other[0]);
    }
    
    /// Explicit conversion from other PixelRGB types.
    template <class OtherT> explicit PixelRGB( PixelRGB<OtherT> const& other ) {
      m_ch[0] = ChannelT(other[0]);
      m_ch[1] = ChannelT(other[1]);
      m_ch[2] = ChannelT(other[2]);
    }
    
    /// Explicit conversion from PixelRGBA types.
    template <class OtherT> explicit PixelRGB( PixelRGBA<OtherT> const& other ) {
      m_ch[0] = ChannelT(other[0]);
      m_ch[1] = ChannelT(other[1]);
      m_ch[2] = ChannelT(other[2]);
    }
    
    /// Explicit conversion from PixelHSV types.
    template <class OtherT> explicit PixelRGB( PixelHSV<OtherT> const& other );

    /// Explicit conversion from PixelXYZ types.
    template <class OtherT> explicit PixelRGB( PixelXYZ<OtherT> const& other );

    /// Channel indexing operator.
    inline ChannelT& operator[](int i) { return m_ch[i]; }
    /// Channel indexing operator (const overload).
    inline ChannelT const& operator[](int i) const { return m_ch[i]; }
    /// Channel indexing operator.
    inline ChannelT& operator()(int i) { return m_ch[i]; }
    /// Channel indexing operator (const overload).
    inline ChannelT const& operator()(int i) const { return m_ch[i]; }

    /// Red channel accessor.
    inline ChannelT& r() { return m_ch[0]; }
    /// Red channel accessor (const overload).
    inline ChannelT const& r() const { return m_ch[0]; }
    /// Green channel accessor.
    inline ChannelT& g() { return m_ch[1]; }
    /// Green channel accessor (const overload).
    inline ChannelT const& g() const { return m_ch[1]; }
    /// Blue channel accessor.
    inline ChannelT& b() { return m_ch[2]; }
    /// Blue channel accessor (const overload).
    inline ChannelT const& b() const { return m_ch[2]; }
  };

  /// Declare the standard 3-channel pixel traits for PixelRGB.
  VW_DECLARE_PIXEL_TYPE(PixelRGB,3);

  /// Print a PixelRGB to a debugging stream.
  template <class ChannelT>
  std::ostream& operator<<( std::ostream& os, PixelRGB<ChannelT> const& pix ) {
    return os << "PixelRGB(" << _numeric(pix.r()) << "," << _numeric(pix.g()) << "," << _numeric(pix.b()) << ")";
  }

  /// Weighted conversion from PixelRGB to PixelGray using user-specified weights.
  template <class ChannelT>
  PixelGray<ChannelT> weighted_rgb_to_gray( PixelRGB<ChannelT> const& rgb, double rw, double gw, double bw ) {
    return PixelGray<ChannelT>( rgb.r()*rw + rgb.g()*gw + rgb.b()*bw );
  }

  /// Weighted conversion from PixelRGB to PixelGray using the default weights.
  template <class ChannelT>
  PixelGray<ChannelT> weighted_rgb_to_gray( PixelRGB<ChannelT> const& rgb ) {
    return weighted_rgb_to_gray( rgb, VW_RGB_TO_GRAY_R_WEIGHT, VW_RGB_TO_GRAY_G_WEIGHT, VW_RGB_TO_GRAY_B_WEIGHT );
  }


  // *******************************************************************
  // The PixelRGBA red/green/blue/alpha pixel type.
  // *******************************************************************

  /// A full-color RGB pixel type with an alpha channel.
  template <class ChannelT>
  class PixelRGBA : public PixelMathBase< PixelRGBA<ChannelT> >
  {
    ChannelT m_ch[4];

  public:
    // Default constructor (zero value, i.e. transparent).
    PixelRGBA() { m_ch[0]=m_ch[1]=m_ch[2]=m_ch[3]=0; }

    /// Constructs an opaque gray pixel from the raw luminance value.
    /// This is marked as explicit to prevent you from accidentially
    /// initializing a pixel to zero when what you really want is the
    /// default-constructed value.
    explicit PixelRGBA( ChannelT const& v ) { 
      m_ch[0]=m_ch[1]=m_ch[2]=v;
      m_ch[3]=ChannelRange<ChannelT>::max();
    }

    /// Constructs a pixel with the given channel values.
    PixelRGBA( ChannelT const& r, ChannelT const& g, ChannelT const& b, ChannelT const& a ) {
      m_ch[0]=r; m_ch[1]=g; m_ch[2]=b; m_ch[3]=a;
    }

    /// Explicit conversion from PixelGray types.
    template <class OtherT> explicit PixelRGBA( PixelGray<OtherT> const& other ) {
      m_ch[0] = m_ch[1] = m_ch[2] = ChannelT(other[0]);
      m_ch[3]=ChannelRange<ChannelT>::max();
    }
    
    /// Explicit conversion from PixelGrayA types.
    template <class OtherT> explicit PixelRGBA( PixelGrayA<OtherT> const& other ) {
      m_ch[0] = m_ch[1] = m_ch[2] = ChannelT(other[0]);
      m_ch[3]=ChannelT(other[1]);
    }
    
    /// Explicit conversion from PixelRGB types.
    template <class OtherT> explicit PixelRGBA( PixelRGB<OtherT> const& other ) {
      m_ch[0] = ChannelT(other[0]);
      m_ch[1] = ChannelT(other[1]);
      m_ch[2] = ChannelT(other[2]);
      m_ch[3] = ChannelRange<ChannelT>::max();
    }

    /// Explicit conversion from other PixelRGBA types.    
    template <class OtherT> explicit PixelRGBA( PixelRGBA<OtherT> const& other ) {
      m_ch[0] = ChannelT(other[0]);
      m_ch[1] = ChannelT(other[1]);
      m_ch[2] = ChannelT(other[2]);
      m_ch[3] = ChannelT(other[3]);
    }

    /// Channel indexing operator.
    inline ChannelT& operator[](int i) { return m_ch[i]; }
    /// Value channel accessor (const overload).
    inline ChannelT const& operator[](int i) const { return m_ch[i]; }
    /// Channel indexing operator.
    inline ChannelT& operator()(int i) { return m_ch[i]; }
    /// Value channel accessor (const overload).
    inline ChannelT const& operator()(int i) const { return m_ch[i]; }

    /// Red channel accessor.
    inline ChannelT& r() { return m_ch[0]; }
    /// Red channel accessor (const overload).
    inline ChannelT const& r() const { return m_ch[0]; }
    /// Green channel accessor.
    inline ChannelT& g() { return m_ch[1]; }
    /// Green channel accessor (const overload).
    inline ChannelT const& g() const { return m_ch[1]; }
    /// Blue channel accessor.
    inline ChannelT& b() { return m_ch[2]; }
    /// Blue channel accessor (const overload).
    inline ChannelT const& b() const { return m_ch[2]; }
    /// Alpha channel accessor.
    inline ChannelT& a() { return m_ch[3]; }
    /// Alpha channel accessor (const overload).
    inline ChannelT const& a() const { return m_ch[3]; }
  };

  /// Declare the standard 2-channel pixel traits for PixelGrayA.
  VW_DECLARE_PIXEL_TYPE(PixelRGBA,4);
  template <class T> struct PixelHasAlpha<PixelRGBA<T> > : true_type {};

  /// Print a PixelRGBA to a debugging stream.
  template <class ChannelT>
  std::ostream& operator<<( std::ostream& os, PixelRGBA<ChannelT> const& pix ) {
    return os << "RGBA(" << _numeric(pix.r()) << "," << _numeric(pix.g()) << "," << _numeric(pix.b()) << "," << _numeric(pix.a()) << ")";
  }

  /// Weighted conversion from PixelRGBA to PixelGrayA using user-specified weights.
  template <class ChannelT>
  PixelGrayA<ChannelT> weighted_rgb_to_gray( PixelRGBA<ChannelT> const& rgba, double rw, double gw, double bw ) {
    return PixelGrayA<ChannelT>( rgba.r()*rw + rgba.g()*gw + rgba.b()*bw, rgba.a() );
  }

  /// Weighted conversion from PixelRGBA to PixelGrayA using the default weights.
  template <class ChannelT>
  PixelGrayA<ChannelT> weighted_rgb_to_gray( PixelRGBA<ChannelT> const& rgba ) {
    return weighted_rgb_to_gray( rgba, VW_RGB_TO_GRAY_R_WEIGHT, VW_RGB_TO_GRAY_G_WEIGHT, VW_RGB_TO_GRAY_B_WEIGHT );
  }


  // *******************************************************************
  // The PixelHSV hue/saturation/value pixel type.
  // *******************************************************************

  /// A full-color HSV (hue, saturation, value) pixel type.
  /// Note that this color space is non-linearly related to 
  /// the usual color spaces such as RGB, and has a different 
  /// topology as well.  Therefore many operations, such as 
  /// pixel math and image filtering, may yield surprising 
  /// results if you use this pixel type.  
  template <class ChannelT>
  class PixelHSV : public PixelMathBase< PixelHSV<ChannelT> >
  {
    ChannelT m_ch[3];

  public:
    // Default constructor (zero value).
    PixelHSV() { m_ch[0]=m_ch[1]=m_ch[2]=0; }

    /// Explicitly constructs a gray pixel with the given luminance value.
    PixelHSV( ChannelT const& v ) { m_ch[0]=m_ch[1]=0; m_ch[2]=v; }

    /// Constructs a pixel with the given channel values.
    PixelHSV( ChannelT const& h, ChannelT const& s, ChannelT const& v ) {
      m_ch[0]=h; m_ch[1]=s; m_ch[2]=v;
    }

    /// Explicit conversion from other PixelHSV types.
    template <class OtherT> explicit PixelHSV( PixelHSV<OtherT> const& other ) {
      m_ch[0] = ChannelT(other[0]);
      m_ch[1] = ChannelT(other[1]);
      m_ch[2] = ChannelT(other[2]);
    }

    /// Explicit conversion from PixelRGB types.
    template <class OtherT> explicit PixelHSV( PixelRGB<OtherT> const& rgb );

    /// Channel indexing operator.
    inline ChannelT& operator[](int i) { return m_ch[i]; }
    /// Channel indexing operator (const overload).
    inline ChannelT const& operator[](int i) const { return m_ch[i]; }
    /// Channel indexing operator.
    inline ChannelT& operator()(int i) { return m_ch[i]; }
    /// Channel indexing operator (const overload).
    inline ChannelT const& operator()(int i) const { return m_ch[i]; }

    /// Hue channel accessor.
    inline ChannelT& h() { return m_ch[0]; }
    /// Hue channel accessor (const overload).
    inline ChannelT const& h() const { return m_ch[0]; }
    /// Saturation channel accessor.
    inline ChannelT& s() { return m_ch[1]; }
    /// Saturation channel accessor (const overload).
    inline ChannelT const& s() const { return m_ch[1]; }
    /// Value channel accessor.
    inline ChannelT& v() { return m_ch[2]; }
    /// Value channel accessor (const overload).
    inline ChannelT const& v() const { return m_ch[2]; }
  };

  /// Declare the standard 3-channel pixel traits for PixelHSV.
  VW_DECLARE_PIXEL_TYPE(PixelHSV, 3);

  /// Print a PixelHSV to a debugging stream.
  template <class ChannelT>
  std::ostream& operator<<( std::ostream& os, PixelHSV<ChannelT> const& pix ) {
    return os << "PixelHSV(" << _numeric(pix.h()) << "," << _numeric(pix.s()) << "," << _numeric(pix.v()) << ")";
  }

  /// HSV->RGB conversion
  template <class ChannelT>
  template <class OtherT>
  PixelRGB<ChannelT>::PixelRGB( PixelHSV<OtherT> const& hsv ) {
    double sscale = ChannelRange<ChannelT>::max();
    double hscale = (boost::is_integral<ChannelT>::value) ? (sscale+1) : sscale;
    double s = hsv.s()/sscale;
    double h = hsv.h()/hscale;
    int i = (int)floor(6*h);
    double f = 6*h - i;
    if(!(i & 1)) f = 1 - f; // if i is even
    ChannelT m = _round_if_needed<ChannelT>( hsv.v() * (1 - s) );
    ChannelT n = _round_if_needed<ChannelT>( hsv.v() * (1 - s*f) );
    switch (i) {
    case 6:
    case 0: r()=hsv.v(); g()=n; b()=m; break;
    case 1: r()=n; g()=hsv.v(); b()=m; break;
    case 2: r()=m; g()=hsv.v(); b()=n; break;
    case 3: r()=m; g()=n; b()=hsv.v(); break;
    case 4: r()=n; g()=m; b()=hsv.v(); break;
    case 5: r()=hsv.v(); g()=m; b()=n; break;
    }
  }

  /// RGB->HSV conversion
  template <class ChannelT>
  template <class OtherT>
  PixelHSV<ChannelT>::PixelHSV( PixelRGB<OtherT> const& rgb ) {
    double sscale = ChannelRange<ChannelT>::max();
    double hscale = (boost::is_integral<ChannelT>::value) ? (sscale+1) : sscale;
    if( rgb.r()==rgb.g() && rgb.r()==rgb.b() ) { // S=0
      m_ch[0]=ChannelT();
      m_ch[1]=ChannelT();
      m_ch[2]=rgb.r();
    }
    else if( rgb.r()<=rgb.g() && rgb.r()<=rgb.b() ) { // H=[1/3,2/3)
      m_ch[2] = (rgb.g()>rgb.b()) ? rgb.g() : rgb.b();
      m_ch[0] = _round_if_needed<ChannelT>( (3 - (rgb.g()-(double)rgb.b())/(m_ch[2]-rgb.r()))*hscale/6.0 );
      m_ch[1] = _round_if_needed<ChannelT>( ((m_ch[2]-rgb.r())/(double)m_ch[2])*sscale );
    }
    else if( rgb.g()<=rgb.r() && rgb.g()<rgb.b() ) { // H=[2/3,1)
      m_ch[2] = (rgb.r()>rgb.b()) ? rgb.r() : rgb.b();
      m_ch[0] = _round_if_needed<ChannelT>( (5 - (rgb.b()-(double)rgb.r())/(m_ch[2]-rgb.g()))*hscale/6.0 );
      m_ch[1] = _round_if_needed<ChannelT>( ((m_ch[2]-rgb.g())/(double)m_ch[2])*sscale );
    }
    else { // H=[0,1/3)
      m_ch[2] = (rgb.r()>rgb.g()) ? rgb.r() : rgb.g();
      m_ch[0] = _round_if_needed<ChannelT>( (1 - (rgb.r()-(double)rgb.g())/(m_ch[2]-rgb.b()))*hscale/6.0 );
      m_ch[1] = _round_if_needed<ChannelT>( ((m_ch[2]-rgb.b())/(double)m_ch[2])*sscale );
    }
  }


  // *******************************************************************
  // The PixelXYZ CIE 1931 XYZ linear perceptual pixel type.
  // *******************************************************************

  /// An full-color XYZ color space pixel type.  This color space is
  /// linearly related to RGB, but the linear transformation (as
  /// defined by the International Commission on Illumination (CIE) in
  /// 1931) is nonuniform.  When converting between XYZ and RGB,
  /// values are clamped to the usual ranges when using integral pixel
  /// types to avoid wrap-around.  Values are not clamped when using
  /// floating-point pixel types to avoid color distortion.
  template <class ChannelT>
  class PixelXYZ : public PixelMathBase< PixelXYZ<ChannelT> >
  {
    ChannelT m_ch[3];

  public:
    // Default constructor (zero value).
    PixelXYZ() { m_ch[0]=m_ch[1]=m_ch[2]=0; }

    /// Explicitly constructs a gray pixel with the given Y value.
    PixelXYZ( ChannelT const& y ) { m_ch[0]=m_ch[2]=0; m_ch[1]=y; }

    /// Constructs a pixel with the given channel values.
    PixelXYZ( ChannelT const& x, ChannelT const& y, ChannelT const& z ) {
      m_ch[0]=x; m_ch[1]=y; m_ch[2]=z;
    }

    /// Explicit conversion from other PixelXYZ types.
    template <class OtherT> explicit PixelXYZ( PixelXYZ<OtherT> const& other ) {
      m_ch[0] = ChannelT(other[0]);
      m_ch[1] = ChannelT(other[1]);
      m_ch[2] = ChannelT(other[2]);
    }

    /// Explicit conversion from PixelRGB types.
    template <class OtherT> explicit PixelXYZ( PixelRGB<OtherT> const& rgb );

    /// Channel indexing operator.
    inline ChannelT& operator[](int i) { return m_ch[i]; }
    /// Channel indexing operator (const overload).
    inline ChannelT const& operator[](int i) const { return m_ch[i]; }
    /// Channel indexing operator.
    inline ChannelT& operator()(int i) { return m_ch[i]; }
    /// Channel indexing operator (const overload).
    inline ChannelT const& operator()(int i) const { return m_ch[i]; }

    /// X channel accessor.
    inline ChannelT& x() { return m_ch[0]; }
    /// X channel accessor (const overload).
    inline ChannelT const& x() const { return m_ch[0]; }
    /// Y channel accessor.
    inline ChannelT& y() { return m_ch[1]; }
    /// Y channel accessor (const overload).
    inline ChannelT const& y() const { return m_ch[1]; }
    /// Z channel accessor.
    inline ChannelT& z() { return m_ch[2]; }
    /// Z channel accessor (const overload).
    inline ChannelT const& z() const { return m_ch[2]; }
  };

  /// Declare the standard 3-channel pixel traits for PixelXYZ.
  VW_DECLARE_PIXEL_TYPE(PixelXYZ, 3);

  /// Print a PixelXYZ to a debugging stream.
  template <class ChannelT>
  std::ostream& operator<<( std::ostream& os, PixelXYZ<ChannelT> const& pix ) {
    return os << "PixelXYZ(" << _numeric(pix.x()) << "," << _numeric(pix.y()) << "," << _numeric(pix.z()) << ")";
  }

  /// XYZ->RGB conversion
  template <class ChannelT>
  template <class OtherT>
  PixelRGB<ChannelT>::PixelRGB( PixelXYZ<OtherT> const& xyz ) {
    r() = _round_if_needed<ChannelT>( VW_XYZ_TO_RGB_XR_WEIGHT*xyz.x() +
                                      VW_XYZ_TO_RGB_YR_WEIGHT*xyz.y() +
                                      VW_XYZ_TO_RGB_ZR_WEIGHT*xyz.z() );
    g() = _round_if_needed<ChannelT>( VW_XYZ_TO_RGB_XG_WEIGHT*xyz.x() +
                                      VW_XYZ_TO_RGB_YG_WEIGHT*xyz.y() +
                                      VW_XYZ_TO_RGB_ZG_WEIGHT*xyz.z() );
    b() = _round_if_needed<ChannelT>( VW_XYZ_TO_RGB_XB_WEIGHT*xyz.x() +
                                      VW_XYZ_TO_RGB_YB_WEIGHT*xyz.y() +
                                      VW_XYZ_TO_RGB_ZB_WEIGHT*xyz.z() );
  }

  /// RGB->XYZ conversion
  template <class ChannelT>
  template <class OtherT>
  PixelXYZ<ChannelT>::PixelXYZ( PixelRGB<OtherT> const& rgb ) {
    x() = _round_if_needed<ChannelT>( VW_RGB_TO_XYZ_RX_WEIGHT*rgb.r() +
                                      VW_RGB_TO_XYZ_GX_WEIGHT*rgb.g() +
                                      VW_RGB_TO_XYZ_BX_WEIGHT*rgb.b() );
    y() = _round_if_needed<ChannelT>( VW_RGB_TO_XYZ_RY_WEIGHT*rgb.r() +
                                      VW_RGB_TO_XYZ_GY_WEIGHT*rgb.g() +
                                      VW_RGB_TO_XYZ_BY_WEIGHT*rgb.b() );
    z() = _round_if_needed<ChannelT>( VW_RGB_TO_XYZ_RZ_WEIGHT*rgb.r() +
                                      VW_RGB_TO_XYZ_GZ_WEIGHT*rgb.g() +
                                      VW_RGB_TO_XYZ_BZ_WEIGHT*rgb.b() );
  }


  // *******************************************************************
  // The Vector mathemaical vector pixel type.
  // *******************************************************************

  /// Declare the standard variable-channel pixel traits for the 
  /// Vector class when used as a pixel type.
  VW_DECLARE_PIXEL_TYPE_NCHANNELS(Vector);

}

#endif // __VW_IMAGE_PIXELTYPES_H__
