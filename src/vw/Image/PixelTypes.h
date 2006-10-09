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
/// Defines the five standard pixel types: PixelGray, PixelGrayA,
/// PixelRGB, PixelRGBA, and PixelHSV.
///
/// The alpha pixel types (PixelGrayA and PixelRGBA) have 
/// arithmetic operators defined that assume the pixels 
/// are represented in pre-multiplied form.
///
#ifndef __VW_IMAGE_PIXELTYPES_H__
#define __VW_IMAGE_PIXELTYPES_H__

#include <ostream>

#include <vw/Image/PixelTypeInfo.h>
#include <vw/Image/PixelMath.h>
#include <vw/Math/Vector.h>

namespace vw {

  // *******************************************************************
  // The PixelGray grayscale pixel type.
  // *******************************************************************

  /// Gray pixel type
  template <class ChannelT>
  class PixelGray : public PixelMathBase< PixelGray<ChannelT> >
  {
  private:
    ChannelT m_ch[1];
  public:
    PixelGray() { m_ch[0]=0; }
    PixelGray( ChannelT v ) { m_ch[0]=v; }

    template <class OtherT> explicit PixelGray( PixelGray<OtherT> other ) {
      m_ch[0] = ChannelT(other[0]);
    }
    
    template <class OtherT> explicit PixelGray( PixelGrayA<OtherT> other ) {
      m_ch[0]=ChannelT(other[0]);
    }

    template <class OtherT> explicit PixelGray( PixelRGB<OtherT> other ) {
      typedef typename AccumulatorType<ChannelT>::type a_t;
      m_ch[0] = ChannelT( ( a_t(other[0]) + a_t(other[1]) + a_t(other[2]) ) / 3 );
    }
    template <class OtherT> explicit PixelGray( PixelRGBA<OtherT> other ) {
      typedef typename AccumulatorType<ChannelT>::type a_t;
      m_ch[0] = ChannelT( ( a_t(other[0]) + a_t(other[1]) + a_t(other[2]) ) / 3 );
    }

    inline ChannelT& operator[](int i) { return m_ch[i]; }
    inline ChannelT const& operator[](int i) const { return m_ch[i]; }
    inline ChannelT& operator()(int i) { return m_ch[i]; }
    inline ChannelT const& operator()(int i) const { return m_ch[i]; }

    inline ChannelT& v() { return m_ch[0]; }
    inline ChannelT const& v() const { return m_ch[0]; }
  };

  /// \cond INTERNAL
  VW_DECLARE_PIXEL_TYPE(PixelGray,1);
  /// \endcond

  // Promote scalars for addition on the left
  template <class ScalarT, class ChannelT>
  typename boost::enable_if< IsScalar<ScalarT>, PixelGray<typename SumType<ScalarT,ChannelT>::type> >::type
  inline operator+( ScalarT scalar, PixelGray<ChannelT> pixel ) {
    return PixelGray<ScalarT>(scalar) + pixel;
  }

  // Promote scalars for addition on the right
  template <class ScalarT, class ChannelT>
  typename boost::enable_if< IsScalar<ScalarT>, PixelGray<typename SumType<ChannelT,ScalarT>::type> >::type
  inline operator+( PixelGray<ChannelT> pixel, ScalarT scalar ) {
    return pixel + PixelGray<ScalarT>(scalar);
  }

  // Promote scalars for subtraction on the left
  template <class ScalarT, class ChannelT>
  typename boost::enable_if< IsScalar<ScalarT>, PixelGray<typename DifferenceType<ScalarT,ChannelT>::type> >::type
  inline operator-( ScalarT scalar, PixelGray<ChannelT> pixel ) {
    return PixelGray<ScalarT>(scalar) - pixel;
  }

  // Promote scalars for subtraction on the right
  template <class ScalarT, class ChannelT>
  typename boost::enable_if< IsScalar<ScalarT>, PixelGray<typename DifferenceType<ChannelT,ScalarT>::type> >::type
  inline operator-( PixelGray<ChannelT> pixel, ScalarT scalar ) {
    return pixel - PixelGray<ScalarT>(scalar);
  }

  template <class ChannelT>
  std::ostream& operator<<( std::ostream& os, PixelGray<ChannelT> const& pix ) {
    return os << "Gray(" << _numeric(pix.v()) << ")";
  }


  // *******************************************************************
  // The PixelGrayA grayscale/alpha pixel type.
  // *******************************************************************

  /// Gray pixel type with alpha channel
  template <class ChannelT>
  class PixelGrayA : public PixelMathBase< PixelGrayA<ChannelT> >
  {
  private:
    ChannelT m_ch[2];
  public:
    PixelGrayA() { m_ch[0]=m_ch[1]=0; }
    PixelGrayA( ChannelT v ) { m_ch[0]=m_ch[1]=v; }
    PixelGrayA( ChannelT v, ChannelT a ) { m_ch[0]=v; m_ch[1]=a; }

    template <class OtherT> explicit PixelGrayA( PixelGray<OtherT> other ) {
      m_ch[0] = ChannelT(other[0]);
      m_ch[1] = ChannelRange<ChannelT>::max();
    }

    template <class OtherT> explicit PixelGrayA( PixelGrayA<OtherT> other ) {
      m_ch[0] = ChannelT(other[0]);
      m_ch[1] = ChannelT(other[1]);
    }

    template <class OtherT> explicit PixelGrayA( PixelRGB<OtherT> other ) {
      typedef typename AccumulatorType<ChannelT>::type a_t;
      m_ch[0] = ChannelT( ( a_t(other[0]) + a_t(other[1]) + a_t(other[2]) ) / 3 );
      m_ch[1] = ChannelRange<ChannelT>::max();
    }

    template <class OtherT> explicit PixelGrayA( PixelRGBA<OtherT> other ) {
      typedef typename AccumulatorType<ChannelT>::type a_t;
      m_ch[0] = ChannelT( ( a_t(other[0]) + a_t(other[1]) + a_t(other[2]) ) / 3 );
      m_ch[1] = ChannelT(other[3]);
    }

    inline ChannelT& operator[](int i) { return m_ch[i]; }
    inline ChannelT const& operator[](int i) const { return m_ch[i]; }
    inline ChannelT& operator()(int i) { return m_ch[i]; }
    inline ChannelT const& operator()(int i) const { return m_ch[i]; }

    inline ChannelT& v() { return m_ch[0]; }
    inline ChannelT const& v() const { return m_ch[0]; }
    inline ChannelT& a() { return m_ch[1]; }
    inline ChannelT const& a() const { return m_ch[1]; }
  };

  /// \cond INTERNAL
  VW_DECLARE_PIXEL_TYPE(PixelGrayA,2);
  /// \endcond

  template <class ChannelT>
  std::ostream& operator<<( std::ostream& os, PixelGrayA<ChannelT> const& pix ) {
    return os << "GrayA(" << _numeric(pix.v()) << "," << _numeric(pix.a()) << ")";
  }


  // *******************************************************************
  // The PixelRGB red/green/blue pixel type.
  // *******************************************************************

  /// RGB pixel type
  template <class ChannelT>
  class PixelRGB : public PixelMathBase< PixelRGB<ChannelT> >
  {
  private:
    ChannelT m_ch[3];
  public:
    PixelRGB() { m_ch[0]=m_ch[1]=m_ch[2]=0; }
    PixelRGB( ChannelT const& v ) { m_ch[0]=m_ch[1]=m_ch[2]=v; };
    PixelRGB( ChannelT const& r, ChannelT const& g, ChannelT const& b ) { m_ch[0]=r; m_ch[1]=g; m_ch[2]=b; }

    template <class OtherT> explicit PixelRGB( PixelGray<OtherT> const& other ) {
      m_ch[0] = m_ch[1] = m_ch[2] = ChannelT(other[0]);
    }
    
    template <class OtherT> explicit PixelRGB( PixelGrayA<OtherT> const& other ) {
      m_ch[0] = m_ch[1] = m_ch[2] = ChannelT(other[0]);
    }
    
    template <class OtherT> explicit PixelRGB( PixelRGB<OtherT> const& other ) {
      m_ch[0] = ChannelT(other[0]);
      m_ch[1] = ChannelT(other[1]);
      m_ch[2] = ChannelT(other[2]);
    }
    
    template <class OtherT> explicit PixelRGB( PixelRGBA<OtherT> const& other ) {
      m_ch[0] = ChannelT(other[0]);
      m_ch[1] = ChannelT(other[1]);
      m_ch[2] = ChannelT(other[2]);
    }
    
    template <class OtherT> explicit PixelRGB( PixelHSV<OtherT> const& other );

    inline ChannelT& operator[](int i) { return m_ch[i]; }
    inline ChannelT const& operator[](int i) const { return m_ch[i]; }
    inline ChannelT& operator()(int i) { return m_ch[i]; }
    inline ChannelT const& operator()(int i) const { return m_ch[i]; }

    inline ChannelT& r() { return m_ch[0]; }
    inline ChannelT const& r() const { return m_ch[0]; }
    inline ChannelT& g() { return m_ch[1]; }
    inline ChannelT const& g() const { return m_ch[1]; }
    inline ChannelT& b() { return m_ch[2]; }
    inline ChannelT const& b() const { return m_ch[2]; }
  };

  /// \cond INTERNAL
  VW_DECLARE_PIXEL_TYPE(PixelRGB,3);
  /// \endcond

  template <class ChannelT>
  std::ostream& operator<<( std::ostream& os, PixelRGB<ChannelT> const& pix ) {
    return os << "RGB(" << _numeric(pix.r()) << "," << _numeric(pix.g()) << "," << _numeric(pix.b()) << ")";
  }


  // *******************************************************************
  // The PixelRGB red/green/blue/alpha pixel type.
  // *******************************************************************

  /// RGB pixel type with alpha channel
  template <class ChannelT>
  class PixelRGBA : public PixelMathBase< PixelRGBA<ChannelT> >
  {
  private:
    ChannelT m_ch[4];
  public:
    PixelRGBA() { m_ch[0]=m_ch[1]=m_ch[2]=m_ch[3]=0; }
    PixelRGBA( ChannelT const& v ) { m_ch[0]=m_ch[1]=m_ch[2]=m_ch[3]=v; }
    PixelRGBA( ChannelT const& r, ChannelT const& g, ChannelT const& b, ChannelT const& a ) { m_ch[0]=r; m_ch[1]=g; m_ch[2]=b; m_ch[3]=a; }

    template <class OtherT> explicit PixelRGBA( PixelGray<OtherT> const& other ) {
      m_ch[0] = m_ch[1] = m_ch[2] = ChannelT(other[0]);
      m_ch[3]=ChannelRange<ChannelT>::max();
    }
    
    template <class OtherT> explicit PixelRGBA( PixelGrayA<OtherT> const& other ) {
      m_ch[0] = m_ch[1] = m_ch[2] = ChannelT(other[0]);
      m_ch[3]=ChannelT(other[1]);
    }
    
    template <class OtherT> explicit PixelRGBA( PixelRGB<OtherT> const& other ) {
      m_ch[0] = ChannelT(other[0]);
      m_ch[1] = ChannelT(other[1]);
      m_ch[2] = ChannelT(other[2]);
      m_ch[3] = ChannelRange<ChannelT>::max();
    }
    
    template <class OtherT> explicit PixelRGBA( PixelRGBA<OtherT> const& other ) {
      m_ch[0] = ChannelT(other[0]);
      m_ch[1] = ChannelT(other[1]);
      m_ch[2] = ChannelT(other[2]);
      m_ch[3] = ChannelT(other[3]);
    }

    inline ChannelT& operator[](int i) { return m_ch[i]; }
    inline ChannelT const& operator[](int i) const { return m_ch[i]; }
    inline ChannelT& operator()(int i) { return m_ch[i]; }
    inline ChannelT const& operator()(int i) const { return m_ch[i]; }

    inline ChannelT& r() { return m_ch[0]; }
    inline ChannelT const& r() const { return m_ch[0]; }
    inline ChannelT& g() { return m_ch[1]; }
    inline ChannelT const& g() const { return m_ch[1]; }
    inline ChannelT& b() { return m_ch[2]; }
    inline ChannelT const& b() const { return m_ch[2]; }
    inline ChannelT& a() { return m_ch[3]; }
    inline ChannelT const& a() const { return m_ch[3]; }
  };

  /// \cond INTERNAL
  VW_DECLARE_PIXEL_TYPE(PixelRGBA,4);
  /// \endcond

  template <class ChannelT>
  std::ostream& operator<<( std::ostream& os, PixelRGBA<ChannelT> const& pix ) {
    return os << "RGBA(" << _numeric(pix.r()) << "," << _numeric(pix.g()) << "," << _numeric(pix.b()) << "," << _numeric(pix.a()) << ")";
  }


  // *******************************************************************
  // The PixelHSV hue/saturation/value pixel type.
  // *******************************************************************

  /// HSV Pixel Type (Hue, Saturation, Value).
  template <class ChannelT>
  class PixelHSV : public PixelMathBase< PixelRGBA<ChannelT> >
  {
  private:
    ChannelT m_ch[3];
  public:
    PixelHSV() { m_ch[0]=m_ch[1]=m_ch[2]=0; }
    PixelHSV( ChannelT const& v ) { m_ch[0]=m_ch[1]=0; m_ch[2]=v; }
    PixelHSV( ChannelT const& h, ChannelT const& s, ChannelT const& v ) { m_ch[0]=h; m_ch[1]=s; m_ch[2]=v; }

    template <class OtherT> explicit PixelHSV( PixelHSV<OtherT> const& other ) {
      m_ch[0] = ChannelT(other[0]);
      m_ch[1] = ChannelT(other[1]);
      m_ch[2] = ChannelT(other[2]);
    }

    template <class OtherT> explicit PixelHSV( PixelRGB<OtherT> const& rgb );

    inline ChannelT& operator[](int i) { return m_ch[i]; }
    inline ChannelT const& operator[](int i) const { return m_ch[i]; }
    inline ChannelT& operator()(int i) { return m_ch[i]; }
    inline ChannelT const& operator()(int i) const { return m_ch[i]; }
    inline ChannelT& h() { return m_ch[0]; }
    inline ChannelT const& h() const { return m_ch[0]; }
    inline ChannelT& s() { return m_ch[1]; }
    inline ChannelT const& s() const { return m_ch[1]; }
    inline ChannelT& v() { return m_ch[2]; }
    inline ChannelT const& v() const { return m_ch[2]; }
  };

  // Here we override the usual PixelMathBase operators so that HSV pixel 
  // math agrees with RGB pixel math.  For addition this means converting 
  // to RGB and back.  For scalar multiplicaion it just means scaling the 
  // value channel.  If we end up doing this a lot we should probably 
  // factor this out into its own alternate CRTP base class.

  template <class ChannelT>
  inline PixelHSV<ChannelT> operator-( PixelHSV<ChannelT> const& pixel ) {
    return PixelHSV<ChannelT>( -PixelRGB<ChannelT>(pixel) );
  }

  template <class Channel1T, class Channel2T>
  PixelHSV<typename SumType<Channel1T,Channel2T>::type>
  inline operator+( PixelHSV<Channel1T> const& pixel1, PixelHSV<Channel2T> const& pixel2 ) {
    typedef typename SumType<Channel1T,Channel2T>::type ch_type;
    return PixelHSV<ch_type>( PixelRGB<ch_type>(pixel1) + PixelRGB<ch_type>(pixel2) );
  }

  template <class Channel1T, class Channel2T>
  PixelHSV<typename DifferenceType<Channel1T,Channel2T>::type>
  inline operator-( PixelHSV<Channel1T> const& pixel1, PixelHSV<Channel2T> const& pixel2 ) {
    typedef typename DifferenceType<Channel1T,Channel2T>::type ch_type;
    return PixelHSV<ch_type>( PixelRGB<ch_type>(pixel1) - PixelRGB<ch_type>(pixel2) );
  }

  template <class ChannelT, class ScalarT>
  typename boost::enable_if< IsScalar<ScalarT>, typename CompoundResult<ArgValProductFunctor<ScalarT>, PixelHSV<ChannelT> >::type >::type
  inline operator*( PixelHSV<ChannelT> const& pixel, ScalarT scalar ) {
    typename CompoundResult<ArgValProductFunctor<ScalarT>, PixelHSV<ChannelT> >::type result = pixel;
    result.v() *= scalar;
    return result;
  }

  template <class ScalarT, class ChannelT>
  typename boost::enable_if< IsScalar<ScalarT>, typename CompoundResult<ValArgProductFunctor<ScalarT>, PixelHSV<ChannelT> >::type >::type
  inline operator*( ScalarT scalar, PixelHSV<ChannelT> const& pixel ) {
    typename CompoundResult<ValArgProductFunctor<ScalarT>, PixelHSV<ChannelT> >::type result = pixel;
    result.v() *= scalar;
    return result;
  }

  template <class ChannelT, class ScalarT>
  typename boost::enable_if< IsScalar<ScalarT>, typename CompoundResult<ArgValQuotientFunctor<ScalarT>, PixelHSV<ChannelT> >::type >::type
  inline operator/( PixelHSV<ChannelT> const& pixel, ScalarT scalar ) {
    typename CompoundResult<ArgValQuotientFunctor<ScalarT>, PixelHSV<ChannelT> >::type result = pixel;
    pixel.v() /= scalar;
    return result;
  }

  /// \cond INTERNAL
  VW_DECLARE_PIXEL_TYPE(PixelHSV, 3);
  /// \endcond

  template <class ChannelT>
  std::ostream& operator<<( std::ostream& os, PixelHSV<ChannelT> const& pix ) {
    return os << "HSV(" << _numeric(pix.h()) << "," << _numeric(pix.s()) << "," << _numeric(pix.v()) << ")";
  }


  /// HSV->RGB conversion
  /// FIXME This really only supports floating-point pixel types
  template <class ChannelT>
  template <class OtherT>
  PixelRGB<ChannelT>::PixelRGB( PixelHSV<OtherT> const& hsv ) {
    int i = (int)floor(6*hsv.h());
    ChannelT f = 6*hsv.h() - i;
    if(!(i & 1)) f = 1 - f; // if i is even
    ChannelT m = hsv.v() * (1 - hsv.s());
    ChannelT n = hsv.v() * (1 - hsv.s() * f);
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
  /// FIXME This really only supports floating-point pixel types
  template <class ChannelT>
  template <class OtherT>
  PixelHSV<ChannelT>::PixelHSV( PixelRGB<OtherT> const& rgb ) {
    if( rgb.r()==rgb.g() && rgb.r()==rgb.b() ) { // H=1
      m_ch[0]=1;
      m_ch[1]=0;
      m_ch[2]=rgb.r();
    }
    else if( rgb.r()<=rgb.g() && rgb.r()<=rgb.b() ) { // H=[1/3,2/3)
      m_ch[2] = (rgb.g()>rgb.b()) ? rgb.g() : rgb.b();
      m_ch[0] = (3 - (rgb.g()-rgb.b())/(m_ch[2]-rgb.r()))/6.0;
      m_ch[1] = (m_ch[2]-rgb.r())/m_ch[2];
    }
    else if( rgb.g()<=rgb.r() && rgb.g()<rgb.b() ) { // H=[2/3,1)
      m_ch[2] = (rgb.r()>rgb.b()) ? rgb.r() : rgb.b();
      m_ch[0] = (5 - (rgb.b()-rgb.r())/(m_ch[2]-rgb.g()))/6.0;
      m_ch[1] = (m_ch[2]-rgb.g())/m_ch[2];
    }
    else { // H=[0,1/3) and H=1
      m_ch[2] = (rgb.r()>rgb.g()) ? rgb.r() : rgb.g();
      m_ch[0] = (1 - (rgb.r()-rgb.g())/(m_ch[2]-rgb.b()))/6.0;
      m_ch[1] = (m_ch[2]-rgb.b())/m_ch[2];
    }
  }


  // *******************************************************************
  // The Vector mathemaical vector pixel type.
  // *******************************************************************

  template <class ElemT, int SizeN>                                  
  struct CompoundChannelType<Vector<ElemT, SizeN> > {            
    typedef ElemT type;                                   
  };                                                         

  template <class ElemT, int SizeN>                                  
  struct CompoundNumChannels<Vector<ElemT, SizeN> > {           
    static const unsigned value = SizeN;                 
  };                                                         

  template <class OldChT, class NewChT, int SizeN>                      
  struct CompoundChannelCast<Vector<OldChT,SizeN>, NewChT> {       
    typedef Vector<NewChT, SizeN> type;                             
  };                                                         

  template <class OldChT, class NewChT, int SizeN>                      
  struct CompoundChannelCast<Vector<OldChT,SizeN>, const NewChT> {
    typedef const Vector<NewChT, SizeN> type;                       
  };

}

#endif // __VW_IMAGE_PIXELTYPES_H__
