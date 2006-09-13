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

/// \file ImageOperators.h
/// 
/// Arithmetic operators for operating on images.
/// 
/// This header defines the standard set of arithmetic operators for 
/// images.  The following are currently supported:
/// 
///  - <TT>-</TT> image
///  - image <TT>+</TT> image
///  - image <TT>+=</TT> image
///  - image <TT>+ </TT>scalar
///  - image <TT>+= </TT>scalar
///  - scalar <TT>+</TT> image
///  - image <TT>-</TT> image
///  - image <TT>-=</TT> image
///  - image <TT>- </TT>scalar
///  - image <TT>-= </TT>scalar
///  - scalar <TT>-</TT> image
///  - image <TT>*</TT> image
///  - image <TT>*=</TT> image
///  - image <TT>* </TT>scalar
///  - image <TT>*= </TT>scalar
///  - scalar <TT>*</TT> image
///  - image <TT>/</TT> image
///  - image <TT>/=</TT> image
///  - image <TT>/ </TT>scalar
///  - image <TT>/= </TT>scalar
///  - scalar <TT>/</TT> image
///
/// Each operator returns either an UnaryPerPixelImageView or a 
/// BinaryPerPixelView, as appropriate, so they can be combined 
/// efficiently without introducing intermediate image buffers.  
/// In order for the image operators to work, the corresponding 
/// operators for the underlying pixel types must exist.
///
/// We currently don't support the comparison operators here 
/// 
#ifndef __VW_IMAGE__OPERATORS_H__
#define __VW_IMAGE__OPERATORS_H__

#include <boost/utility/enable_if.hpp>

#include <vw/Core/Functors.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/PerPixelViews.h>

namespace vw {

  // ********************************************************************
  // ImageView Arithmetic Operations
  //
  // The image view operators themselves return PerPixelImageViews 
  // using the above PixelType aruthmetic operation functors.
  //
  // ********************************************************************

  /// Negation of an image.
  template <class ImageT>
  UnaryPerPixelView<ImageT,ArgNegationFunctor>
  inline operator-( ImageViewBase<ImageT> const& image ) {
    return UnaryPerPixelView<ImageT,ArgNegationFunctor>( image.impl() );
  }


  /// Sum of two images.
  template <class Image1T, class Image2T>
  BinaryPerPixelView<Image1T,Image2T,ArgArgSumFunctor>
  inline operator+( ImageViewBase<Image1T> const& image1, ImageViewBase<Image2T> const& image2 ) {
    return BinaryPerPixelView<Image1T,Image2T,ArgArgSumFunctor>( image1.impl(), image2.impl() );
  }

  /// Sum of an image and a scalar.
  template <class ImageT, class ScalarT>
  typename boost::disable_if< IsImageView<ScalarT>,
                              UnaryPerPixelView<ImageT,ArgValSumFunctor<ScalarT> >
                              >::type
  inline operator+( ImageViewBase<ImageT> const& image, ScalarT const& scalar ) {
    return UnaryPerPixelView<ImageT,ArgValSumFunctor<ScalarT> >( image.impl(), ArgValSumFunctor<ScalarT>(scalar) );
  }

  /// Sum of a scalar and an image.
  template <class ImageT, class ScalarT>
  typename boost::disable_if< IsImageView<ScalarT>,
                              UnaryPerPixelView<ImageT,ValArgSumFunctor<ScalarT> >
                              >::type
  inline operator+( ScalarT const& scalar, ImageViewBase<ImageT> const& image ) {
    return UnaryPerPixelView<ImageT,ValArgSumFunctor<ScalarT> >( image.impl(), ValArgSumFunctor<ScalarT>(scalar) );
  }

  /// Sum-assignment of an image.
  template <class ImageT, class ArgT>
  inline ImageT const& operator+=( ImageViewBase<ImageT> const& image, ArgT const& arg ) {
    return (image.impl() = image + arg);
  }


  /// Difference of two images.
  template <class Image1T, class Image2T>
  BinaryPerPixelView<Image1T,Image2T,ArgArgDifferenceFunctor>
  inline operator-( ImageViewBase<Image1T> const& image1, ImageViewBase<Image2T> const& image2 ) {
    return BinaryPerPixelView<Image1T,Image2T,ArgArgDifferenceFunctor>( image1.impl(), image2.impl() );
  }

  /// Difference of an image and a scalar.
  template <class ImageT, class ScalarT>
  typename boost::disable_if< IsImageView<ScalarT>,
                              UnaryPerPixelView<ImageT,ArgValDifferenceFunctor<ScalarT> >
                              >::type
  inline operator-( ImageViewBase<ImageT> const& image, ScalarT const& scalar ) {
    return UnaryPerPixelView<ImageT,ArgValDifferenceFunctor<ScalarT> >( image.impl(), ArgValDifferenceFunctor<ScalarT>(scalar) );
  }

  /// Difference of a scalar and an image.
  template <class ImageT, class ScalarT>
  typename boost::disable_if< IsImageView<ScalarT>,
                              UnaryPerPixelView<ImageT,ValArgDifferenceFunctor<ScalarT> >
                              >::type
  inline operator-( ScalarT const& scalar, ImageViewBase<ImageT> const& image ) {
    return UnaryPerPixelView<ImageT,ValArgDifferenceFunctor<ScalarT> >( image.impl(), ValArgDifferenceFunctor<ScalarT>(scalar) );
  }

  /// Difference-assignment of an image.
  template <class ImageT, class ArgT>
  inline ImageT const& operator-=( ImageViewBase<ImageT> const& image, ArgT const& arg ) {
    return (image.impl() = image - arg);
  }


  /// Product of two images.
  template <class Image1T, class Image2T>
  BinaryPerPixelView<Image1T,Image2T,ArgArgProductFunctor>
  inline operator*( ImageViewBase<Image1T> const& image1, ImageViewBase<Image2T> const& image2 ) {
    return BinaryPerPixelView<Image1T,Image2T,ArgArgProductFunctor>( image1.impl(), image2.impl() );
  }
  
  /// Product of an image and a scalar.
  template <class ImageT, class ScalarT>
  typename boost::disable_if< IsImageView<ScalarT>,
                             UnaryPerPixelView<ImageT,ArgValProductFunctor<ScalarT> >
                             >::type
  inline operator*( ImageViewBase<ImageT> const& image, ScalarT const& scalar ) {
    return UnaryPerPixelView<ImageT,ArgValProductFunctor<ScalarT> >( image.impl(), ArgValProductFunctor<ScalarT>(scalar) );
  }

  /// Product of a scalar and an image.
  template <class ImageT, class ScalarT>
  typename boost::disable_if< IsImageView<ScalarT>,
                              UnaryPerPixelView<ImageT,ValArgProductFunctor<ScalarT> >
                              >::type
  inline operator*( ScalarT const& scalar, ImageViewBase<ImageT> const& image ) {
    return UnaryPerPixelView<ImageT,ValArgProductFunctor<ScalarT> >( image.impl(), ValArgProductFunctor<ScalarT>(scalar) );
  }

  /// Product-assignment of an image.
  template <class ImageT, class ArgT>
  inline ImageT const& operator*=( ImageViewBase<ImageT> const& image, ArgT const& arg ) {
    return (image.impl() = image * arg);
  }


  /// Quotient of two images.
  template <class Image1T, class Image2T>
  BinaryPerPixelView<Image1T,Image2T,ArgArgQuotientFunctor>
  inline operator/( ImageViewBase<Image1T> const& image1, ImageViewBase<Image2T> const& image2 ) {
    return BinaryPerPixelView<Image1T,Image2T,ArgArgQuotientFunctor>( image1.impl(), image2.impl() );
  }

  /// Quotient of an image and a scalar.
  template <class ImageT, class ScalarT>
  typename boost::disable_if< IsImageView<ScalarT>,
                              UnaryPerPixelView<ImageT,ArgValQuotientFunctor<ScalarT> >
                              >::type
  inline operator/( ImageViewBase<ImageT> const& image, ScalarT const& scalar ) {
    return UnaryPerPixelView<ImageT,ArgValQuotientFunctor<ScalarT> >( image.impl(), ArgValQuotientFunctor<ScalarT>(scalar) );
  }
  
  /// Quotient of a scalar and an image.
  template <class ImageT, class ScalarT>
  typename boost::disable_if< IsImageView<ScalarT>,
                              UnaryPerPixelView<ImageT,ValArgQuotientFunctor<ScalarT> >
                              >::type
  inline operator/( ScalarT const& scalar, ImageViewBase<ImageT> const& image ) {
    return UnaryPerPixelView<ImageT,ValArgQuotientFunctor<ScalarT> >( image.impl(), ValArgQuotientFunctor<ScalarT>(scalar) );
  }

  /// Quotient-assignment of an image.
  template <class ImageT, class ArgT>
  inline ImageT const& operator/=( ImageViewBase<ImageT> const& image, ArgT const& arg ) {
    return (image.impl() = image / arg);
  }

} // namespace vw

#endif // __VW_IMAGE__OPERATORS_H__
