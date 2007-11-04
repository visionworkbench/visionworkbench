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

/// \file InterestTraits.h
/// 
/// Base type trait definitions for interest views.
/// 
#ifndef _INTEREST_POINT_INTEREST_TRAITS_H_
#define _INTEREST_POINT_INTEREST_TRAITS_H_

#include <vw/Image.h>

namespace vw { namespace ip {

  /// Template functions to identify the default processed view types.
  template <class SrcT>
  struct DefaultRasterizeT {
    typedef ImageView<typename SrcT::pixel_type> type;
  };

  template <class SrcT>
  struct DefaultGradientT {
    typedef SeparableConvolutionView<SrcT, typename DefaultKernelT<typename SrcT::pixel_type>::type, ConstantEdgeExtension> type;
  };
  
  template <class SrcT>
  struct DefaultMagT {
    typedef typename DefaultGradientT<SrcT>::type gradient_type;
    typedef BinaryPerPixelView<gradient_type, gradient_type, vw::math::ArgArgHypotFunctor> type;
  };

  template <class SrcT>
  struct DefaultOriT {
    typedef typename DefaultGradientT<SrcT>::type gradient_type;
    typedef BinaryPerPixelView<gradient_type, gradient_type, vw::math::ArgArgAtan2Functor> type;
  };

  template <class SrcT, class InterestT>
  struct DefaultInterestT {
    typedef typename InterestT::template ViewType<SrcT> type;
  };

  // TODO: how should we handle ImageResourceViews?
  // TODO: should enforce defaults only when SrcT is ImageViewRef instead of checking
  //       IsMultiplyAccessible (also solves above)

  /// By default, we do not rasterize any of the processed views to conserve
  /// memory and computation.
  template <class SrcT, class InterestT>
  struct InterestDefaultTraits {
    typedef typename DefaultRasterizeT<SrcT>::type          rasterize_type;
    typedef typename DefaultGradientT<SrcT>::type           gradient_type;
    typedef typename DefaultMagT<SrcT>::type                mag_type;
    typedef typename DefaultOriT<SrcT>::type                ori_type;
    typedef typename DefaultInterestT<SrcT,InterestT>::type interest_type;
  };

  /// InterestOperatorTraits can be partially specialized on the
  /// interest view to specify which processed views to fully
  /// rasterize. By default, it only rasterizes the interest image.
  template <class SrcT, class InterestT>
  struct InterestOperatorTraits {
    typedef typename DefaultRasterizeT<SrcT>::type          rasterize_type;
    typedef typename DefaultGradientT<SrcT>::type           gradient_type;
    typedef typename DefaultMagT<SrcT>::type                mag_type;
    typedef typename DefaultOriT<SrcT>::type                ori_type;
    typedef rasterize_type                                  interest_type;
  };

  /// This template function decides whether to use the memory-optimized
  /// default types or the potentially speed-optimized types associated with
  /// the interest view class. If the source type supports efficient access,
  /// we infer that speed is desired; otherwise we use the default types.
  template <class SrcT, class InterestT>
  struct InterestTraitsHelper {
    //typedef typename boost::mpl::if_< IsMultiplyAccessible<SrcT>, InterestOperatorTraits<SrcT, InterestT>, InterestDefaultTraits<SrcT, InterestT> >::type type;
    typedef InterestOperatorTraits<SrcT,InterestT> type;
  };

  template <class PixelT, class InterestT>
  struct InterestTraitsHelper <ImageViewRef<PixelT>, InterestT> {
    typedef InterestDefaultTraits<ImageViewRef<PixelT>,InterestT> type;
  };

  /// Inherits the inferred appropriate type traits.
  template <class SrcT, class InterestT>
  struct InterestTraits : public InterestTraitsHelper<SrcT,InterestT>::type {};

  /// Type(s) of peak in the interest image that indicate a feature.
  enum { IP_MAX, IP_MIN, IP_MINMAX };

  /// Peak type defaults to maxima.
  template <class InterestT>
  struct InterestPeakType { static const int peak_type = IP_MAX; };

} } //namespace vw::ip

#endif
