// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file InterestTraits.h
/// 
/// Base type trait definitions for interest views.
/// 
#ifndef _INTEREST_POINT_INTEREST_TRAITS_H_
#define _INTEREST_POINT_INTEREST_TRAITS_H_

#include <vw/Image/ImageView.h>

namespace vw { 
namespace ip {

  /// This template function decides whether to use the memory-optimized
  /// default types or the potentially speed-optimized types associated with
  /// the interest view class. If the source type supports efficient access,
  /// we infer that speed is desired; otherwise we use the default types.
  ///
  /// InterestOperatorTraits can be partially specialized on the
  /// interest view to specify which processed views to fully
  /// rasterize. By default, it only rasterizes the interest image.
  template <class SrcT, class InterestT>
  struct InterestOperatorTraits {
    typedef ImageView<typename SrcT::pixel_type>                rasterize_type;
    typedef ImageView<typename SrcT::pixel_type>                gradient_type;
    typedef ImageView<typename SrcT::pixel_type>                mag_type;
    typedef ImageView<typename SrcT::pixel_type>                ori_type;
    typedef ImageView<typename SrcT::pixel_type>                interest_type;
  };

  // This should speed things up by delaying the rasterization of some
  // of these views, but instead it slows things down.  Why??
//   struct InterestOperatorTraits {
//     typedef ImageView<typename SrcT::pixel_type>                                                          rasterize_type;
//     typedef SeparableConvolutionView<SrcT, 
//                                      typename DefaultKernelT<typename SrcT::pixel_type>::type, 
//                                      ConstantEdgeExtension>                                               gradient_type;
//     typedef BinaryPerPixelView<gradient_type, gradient_type, vw::math::ArgArgHypotFunctor>                mag_type;
//     typedef BinaryPerPixelView<gradient_type, gradient_type, vw::math::ArgArgAtan2Functor>                ori_type;
//     typedef ImageView<typename SrcT::pixel_type>                                                          interest_type;
//   };
  
  /// Type(s) of peak in the interest image that indicate a feature.
  enum { IP_MAX, IP_MIN, IP_MINMAX };

  /// Peak type defaults to maxima.
  template <class InterestT>
  struct InterestPeakType { static const int peak_type = IP_MAX; };

} } //namespace vw::ip

#endif
