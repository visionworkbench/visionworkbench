// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
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
    typedef ImageView<typename SrcT::pixel_type>  rasterize_type;
    typedef ImageView<typename SrcT::pixel_type>  gradient_type;
    typedef ImageView<typename SrcT::pixel_type>  mag_type;
    typedef ImageView<typename SrcT::pixel_type>  ori_type;
    typedef ImageView<typename SrcT::pixel_type>  interest_type;
    typedef ImageView<typename SrcT::pixel_type>  integral_type;
  };

  /// Type(s) of peak in the interest image that indicate a feature.
  enum { IP_MAX, IP_MIN, IP_MINMAX };

  /// Peak type defaults to maxima.
  template <class InterestT>
  struct InterestPeakType { static const int peak_type = IP_MAX; };

} } //namespace vw::ip

#endif
