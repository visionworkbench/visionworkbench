// __BEGIN_LICENSE__
//  Copyright (c) 2006-2026, United States Government as represented by the
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

/// \file AutoNormalize.h
///
/// Normalize overloads that auto-detect the input range via
/// min_max_channel_values(). These pull in Statistics.h, so they
/// are split out from Algorithms.h to keep that header lightweight.

#ifndef __VW_IMAGE_AUTO_NORMALIZE_H__
#define __VW_IMAGE_AUTO_NORMALIZE_H__

#include <vw/Image/Algorithms.h>
#include <vw/Image/Statistics.h>

namespace vw {

/// Renormalize the values in an image to fall within the range [low,high).
template <class ImageT>
UnaryPerPixelView<ImageT,
  UnaryCompoundFunctor<
    ChannelNormalizeFunctor<typename ImageT::pixel_type>,
    typename ImageT::pixel_type>>
normalize(ImageViewBase<ImageT> const& image,
          typename ImageChannelType<ImageT>::type low,
          typename ImageChannelType<ImageT>::type high) {
  typedef UnaryCompoundFunctor<
    ChannelNormalizeFunctor<typename ImageT::pixel_type>,
    typename ImageT::pixel_type> func_type;
  typename ImageChannelType<ImageT>::type old_min, old_max;
  min_max_channel_values(image, old_min, old_max);
  func_type func(ChannelNormalizeFunctor<typename ImageT::pixel_type>(
    old_min, old_max, low, high));
  return UnaryPerPixelView<ImageT, func_type>(image.impl(), func);
}

/// Renormalize the values in an image to fall within the range [0,high).
template <class ImageT>
UnaryPerPixelView<ImageT,
  UnaryCompoundFunctor<
    ChannelNormalizeFunctor<typename ImageT::pixel_type>,
    typename ImageT::pixel_type>>
normalize(ImageViewBase<ImageT> const& image,
          typename ImageChannelType<ImageT>::type high) {
  typedef UnaryCompoundFunctor<
    ChannelNormalizeFunctor<typename ImageT::pixel_type>,
    typename ImageT::pixel_type> func_type;
  typedef ChannelRange<typename ImageChannelType<ImageT>::type> range_type;
  typename ImageChannelType<ImageT>::type old_min, old_max;
  min_max_channel_values(image, old_min, old_max);
  func_type func(ChannelNormalizeFunctor<typename ImageT::pixel_type>(
    old_min, old_max, range_type::min(), high));
  return UnaryPerPixelView<ImageT, func_type>(image.impl(), func);
}

/// Renormalize the values in an image to fall within the range [min,max).
template <class ImageT>
UnaryPerPixelView<ImageT,
  UnaryCompoundFunctor<
    ChannelNormalizeFunctor<typename ImageT::pixel_type>,
    typename ImageT::pixel_type>>
normalize(ImageViewBase<ImageT> const& image) {
  typedef UnaryCompoundFunctor<
    ChannelNormalizeFunctor<typename ImageT::pixel_type>,
    typename ImageT::pixel_type> func_type;
  typedef ChannelRange<typename ImageChannelType<ImageT>::type> range_type;
  typename ImageChannelType<ImageT>::type old_min, old_max;
  min_max_channel_values(image, old_min, old_max);
  func_type func(ChannelNormalizeFunctor<typename ImageT::pixel_type>(
    old_min, old_max, range_type::min(), range_type::max()));
  return UnaryPerPixelView<ImageT, func_type>(image.impl(), func);
}

} // namespace vw

#endif // __VW_IMAGE_AUTO_NORMALIZE_H__
