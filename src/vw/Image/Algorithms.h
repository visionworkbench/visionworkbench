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

/// \file Algorithms.h
///
/// Lazy per-pixel algorithms: clamp, normalize, threshold.
/// Surface analysis views (MeanFillTransparent, ComputeNormals,
/// DotProd, TwoThresholdFill) are in ImageSurface.h.

#ifndef __VW_IMAGE_ALGORITHMS_H__
#define __VW_IMAGE_ALGORITHMS_H__

#include <vw/Image/AlgorithmFunctions.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Image/PerPixelAccessorViews.h>

// Backward compatibility
#include <vw/Image/ImageSurface.h>

namespace vw {

// clamp()

template <class PixelT>
class ChannelClampFunctor: public UnaryReturnSameType {
  typedef typename CompoundChannelType<PixelT>::type channel_type;
  channel_type m_low, m_high;
public:
  ChannelClampFunctor(channel_type low, channel_type high):
    m_low(low), m_high(high) {}

  channel_type operator()(channel_type value) const {
    if      (value > m_high) return m_high;
    else if (value < m_low)  return m_low;
    else                     return value;
  }
};

/// Clamp the values in an image to fall within the range [low,high].
template <class ImageT, class LowT, class HighT>
UnaryPerPixelView<ImageT,
  UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type>,
                       typename ImageT::pixel_type>>
clamp(ImageViewBase<ImageT> const& image, LowT low, HighT high) {
  typedef UnaryCompoundFunctor<
    ChannelClampFunctor<typename ImageT::pixel_type>,
    typename ImageT::pixel_type> func_type;
  func_type func(ChannelClampFunctor<typename ImageT::pixel_type>(low, high));
  return UnaryPerPixelView<ImageT, func_type>(image.impl(), func);
}

/// Clamp the values in an image to fall within the range [0,high].
template <class ImageT, class HighT>
UnaryPerPixelView<ImageT,
  UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type>,
                       typename ImageT::pixel_type>>
clamp(ImageViewBase<ImageT> const& image, HighT high) {
  typedef UnaryCompoundFunctor<
    ChannelClampFunctor<typename ImageT::pixel_type>,
    typename ImageT::pixel_type> func_type;
  typedef ChannelRange<
    typename CompoundChannelType<typename ImageT::pixel_type>::type>
    range_type;
  typename CompoundChannelType<typename ImageT::pixel_type>::type
    min_val = range_type::min();
  func_type func(
    ChannelClampFunctor<typename ImageT::pixel_type>(min_val, high));
  return UnaryPerPixelView<ImageT, func_type>(image.impl(), func);
}

/// Clamp the values in an image to fall within the range [min,max].
template <class ImageT>
UnaryPerPixelView<ImageT,
  UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type>,
                       typename ImageT::pixel_type>>
clamp(ImageViewBase<ImageT> const& image) {
  typedef UnaryCompoundFunctor<
    ChannelClampFunctor<typename ImageT::pixel_type>,
    typename ImageT::pixel_type> func_type;
  typedef ChannelRange<
    typename CompoundChannelType<typename ImageT::pixel_type>::type>
    range_type;
  typename CompoundChannelType<typename ImageT::pixel_type>::type
    min_val = range_type::min();
  typename CompoundChannelType<typename ImageT::pixel_type>::type
    max_val = range_type::max();
  func_type func(
    ChannelClampFunctor<typename ImageT::pixel_type>(min_val, max_val));
  return UnaryPerPixelView<ImageT, func_type>(image.impl(), func);
}

// normalize()

template <class PixelT>
class ChannelNormalizeFunctor: public UnaryReturnSameType {
  typedef typename CompoundChannelType<PixelT>::type channel_type;
  channel_type m_old_min, m_new_min;
  double m_old_to_new_ratio;
public:
  ChannelNormalizeFunctor(channel_type old_min, channel_type old_max,
                          channel_type new_min, channel_type new_max):
    m_old_min(old_min), m_new_min(new_min) {
    if (old_max == old_min)
      m_old_to_new_ratio = 0.0;
    else
      m_old_to_new_ratio =
        (new_max - new_min) / (double)(old_max - old_min);
  }

  template <class ChannelT>
  ChannelT operator()(ChannelT value) const {
    return (ChannelT)((value - m_old_min) * m_old_to_new_ratio +
                      m_new_min);
  }
};

template <class PixelT>
class ChannelNormalizeRetainAlphaFunctor: public UnaryReturnSameType {
  typedef typename CompoundChannelType<PixelT>::type channel_type;
  typedef typename PixelWithoutAlpha<PixelT>::type non_alpha_type;
  typedef ChannelNormalizeFunctor<non_alpha_type> norm_func_type;
  UnaryCompoundFunctor<norm_func_type, non_alpha_type> m_compound_func;
public:
  ChannelNormalizeRetainAlphaFunctor(channel_type old_min,
                                     channel_type old_max,
                                     channel_type new_min,
                                     channel_type new_max):
    m_compound_func(
      norm_func_type(old_min, old_max, new_min, new_max)) {}

  PixelT operator()(PixelT value) const {
    if (is_transparent(value))
      return value;
    else {
      PixelT result;
      non_alpha_channels(result) =
        m_compound_func(non_alpha_channels(value));
      alpha_channel(result) = alpha_channel(value);
      return result;
    }
  }
};

/// Renormalize the values in an image to fall within the range
/// [low,high), but leave the values in the alpha channel untouched.
template <class ImageT>
UnaryPerPixelView<ImageT,
  ChannelNormalizeRetainAlphaFunctor<typename ImageT::pixel_type>>
normalize_retain_alpha(ImageViewBase<ImageT> const& image,
                       typename ImageChannelType<ImageT>::type old_low,
                       typename ImageChannelType<ImageT>::type old_high,
                       typename ImageChannelType<ImageT>::type new_low,
                       typename ImageChannelType<ImageT>::type new_high) {
  typedef ChannelNormalizeRetainAlphaFunctor<
    typename ImageT::pixel_type> func_type;
  func_type func(old_low, old_high, new_low, new_high);
  return UnaryPerPixelView<ImageT, func_type>(image.impl(), func);
}

/// Renormalize the values in an image to fall within the range [low,high).
template <class ImageT>
UnaryPerPixelView<ImageT,
  UnaryCompoundFunctor<
    ChannelNormalizeFunctor<typename ImageT::pixel_type>,
    typename ImageT::pixel_type>>
normalize(ImageViewBase<ImageT> const& image,
          typename ImageChannelType<ImageT>::type old_low,
          typename ImageChannelType<ImageT>::type old_high,
          typename ImageChannelType<ImageT>::type new_low,
          typename ImageChannelType<ImageT>::type new_high) {
  typedef UnaryCompoundFunctor<
    ChannelNormalizeFunctor<typename ImageT::pixel_type>,
    typename ImageT::pixel_type> func_type;
  func_type func(ChannelNormalizeFunctor<typename ImageT::pixel_type>(
    old_low, old_high, new_low, new_high));
  return UnaryPerPixelView<ImageT, func_type>(image.impl(), func);
}

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

// threshold()

template <class PixelT>
class ChannelThresholdFunctor {
  typedef typename CompoundChannelType<PixelT>::type channel_type;
  channel_type m_thresh, m_low, m_high;
public:
  ChannelThresholdFunctor(channel_type thresh, channel_type low,
                          channel_type high):
    m_thresh(thresh), m_low(low), m_high(high) {}

  template <class Args> struct result {
    typedef channel_type type;
  };

  inline channel_type operator()(channel_type const& val) const {
    return (val > m_thresh) ? m_high : m_low;
  }
};

/// Threshold the values in an image, generating a two-valued output
/// image with values low and high.
template <class ImageT, class ThreshT, class LowT, class HighT>
UnaryPerPixelView<ImageT,
  UnaryCompoundFunctor<
    ChannelThresholdFunctor<typename ImageT::pixel_type>,
    typename ImageT::pixel_type>>
threshold(ImageViewBase<ImageT> const& image,
          ThreshT thresh, LowT low, HighT high) {
  typedef UnaryCompoundFunctor<
    ChannelThresholdFunctor<typename ImageT::pixel_type>,
    typename ImageT::pixel_type> func_type;
  func_type func(
    ChannelThresholdFunctor<typename ImageT::pixel_type>(thresh, low, high));
  return UnaryPerPixelView<ImageT, func_type>(image.impl(), func);
}

/// Threshold the values in an image with values 0 and high.
template <class ImageT, class ThreshT, class HighT>
UnaryPerPixelView<ImageT,
  UnaryCompoundFunctor<
    ChannelThresholdFunctor<typename ImageT::pixel_type>,
    typename ImageT::pixel_type>>
threshold(ImageViewBase<ImageT> const& image,
          ThreshT thresh, HighT high) {
  typedef UnaryCompoundFunctor<
    ChannelThresholdFunctor<typename ImageT::pixel_type>,
    typename ImageT::pixel_type> func_type;
  typedef ChannelRange<typename ImageChannelType<ImageT>::type> range_type;
  func_type func(ChannelThresholdFunctor<typename ImageT::pixel_type>(
    thresh, range_type::min(), high));
  return UnaryPerPixelView<ImageT, func_type>(image.impl(), func);
}

/// Threshold the values in an image against a threshold.
template <class ImageT, class ThreshT>
UnaryPerPixelView<ImageT,
  UnaryCompoundFunctor<
    ChannelThresholdFunctor<typename ImageT::pixel_type>,
    typename ImageT::pixel_type>>
threshold(ImageViewBase<ImageT> const& image, ThreshT thresh) {
  typedef UnaryCompoundFunctor<
    ChannelThresholdFunctor<typename ImageT::pixel_type>,
    typename ImageT::pixel_type> func_type;
  typedef ChannelRange<typename ImageChannelType<ImageT>::type> range_type;
  func_type func(ChannelThresholdFunctor<typename ImageT::pixel_type>(
    thresh, range_type::min(), range_type::max()));
  return UnaryPerPixelView<ImageT, func_type>(image.impl(), func);
}

/// Threshold the values in an image against zero.
template <class ImageT>
UnaryPerPixelView<ImageT,
  UnaryCompoundFunctor<
    ChannelThresholdFunctor<typename ImageT::pixel_type>,
    typename ImageT::pixel_type>>
threshold(ImageViewBase<ImageT> const& image) {
  typedef UnaryCompoundFunctor<
    ChannelThresholdFunctor<typename ImageT::pixel_type>,
    typename ImageT::pixel_type> func_type;
  typedef ChannelRange<typename ImageChannelType<ImageT>::type> range_type;
  func_type func(ChannelThresholdFunctor<typename ImageT::pixel_type>(
    0, range_type::min(), range_type::max()));
  return UnaryPerPixelView<ImageT, func_type>(image.impl(), func);
}

} // namespace vw

#endif // __VW_IMAGE_ALGORITHMS_H__
