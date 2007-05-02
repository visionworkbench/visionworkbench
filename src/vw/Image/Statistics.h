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

/// \file Statistics.h
/// 
/// These functions compute the min, max, minmax, and
/// mean value of an image on a per-channel basis.  
///

// FIXME: How do we deal with alpha channels here?
#ifndef __VW_IMAGE_STATISTICS_H__
#define __VW_IMAGE_STATISTICS_H__

#include <boost/type_traits.hpp>

#include <vw/Image/ImageView.h>
#include <vw/Image/Manipulation.h>

namespace vw {

  /// Compute the minimum value stored in all of the channels of all of the planes of the images.
  template <class ViewT>
  typename CompoundChannelType<typename ViewT::pixel_type>::type
  min_channel_value( const ImageViewBase<ViewT>& view_ ) {
    const ViewT& view = view_.impl();
 
    typedef typename ViewT::pixel_accessor pixel_accessor;
    typedef typename ViewT::pixel_type pixel_type;
    typedef typename CompoundChannelType<typename ViewT::pixel_type>::type channel_type;

    channel_type min = ChannelLimits<channel_type>::max;
    bool valid = false;

    int num_channels = view.channels();
    if (PixelHasAlpha<pixel_type>::value)
      num_channels = view.channels() - 1;

    pixel_accessor plane_acc = view.origin();
    for (int32 p = 0; p < view.planes(); p++, plane_acc.next_plane()) { 
      pixel_accessor col_acc = plane_acc;
      for (int32 i = 0; i < view.cols(); i++, col_acc.next_col()) {
        pixel_accessor row_acc = col_acc;
        for (int32 j = 0; j < view.rows(); j++, row_acc.next_row()) {
          typename ViewT::result_type pix = *row_acc;
          for (int32 channel = 0; channel < num_channels; channel++) {
            if (!is_transparent(pix)) {
              valid = true;
              channel_type channel_value = compound_select_channel<channel_type>(pix,channel);
              if( channel_value < min ) min = channel_value;
            }
          }
        }
      }  
    }
    if (!valid) 
      vw_throw(ArgumentErr() << "min_channel_value(): the image contained zero valid pixels.");
    return min;
  }


  /// Compute the maximum value stored in all of the channels of all of the planes of the images.
  template <class ViewT>
  typename CompoundChannelType<typename ViewT::pixel_type>::type
  max_channel_value( const ImageViewBase<ViewT>& view_ ) {
    const ViewT& view = view_.impl();

    typedef typename ViewT::pixel_accessor pixel_accessor;
    typedef typename ViewT::pixel_type pixel_type;
    typedef typename CompoundChannelType<typename ViewT::pixel_type>::type channel_type;

    channel_type max = ChannelLimits<channel_type>::min;
    bool valid = false;

    int num_channels = view.channels();
    if (PixelHasAlpha<pixel_type>::value)
      num_channels = view.channels() - 1;

    pixel_accessor plane_acc = view.origin();
    for (int32 p = 0; p < view.planes(); p++, plane_acc.next_plane()) { 
      pixel_accessor col_acc = plane_acc;
      for (int32 i = 0; i < view.cols(); i++, col_acc.next_col()) {
        pixel_accessor row_acc = col_acc;
        for (int32 j = 0; j < view.rows(); j++, row_acc.next_row()) {
          typename ViewT::result_type pix = *row_acc;
          for (int32 channel = 0; channel < num_channels; channel++) {
            if (!is_transparent(pix)) {
              valid = true;
              channel_type channel_value = compound_select_channel<channel_type>(pix,channel);
              if( channel_value > max ) max = channel_value;
            }
          }
        }
      } 
    }
    if (!valid) 
      vw_throw(ArgumentErr() << "max_channel_value(): the image contained zero valid pixels.");
    return max;
  }


  /// Simultaneously compute the min and max value in all of the
  /// channels of all of the planes of the image.
  template <class ViewT>
  void min_max_channel_values( const ImageViewBase<ViewT> &view_, 
                               typename CompoundChannelType<typename ViewT::pixel_type>::type &min, 
                               typename CompoundChannelType<typename ViewT::pixel_type>::type &max ) {
    const ViewT& view = view_.impl();

    typedef typename ViewT::pixel_accessor pixel_accessor;
    typedef typename ViewT::pixel_type pixel_type;
    typedef typename CompoundChannelType<typename ViewT::pixel_type>::type channel_type;

    min = ChannelLimits<channel_type>::max;
    max = ChannelLimits<channel_type>::min;
    
    int num_channels = view.channels();
    if (PixelHasAlpha<pixel_type>::value)
      num_channels = view.channels() - 1;

    pixel_accessor plane_acc = view.origin();
    for (int32 p = 0; p < view.planes(); p++, plane_acc.next_plane()) { 
      pixel_accessor col_acc = plane_acc;
      for (int32 i = 0; i < view.cols(); i++, col_acc.next_col()) {
        pixel_accessor row_acc = col_acc;
        for (int32 j = 0; j < view.rows(); j++, row_acc.next_row()) {
          typename ViewT::result_type pix = *row_acc;
          for (int32 channel = 0; channel < num_channels; channel++) {
            if (!is_transparent(pix)) {
              channel_type channel_value = compound_select_channel<channel_type>(pix,channel);
              if( channel_value < min ) min = channel_value;
              if( channel_value > max ) max = channel_value;
            }
          }
        }
      }  
    }

    // Before we allow this function to return, we check to make sure
    // that there was at least one valid pixel in the image.
    if (min > max) 
      vw_throw(ArgumentErr() << "min_max_channel_values(): the image contained zero valid pixels.");
  }


  /// Compute the mean value stored in all of the channels of all of the planes of the image.
  template <class ViewT>
  typename CompoundChannelType<typename ViewT::pixel_type>::type
  mean_channel_value( const ImageViewBase<ViewT> &view_ ) {
    const ViewT& view = view_.impl();

    typedef typename ViewT::pixel_accessor pixel_accessor;
    typedef typename CompoundChannelType<typename ViewT::pixel_type>::type channel_type;

    double accum = 0;
    bool valid = false;

    int num_channels = view.channels();
    if (PixelHasAlpha<typename ViewT::pixel_type>::value)
      num_channels = view.channels() - 1;

    pixel_accessor plane_acc = view.origin();
    for (int32 p = 0; p < view.planes(); p++, plane_acc.next_plane()) { 
      pixel_accessor col_acc = plane_acc;
      for (int32 i = 0; i < view.cols(); i++, col_acc.next_col()) {
        pixel_accessor row_acc = col_acc;
        for (int32 j = 0; j < view.rows(); j++, row_acc.next_row()) {
          typename ViewT::result_type pix = *row_acc;
          for (int32 channel = 0; channel < num_channels; channel++) {
            if (!is_transparent(pix)) {
              valid = true;
              channel_type channel_value = compound_select_channel<channel_type>(pix,channel);
              accum += channel_value;
            }
          }
        }
      }  
    }
    if (!valid) 
      vw_throw(ArgumentErr() << "mean_channel_value(): the image contained zero valid pixels.");
    return channel_type(accum / (view.planes() * view.rows() * view.cols() * view.channels()));
  }

  /// Compute the standard deviation of the values stored in all of the channels of all of the planes of the image.
  template <class ViewT>
  typename CompoundChannelType<typename ViewT::pixel_type>::type
  stddev_channel_value( const ImageViewBase<ViewT> &view_ ) {
    const ViewT& view = view_.impl();

    typedef typename ViewT::pixel_accessor pixel_accessor;
    typedef typename CompoundChannelType<typename ViewT::pixel_type>::type channel_type;

    double accum = 0;
    double mean = mean_channel_value(view);
    bool valid = false;

    int num_channels = view.channels();
    if (PixelHasAlpha<typename ViewT::pixel_type>::value)
      num_channels = view.channels() - 1;

    pixel_accessor plane_acc = view.origin();
    for (int32 p = 0; p < view.planes(); p++, plane_acc.next_plane()) { 
      pixel_accessor col_acc = plane_acc;
      for (int32 i = 0; i < view.cols(); i++, col_acc.next_col()) {
        pixel_accessor row_acc = col_acc;
        for (int32 j = 0; j < view.rows(); j++, row_acc.next_row()) {
          typename ViewT::result_type pix = *row_acc;
          for (int32 channel = 0; channel < num_channels; channel++) {
            if (!is_transparent(pix)) {
              valid = true;
              channel_type channel_value = compound_select_channel<channel_type>(pix,channel);
              accum += pow(channel_value - mean, 2);
            }
          }
        }
      }  
    }
    if (!valid) 
      vw_throw(ArgumentErr() << "stddev_channel_value(): the image contained zero valid pixels.");
    return channel_type(sqrt(accum) / ((view.planes() * view.rows() * view.cols() * view.channels()) - 1));
  }

  /// Compute the sum of the pixels of the channels of all of the planes of the image.
  template <class ViewT>
  typename AccumulatorType<typename PixelChannelType<typename ViewT::pixel_type>::type>::type
  sum_of_channel_values( const ImageViewBase<ViewT> &view_ ) {
    typedef typename AccumulatorType<typename CompoundChannelType<typename ViewT::pixel_type>::type>::type accum_type;
    const ViewT& view = view_.impl();

    typedef typename ViewT::pixel_accessor pixel_accessor;
    typedef typename CompoundChannelType<typename ViewT::pixel_type>::type channel_type;

    accum_type accum = 0;
    bool valid = false;

    int num_channels = view.channels();
    if (PixelHasAlpha<typename ViewT::pixel_type>::value)
      num_channels = view.channels() - 1;

    pixel_accessor plane_acc = view.origin();
    for (int32 p = 0; p < view.planes(); p++, plane_acc.next_plane()) { 
      pixel_accessor col_acc = plane_acc;
      for (int32 i = 0; i < view.cols(); i++, col_acc.next_col()) {
        pixel_accessor row_acc = col_acc;
        for (int32 j = 0; j < view.rows(); j++, row_acc.next_row()) {
          typename ViewT::result_type pix = *row_acc;
          for (int32 channel = 0; channel < num_channels; channel++) {
            if (!is_transparent(pix)) {
              valid = true;
              channel_type channel_value = compound_select_channel<channel_type>(pix,channel);
              accum += channel_value;
            }
          }
        }
      }  
    }
    if (!valid) 
      vw_throw(ArgumentErr() << "sum_of_channel_values(): the image contained zero valid pixels.");
    return accum;
  }


  /// Compute the median pixel value of an image.  This routine
  /// computes the exact median by sorting the entire image, so there
  /// must be a ordering defined for the pixel type in the image.
  /// That is, you must have defined a operator<, operator> and
  /// operator== for the pixel class you are using.  If you are using
  /// a built-in numerical type for you pixel type, you get this for
  /// free. Sorting the image is time consuming, so this operation is
  /// not recommended if performance is important.
  template <class ViewT>
  typename ViewT::pixel_type median_pixel_value(const ImageViewBase<ViewT> &view_) {
    const ViewT& view = copy(view_).impl();

    typedef typename ViewT::iterator iterator;
    typedef typename ViewT::pixel_type pixel_type;
    typedef typename PixelChannelType<typename ViewT::pixel_type>::type channel_type;

    sort(view.begin(), view.end());
    
    iterator median_pixel = view.begin() + (view.rows() * view.cols() / 2);
    return *median_pixel;
  }

}  // namespace vw


#endif // __VW_IMAGE_STATISTICS_H__
