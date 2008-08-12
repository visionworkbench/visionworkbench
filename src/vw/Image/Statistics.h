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
/// These functions compute a variety of statistics on a per-pixel 
/// or per-channel basis.  At the moment we provide the following 
/// functions:
///
/// - min_pixel_value
/// - max_pixel_value
/// - min_max_pixel_values
/// - min_channel_value
/// - max_channel_value
/// - min_max_channel_values
/// - sum_of_pixel_values
/// - sum_of_channel_values
/// - mean_pixel_value
/// - mean_channel_value
/// - stddev_pixel_value
/// - stddev_channel_value
/// - median_pixel_value
/// - median_channel_value
///

#ifndef __VW_IMAGE_STATISTICS_H__
#define __VW_IMAGE_STATISTICS_H__

#include <boost/type_traits.hpp>

#include <vw/Core/FundamentalTypes.h>
#include <vw/Math/Functors.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/PixelMask.h>

namespace vw {

  /// An adapter to help applying an accumulator to all valid pixels
  /// in an image.
  template <class AccumT>
  class PixelAccumulator : public AccumT {
  public:
    template <class ArgT>
    void operator()( ArgT const& pix ) {
      if ( is_valid(pix) ) {
	AccumT::operator()( remove_mask( pix ) );
      }
    }
  };

  /// An adapter to help applying an accumulator to all channels of
  /// all valid pixels in an image.
  template <class AccumT>
  class ChannelAccumulator : public AccumT {
  public:
    template <class ArgT>
    void operator()( ArgT const& pix ) {
      if ( is_valid(pix) ) {
        compound_apply_in_place( (AccumT&)*this, remove_mask( pix ) );
      }
    }
  };


  /// Compute the minimum value of all valid pixels in the image.
  template <class ViewT>
  typename UnmaskedPixelType<typename ViewT::pixel_type>::type
  min_pixel_value( const ImageViewBase<ViewT>& view ) {
    typedef typename UnmaskedPixelType<typename ViewT::pixel_type>::type accum_type;
    PixelAccumulator<MinMaxAccumulator<accum_type> > accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.minimum();
  }

  /// Compute the maximum value of all valid pixels in the image.
  template <class ViewT>
  typename UnmaskedPixelType<typename ViewT::pixel_type>::type
  max_pixel_value( const ImageViewBase<ViewT>& view ) {
    typedef typename UnmaskedPixelType<typename ViewT::pixel_type>::type accum_type;
    PixelAccumulator<MinMaxAccumulator<accum_type> > accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.maximum();
  }

  /// Simultaneously compute the minimum and maximum values of all
  /// valid pixels in the image.
  template <class ViewT>
  void min_max_pixel_values( const ImageViewBase<ViewT> &view, 
                             typename UnmaskedPixelType<typename ViewT::pixel_type>::type &min, 
                             typename UnmaskedPixelType<typename ViewT::pixel_type>::type &max )
  {
    typedef typename UnmaskedPixelType<typename ViewT::pixel_type>::type accum_type;
    PixelAccumulator<MinMaxAccumulator<accum_type> > accumulator;
    for_each_pixel( view, accumulator );
    min = accumulator.minimum();
    max = accumulator.maximum();
  }

  /// Compute the minimum value stored in all of the channels of all
  /// of the planes of the images.
  template <class ViewT>
  typename PixelChannelType<typename ViewT::pixel_type>::type
  min_channel_value( const ImageViewBase<ViewT>& view ) {
    typedef typename PixelChannelType<typename ViewT::pixel_type>::type accum_type;
    ChannelAccumulator<MinMaxAccumulator<accum_type> > accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.minimum();
  }

  /// Compute the maximum value stored in all of the channels of all
  /// of the planes of the images.
  template <class ViewT>
  typename PixelChannelType<typename ViewT::pixel_type>::type
  max_channel_value( const ImageViewBase<ViewT>& view ) {
    typedef typename PixelChannelType<typename ViewT::pixel_type>::type accum_type;
    ChannelAccumulator<MinMaxAccumulator<accum_type> > accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.maximum();
  }

  /// Simultaneously compute the min and max value in all of the
  /// channels of all of the planes of the image.
  template <class ViewT>
  void min_max_channel_values( const ImageViewBase<ViewT> &view, 
                               typename PixelChannelType<typename ViewT::pixel_type>::type &min, 
                               typename PixelChannelType<typename ViewT::pixel_type>::type &max )
  {
    typedef typename PixelChannelType<typename ViewT::pixel_type>::type accum_type;
    ChannelAccumulator<MinMaxAccumulator<accum_type> > accumulator;
    for_each_pixel( view, accumulator );
    min = accumulator.minimum();
    max = accumulator.maximum();
  }

  /// Compute the sum of all valid pixels in the image.
  template <class ViewT>
  typename PixelChannelCast<typename UnmaskedPixelType<typename ViewT::pixel_type>::type,typename AccumulatorType<typename PixelChannelType<typename ViewT::pixel_type>::type>::type>::type
  sum_of_pixel_values( const ImageViewBase<ViewT>& view ) {
    typedef typename PixelChannelCast<typename UnmaskedPixelType<typename ViewT::pixel_type>::type,typename AccumulatorType<typename PixelChannelType<typename ViewT::pixel_type>::type>::type>::type accum_type;
    PixelAccumulator<Accumulator<accum_type> > accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.value();
  }

  /// Compute the sum of all the channels of all the valid pixels of
  /// the image.
  template <class ViewT>
  typename AccumulatorType<typename PixelChannelType<typename ViewT::pixel_type>::type>::type
  sum_of_channel_values( const ImageViewBase<ViewT>& view ) {
    typedef typename AccumulatorType<typename PixelChannelType<typename ViewT::pixel_type>::type>::type accum_type;    
    ChannelAccumulator<Accumulator<accum_type> > accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.value();
  }


  /// Computes the mean of the values of all the valid pixels of an
  /// image.
  template <class ViewT>
  typename PixelChannelCast<typename UnmaskedPixelType<typename ViewT::pixel_type>::type,double>::type
  mean_pixel_value( const ImageViewBase<ViewT> &view ) {
    typedef typename UnmaskedPixelType<typename ViewT::pixel_type>::type accum_type;
    PixelAccumulator<MeanAccumulator<accum_type> > accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.value();
  }

  /// Computes the weighted mean of the values of all the pixels of an
  /// image, using the alpha channel as a weight and assuming
  /// pre-multiplied channel values.  This function throws an
  /// ArgumentErr() exception if the image has zero size or is
  /// completely transparent.  For images with no alpha or mask
  /// channel this function is identical to mean_pixel_value().
  template <class ViewT>
  typename PixelWithoutAlpha<typename PixelChannelCast<typename UnmaskedPixelType<typename ViewT::pixel_type>::type,double>::type>::type
  weighted_mean_pixel_value( const ImageViewBase<ViewT> &view ) {
    typedef typename PixelChannelCast<typename UnmaskedPixelType<typename ViewT::pixel_type>::type,double>::type accum_type;
    accum_type mean = mean_pixel_value( view );
    if ( PixelHasAlpha<typename ViewT::pixel_type>::value ) {
      double weight = alpha_channel( mean ) / ChannelRange<typename PixelChannelType<typename ViewT::pixel_type>::type>::max();
      VW_ASSERT(weight, ArgumentErr() << "weighted_mean_pixel_value(): no weighted samples");
      mean /= weight;
    }
    return non_alpha_channels(mean);
  }

  /// Computes the mean of the values of the channels of all of the
  /// valid (non-masked) pixels of an image (including alpha but
  /// excluding mask channels).
  template <class ViewT>
  double mean_channel_value( const ImageViewBase<ViewT> &view ) {
    return mean_channel_value( mean_pixel_value( view ) );
  }

  /// Computes the weighted mean of the values of the channels of all
  /// of the valid pixels of an image (excluding alpha or mask
  /// channels), using the alpha channel as a weight.  This function
  /// throws an ArgumentErr() exception if the image has zero size or
  /// is completely transparent.  For images with no alpha channel
  /// this function is identical to mean_channel_value().
  template <class ViewT>
  double weighted_mean_channel_value( const ImageViewBase<ViewT> &view ) {
    return mean_channel_value( weighted_mean_pixel_value( view ) );
  }


  /// Computes the standard deviation of the values of all the pixels
  /// of all of the planes of an image.  For images that have an alpha
  /// channel, this function computes the weighted standard deviation,
  /// using the alpha channel as a weight and assuming pre-multiplied
  /// pixel values.  It returns the result as a fully opaque pixel
  /// with a double channel type. This function throws an ArgumentErr() 
  /// exception if the image has zero size or is completely transparent.
  ///
  /// Note: This function computes the total stanadard deviation, not
  /// the sample standard deviation as was computed by previous
  /// versions.  If you need the sample standard deviation, just
  /// multiply the result by sqrt(num_samples/(num_samples-1)), where 
  /// num_samples=channels*cols*rows*planes.  Note that the concept 
  /// of sample standard deviation is not particularly meaningful for 
  /// images with alpha channels.
  ///
  template <class ViewT>
  typename PixelChannelCast<typename UnmaskedPixelType<typename ViewT::pixel_type>::type,double>::type 
  stddev_pixel_value( const ImageViewBase<ViewT> &view ) {
    PixelAccumulator<StdDevAccumulator<typename ViewT::pixel_type> > accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.value();
  }


  /// Computes the standard deviation of the values of all the
  /// channels of all of the planes of an image.  For images that have
  /// an alpha channel, this function computes the weighted standard
  /// deviation, using the alpha channel as a weight and assuming
  /// pre-multiplied pixel values.  This function throws an
  /// ArgumentErr() exception if the image has zero size or is
  /// completely transparent.
  ///
  /// Note: This function computes the total stanadard deviation, not
  /// the sample standard deviation as was computed by previous
  /// versions.  If you need the sample standard deviation, just
  /// multiply the result by sqrt(num_samples/(num_samples-1)), where 
  /// num_samples=channels*cols*rows*planes.  Note that the concept 
  /// of sample standard deviation is not particularly meaningful for 
  /// images with alpha channels.
  ///
  template <class ViewT>
  double stddev_channel_value( const ImageViewBase<ViewT> &view ) {
    typedef typename PixelChannelType<typename ViewT::pixel_type>::type channel_type;
    ChannelAccumulator<StdDevAccumulator<channel_type> > accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.value();
  }


  /// Compute the median pixel value of an image.  Only valid (e.g. 
  /// non-transparent) pixels are considered.  It computes the median 
  /// by sorting all the valid pixels in the image, so there must be 
  /// a ordering defined for the pixel type.  That is, you must have 
  /// defined a operator< (and operator== ???) for the pixel type you 
  /// are using.  If you are using a built-in numerical type for your 
  /// pixel type, you get this for free. Sorting the image is time-
  /// and memory-intensive, so this operation is not recommended for 
  /// large images.
  template <class ViewT>
  typename ViewT::pixel_type
  median_pixel_value( const ImageViewBase<ViewT> &view ) {
    typedef typename ViewT::pixel_type accum_type;
    PixelAccumulator<MedianAccumulator<accum_type> > accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.value();
  }

  /// Computes the median channel value of an image.  Only non-alpha
  /// channels of valid (e.g.  non-transparent) pixels are considered.
  /// This function computes the median by sorting all the channel
  /// values in the image, which is time- and memory-intensive, so
  /// this operation is not recommended for large images.
  template <class ViewT>
  typename PixelChannelType<typename ViewT::pixel_type>::type
  median_channel_value( const ImageViewBase<ViewT> &view ) {
    typedef typename PixelChannelType<typename ViewT::pixel_type>::type accum_type;
    ChannelAccumulator<MedianAccumulator<accum_type> > accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.value();
  }

}  // namespace vw

#endif // __VW_IMAGE_STATISTICS_H__
