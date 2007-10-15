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
#include <vw/Image/ImageView.h>
#include <vw/Image/Manipulation.h>

namespace vw {

  template <class ViewT>
  class MinMaxPixelAccumulator {
    typedef typename ViewT::pixel_type pixel_type;
    typedef typename PixelWithoutAlpha<pixel_type>::type value_type;
    value_type minval, maxval;
    bool valid;
  public:
    MinMaxPixelAccumulator() : valid(false) {}

    void operator()( pixel_type const& pix ) {
      if (is_transparent(pix)) return;
      value_type value(pix);
      if (!valid) {
        minval = maxval = value;
        valid = true;
      }
      else {
        if( value < minval ) minval = value;
        if( value > maxval ) maxval = value;
      }
    }

    value_type minimum() const { return minval; }
    value_type maximum() const { return maxval; }
  };

  /// Compute the minimum value of all valid pixels in the image.
  template <class ViewT>
  typename PixelWithoutAlpha<typename ViewT::pixel_type>::type
  min_pixel_value( const ImageViewBase<ViewT>& view ) {
    MinMaxPixelAccumulator<ViewT> accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.minimum();
  }

  /// Compute the maximum value of all valid pixels in the image.
  template <class ViewT>
  typename PixelWithoutAlpha<typename ViewT::pixel_type>::type
  max_pixel_value( const ImageViewBase<ViewT>& view ) {
    MinMaxPixelAccumulator<ViewT> accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.maximum();
  }

  /// Simultaneously compute the minimum and maximum values of all
  /// valid pixels in the image.
  template <class ViewT>
  void min_max_pixel_values( const ImageViewBase<ViewT> &view, 
                             typename PixelWithoutAlpha<typename ViewT::pixel_type>::type &min, 
                             typename PixelWithoutAlpha<typename ViewT::pixel_type>::type &max )
  {
    MinMaxPixelAccumulator<ViewT> accumulator;
    for_each_pixel( view, accumulator );
    min = accumulator.minimum();
    max = accumulator.maximum();
  }


  template <class ViewT>
  class MinMaxChannelAccumulator {
    typedef typename ViewT::pixel_type pixel_type;
    typedef typename CompoundChannelType<typename ViewT::pixel_type>::type channel_type;
    int num_channels;
    channel_type minval, maxval;
    bool valid;
  public:
    MinMaxChannelAccumulator()
      : num_channels( PixelNumChannels<pixel_type>::value - (PixelHasAlpha<pixel_type>::value ? 1 : 0) ),
        valid(false) {}

    void operator()( pixel_type const& pix ) {
      if (is_transparent(pix)) return;
      if (!valid) {
        minval = maxval = compound_select_channel<channel_type>(pix,0);
        valid = true;
      }
      for (int channel = 0; channel < num_channels; channel++) {
        channel_type channel_value = compound_select_channel<channel_type>(pix,channel);
        if( channel_value < minval ) minval = channel_value;
        if( channel_value > maxval ) maxval = channel_value;
      }
    }

    channel_type minimum() const { 
      VW_ASSERT(valid, ArgumentErr() << "MinMaxChannelAccumulator: no valid pixels");
      return minval; 
    }

    channel_type maximum() const {
      VW_ASSERT(valid, ArgumentErr() << "MinMaxChannelAccumulator: no valid pixels");
      return maxval;
    }
  };

  /// Compute the minimum value stored in all of the channels of all of the planes of the images.
  template <class ViewT>
  typename CompoundChannelType<typename ViewT::pixel_type>::type
  min_channel_value( const ImageViewBase<ViewT>& view ) {
    MinMaxChannelAccumulator<ViewT> accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.minimum();
  }

  /// Compute the maximum value stored in all of the channels of all of the planes of the images.
  template <class ViewT>
  typename CompoundChannelType<typename ViewT::pixel_type>::type
  max_channel_value( const ImageViewBase<ViewT>& view ) {
    MinMaxChannelAccumulator<ViewT> accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.maximum();
  }

  /// Simultaneously compute the min and max value in all of the
  /// channels of all of the planes of the image.
  template <class ViewT>
  void min_max_channel_values( const ImageViewBase<ViewT> &view, 
                               typename CompoundChannelType<typename ViewT::pixel_type>::type &min, 
                               typename CompoundChannelType<typename ViewT::pixel_type>::type &max )
  {
    MinMaxChannelAccumulator<ViewT> accumulator;
    for_each_pixel( view, accumulator );
    min = accumulator.minimum();
    max = accumulator.maximum();
  }


  template <class ViewT>
  class PixelSumAccumulator {
    typedef typename ViewT::pixel_type pixel_type;
    typedef typename PixelChannelCast<pixel_type,double>::type accum_type;
    accum_type accum;
  public:
    PixelSumAccumulator() : accum() {}

    void operator()( pixel_type const& pix ) {
      if (!is_transparent(pix)) accum += pix;
    }

    accum_type sum() const {
      return accum;
    }
  };

  /// Compute the sum of all valid pixels in the image.
  template <class ViewT>
  typename PixelChannelCast<typename ViewT::pixel_type,double>::type
  sum_of_pixel_values( const ImageViewBase<ViewT>& view ) {
    PixelSumAccumulator<ViewT> accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.sum();
  }


  template <class ViewT>
  class ChannelSumAccumulator {
    typedef typename ViewT::pixel_type pixel_type;
    typedef typename PixelChannelType<pixel_type>::type channel_type;
    typedef typename AccumulatorType<channel_type>::type accum_type;
    accum_type accum;
  public:
    ChannelSumAccumulator() : accum() {}

    void operator()( pixel_type const& pix ) {
      if ( ! is_transparent(pix) ) {
        compound_apply_in_place( *this, non_alpha_channels( pix ) );
      }
    }

    void operator()( channel_type const& val ) {
      accum += val;
    }

    accum_type sum() const { return accum; }
  };

  /// Compute the sum of all channels (other than any alpha channel)
  /// of all the valid (nonzero-alpha) pixels of the image.
  template <class ViewT>
  typename AccumulatorType<typename CompoundChannelType<typename ViewT::pixel_type>::type>::type
  sum_of_channel_values( const ImageViewBase<ViewT>& view ) {
    ChannelSumAccumulator<ViewT> accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.sum();
  }


  template <class ViewT>
  class PixelMeanAccumulator {
    typedef typename ViewT::pixel_type pixel_type;
    typedef typename PixelChannelType<pixel_type>::type channel_type;
    typedef typename PixelChannelCast<pixel_type,double>::type accum_type;
    accum_type accum;
    double num_samples;
  public:
    PixelMeanAccumulator( ViewT const& view )
      : accum(), num_samples( (double) view.planes() * view.rows() * view.cols() ) {}

    void operator()( pixel_type const& pix ) {
      if (is_transparent(pix)) return;
      accum += pix;
    }

    accum_type mean() const {
      VW_ASSERT(num_samples!=0, ArgumentErr() << "mean_pixel_value(): the image contained no pixels.");
      return accum / num_samples;
    }

    accum_type weighted_mean() const {
      double denominator = num_samples;
      if (PixelHasAlpha<pixel_type>::value) {
        denominator = alpha_channel(accum) / ChannelRange<channel_type>::max();
      }
      VW_ASSERT(denominator!=0, ArgumentErr() << "mean_pixel_value(): the image contained zero valid pixels.");
      return accum / denominator;
    }
  };

  /// Computes the mean of the values of all the pixels of all of
  /// the planes of an image.  For images that have an alpha channel,
  /// pixels with zero alpha are treated as zero.
  ///
  template <class ViewT>
  typename PixelChannelCast<typename ViewT::pixel_type,double>::type
  mean_pixel_value( const ImageViewBase<ViewT> &view_ ) {
    const ViewT& view = view_.impl();
    PixelMeanAccumulator<ViewT> accumulator( view );
    for_each_pixel( view, accumulator );
    return accumulator.mean();
  }

  /// Computes the weighted mean of the values of all the pixels of
  /// all of the planes of an image, using the alpha channel as a
  /// weight and assuming pre-multiplied pixel values.  This function
  /// throws an ArgumentErr() exception if the image has zero size or
  /// is completely transparent.  For images with no alpha, this 
  /// function is identical to mean_pixel_value().
  template <class ViewT>
  typename PixelChannelCast<typename ViewT::pixel_type,double>::type
  weighted_mean_pixel_value( const ImageViewBase<ViewT> &view_ ) {
    const ViewT& view = view_.impl();
    PixelMeanAccumulator<ViewT> accumulator( view );
    for_each_pixel( view, accumulator );
    return accumulator.weighted_mean();
  }

  /// Computes the mean of the values of all the channels of all of
  /// the pixels of an image, other than any alpha channel.  For 
  /// images that do have an alpha channel, pixels with zero alpha 
  /// are treated as zero.
  ///
  template <class ViewT>
  double mean_channel_value( const ImageViewBase<ViewT> &view ) {
    return mean_channel_value( non_alpha_channels( mean_pixel_value( view ) ) );
  }

  /// Computes the weighted mean of the values of all the channels of
  /// all of the pixels of an image, other than any alpha channel,
  /// using the alpha channel as a weight and assuming pre-multiplied
  /// pixel values.  This function throws an ArgumentErr() exception
  /// if the image has zero size or is completely transparent.  For
  /// images with no alpha, this function is identical to
  /// mean_channel_value().  For images that do have an alpha channel,
  /// pixels with zero alpha are treated as zero.
  ///
  template <class ViewT>
  double weighted_mean_channel_value( const ImageViewBase<ViewT> &view ) {
    return mean_channel_value( non_alpha_channels( weighted_mean_pixel_value( view ) ) );
  }

  template <class ViewT>
  class PixelStdDevAccumulator {
    typedef typename ViewT::pixel_type pixel_type;
    typedef typename PixelChannelType<pixel_type>::type channel_type;
    typedef typename PixelChannelCast<pixel_type,double>::type value_type;
    typedef typename PixelWithoutAlpha<pixel_type>::type accum_type;
    accum_type mom1_accum, mom2_accum;
    double num_samples;
  public:
    PixelStdDevAccumulator( ViewT const& view )
      : mom1_accum(), mom2_accum(), num_samples( 0 )
    {
      if ( ! PixelHasAlpha<pixel_type>::value ) {
        num_samples = (double) view.planes() * view.rows() * view.cols();
      }
    }

    template <class PixelT>
    typename boost::disable_if<PixelHasAlpha<PixelT> >::type operator()( PixelT const& pix ) {
      value_type value(pix);
      mom1_accum += value;
      mom2_accum += value * value;
    }

    template <class PixelT>
    typename boost::enable_if<PixelHasAlpha<PixelT> >::type operator()( PixelT const& pix ) {
      double weight = (double) select_alpha(pix) / ChannelRange<channel_type>::max();
      if( weight == 0 ) return;
      value_type value(pix);
      mom1_accum += value;
      mom2_accum += value * value / weight;
      num_samples += weight;
    }

    value_type stddev() const {
      VW_ASSERT(num_samples!=0, ArgumentErr() << "stddev_pixel_value(): the image contained zero valid pixels.");
      // FIXME: I don't think plain sqrt() works for compound pixel types....
      return sqrt(mom2_accum/num_samples - (mom1_accum/num_samples)*(mom1_accum/num_samples));
    }
  };

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
  double stddev_pixel_value( const ImageViewBase<ViewT> &view ) {
    PixelStdDevAccumulator<ViewT> accumulator( view );
    for_each_pixel( view, accumulator );
    return accumulator.stddev();
  }


  template <class ViewT>
  class ChannelStdDevAccumulator {
    double mom1_accum, mom2_accum, num_samples, num_channels;
    typedef typename ViewT::pixel_type pixel_type;
    typedef typename CompoundChannelType<pixel_type>::type channel_type;
  public:
    ChannelStdDevAccumulator( ViewT const& view )
      : mom1_accum(0), mom2_accum(0), num_samples(0), num_channels(view.channels())
    {
      if( ! PixelHasAlpha<pixel_type>::value ) {
        num_samples = (double) view.planes() * view.rows() * view.cols() * num_channels;
      }
      else {
        num_channels -= 1;
      }
    }

    template <class PixelT>
    typename boost::disable_if<PixelHasAlpha<PixelT> >::type operator()( PixelT const& pix ) {
      for (int32 channel = 0; channel < num_channels; channel++) {
        channel_type channel_value = compound_select_channel<channel_type>(pix,channel);
        mom1_accum += channel_value;
        mom2_accum += (double) channel_value * channel_value;
      }
    }

    template <class PixelT>
    typename boost::enable_if<PixelHasAlpha<PixelT> >::type operator()( PixelT const& pix ) {
      double weight = (double) compound_select_channel<channel_type>(pix,num_channels)
        / ChannelRange<channel_type>::max();
      if( weight == 0 ) return;
      for (int32 channel = 0; channel < num_channels; channel++) {
        channel_type channel_value = compound_select_channel<channel_type>(pix,channel);
        mom1_accum += channel_value;
        mom2_accum += (double) channel_value * channel_value / weight;
      }
      num_samples += weight * num_channels;
    }

    double stddev() const {
      if (num_samples == 0) 
        vw_throw(ArgumentErr() << "stddev_channel_value(): the image contained zero valid pixels.");
      return sqrt(mom2_accum/num_samples - (mom1_accum/num_samples)*(mom1_accum/num_samples));
    }
  };

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
  double stddev_channel_value( const ImageViewBase<ViewT> &view_ ) {
    const ViewT& view = view_.impl();
    ChannelStdDevAccumulator<ViewT> accumulator( view );
    for_each_pixel( view, accumulator );
    return accumulator.stddev();
  }


  template <class ViewT>
  class MedianPixelAccumulator {
    typedef typename ViewT::pixel_type pixel_type;
    typedef typename PixelWithoutAlpha<pixel_type>::type value_type;
    std::vector<value_type> values;
  public:
    void operator()( pixel_type const& pix ) {
      values.push_back( value_type( pix ) );
    }

    value_type median() {
      VW_ASSERT(values.size()!=0, ArgumentErr() << "pixel_median_value(): the image contained zero valid pixels.");
      sort(values.begin(),values.end());
      return values[values.size()/2];
    }
  };
  
  /// Compute the median pixel value of an image.  If the image has an
  /// alpha channel it is used as a valid data mask and the median
  /// value is returned with no alpha channel.  That is, if the image
  /// has pixel type PixelRGBA, this function will return PixelRGB.
  /// It computes the median by sorting all the valid pixels in the
  /// image, so there must be a ordering defined for the pixel type.
  /// That is, you must have defined a operator<, operator> and
  /// operator== for the pixel type you are using.  If you are using a
  /// built-in numerical type for you pixel type, you get this for
  /// free. Sorting the image is time consuming, so this operation is
  /// not recommended if performance is important.
  template <class ViewT>
  typename PixelWithoutAlpha<typename ViewT::pixel_type>::type
  median_pixel_value( const ImageViewBase<ViewT> &view ) {
    MedianPixelAccumulator<ViewT> accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.median();
  }


  template <class ViewT>
  class MedianChannelAccumulator {
    typedef typename ViewT::pixel_type pixel_type;
    typedef typename PixelChannelType<pixel_type>::type channel_type;
    mutable std::vector<channel_type> values;
    int num_channels;
  public:
    MedianChannelAccumulator() : num_channels( PixelNumChannels<pixel_type>::value ) {
      if( PixelHasAlpha<pixel_type>::value ) num_channels -= 1;
    }

    void operator()( pixel_type const& pix ) {
      if (is_transparent(pix)) return;
      for (int i=0; i<num_channels; ++i) {
        values.push_back( compound_select_channel<channel_type>(pix,i) );
      }
    }

    channel_type median() const {
      VW_ASSERT(values.size()!=0, ArgumentErr() << "median_channel_value(): the image contained zero valid pixels.");
      sort(values.begin(),values.end());
      return values[values.size()/2];
    }
  };
  
  /// Computes the median channel value of an image.  If the image has
  /// an alpha channel it is used as a valid data mask.  This function
  /// computes the median by sorting all the channel values of all the
  /// valid pixels in the image, which is time consuming, so this
  /// operation is not recommended if performance is important.
  template <class ViewT>
  typename PixelChannelType<typename ViewT::pixel_type>::type
  median_channel_value( const ImageViewBase<ViewT> &view ) {
    MedianChannelAccumulator<ViewT> accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.median();
  }

}  // namespace vw

#endif // __VW_IMAGE_STATISTICS_H__
