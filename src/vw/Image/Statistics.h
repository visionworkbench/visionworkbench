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


/// \file Statistics.h
///
/// These functions compute a variety of statistics on a per-pixel
/// or per-channel basis.  At the moment we provide the following
/// functions:
///
/// - min_channel_value
/// - max_channel_value
/// - min_max_channel_values
/// - sum_of_channel_values
/// - mean_channel_value
/// - stddev_channel_value
/// - median_channel_value
/// - weighted_mean_channel_value
///
/// - min_pixel_value
/// - max_pixel_value
/// - min_max_pixel_values
/// - sum_of_pixel_values
/// - mean_pixel_value
/// - stddev_pixel_value
/// - median_pixel_value
/// - weighted_mean_pixel_value

#ifndef __VW_IMAGE_STATISTICS_H__
#define __VW_IMAGE_STATISTICS_H__

#include <boost/type_traits.hpp>

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/PixelMask.h>

namespace vw {

  // CHANNEL operations
  //////////////////////////////////////////////////

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

  /// Compute the sum of all the channels of all the valid pixels of the image.
  template <class ViewT>
  typename AccumulatorType<typename PixelChannelType<typename ViewT::pixel_type>::type>::type
  sum_of_channel_values( const ImageViewBase<ViewT>& view ) {
    typedef typename AccumulatorType<typename PixelChannelType<typename ViewT::pixel_type>::type>::type accum_type;
    ChannelAccumulator<Accumulator<accum_type> > accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.value();
  }

  /// Computes the mean of the values of the channels of all of the
  /// valid (non-masked) pixels of an image (including alpha but
  /// excluding mask channels).
  template <class ViewT>
  double mean_channel_value( const ImageViewBase<ViewT> &view ) {
    typedef typename PixelChannelType<typename ViewT::pixel_type>::type accum_type;
    ChannelAccumulator<MeanAccumulator<accum_type> > accumulator;
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

  // PIXEL operations
  //////////////////////////////////

  // Special Element-wise accumulators
  template <class ValT>
  class EWMinMaxAccumulator : public ReturnFixedType<void> {
    ValT m_min, m_max;
    bool m_valid;
  public:
    EWMinMaxAccumulator() : m_valid(false) {}

    void operator()( ValT const& arg ) {
      if ( !m_valid ) {
        m_min = m_max = arg;
        m_valid = true;
      } else
        for ( size_t i = 0; i < CompoundNumChannels<ValT>::value; i++ )
          if ( arg[i] < m_min[i] ) m_min[i] = arg[i];
          else if ( arg[i] > m_max[i] ) m_max[i] = arg[i];
    }

    bool is_valid() const { return m_valid; }

    ValT minimum() const {
      VW_ASSERT(m_valid, ArgumentErr() << "EWMinMaxAccumulator: no valid samples" );
      return m_min;
    }

    ValT maximum() const {
      VW_ASSERT(m_valid, ArgumentErr() << "EWMinMaxAccumulator: no valid samples" );
      return m_max;
    }
  };

  template <class ValT>
  class EWStdDevAccumulator : public ReturnFixedType<void> {
    typedef typename PixelChannelType<ValT>::type channel_type;
    std::vector<channel_type> m_sum;
    std::vector<channel_type> m_sum_2;
    double num_samples;
  public:
    EWStdDevAccumulator () : num_samples(0) {
      m_sum.resize( CompoundNumChannels<ValT>::value );
      m_sum_2.resize( CompoundNumChannels<ValT>::value );
      for ( vw::int32 i = 0; i < CompoundNumChannels<ValT>::value; i++ )
        m_sum[i] = m_sum_2[i] = 0;
    }

    void operator()( ValT const& value ) {
      num_samples++;
      for ( vw::int32 i = 0; i < CompoundNumChannels<ValT>::value; i++ ) {
        m_sum[i] += value[i];
        m_sum_2[i] += value[i]*value[i];
      }
    }

    ValT value() const {
      VW_ASSERT(num_samples, ArgumentErr() << "EWStdDevAccumulator(): no valid samples.");
      ValT result;
      for ( vw::int32 i = 0; i < CompoundNumChannels<ValT>::value; i++ )
        result[i] = sqrt(m_sum_2[i]/num_samples - (m_sum[i]/num_samples)*(m_sum[i]/num_samples));
      return result;
    }
  };

  template <class ValT>
  class EWMedianAccumulator : public ReturnFixedType<void> {
    typedef std::vector<std::vector<typename PixelChannelType<ValT>::type> > storage_type;
    storage_type m_values;
  public:
    EWMedianAccumulator() {
      m_values.resize( CompoundNumChannels<ValT>::value );
    }

    void operator()( ValT const& value ) {
      for ( vw::int32 i = 0; i < CompoundNumChannels<ValT>::value; i++ )
        m_values[i].push_back( value[i] );
    }

    ValT value() {
      VW_ASSERT(m_values[0].size(), ArgumentErr() << "MedianAccumulator: no valid samples");
      ValT result;
      for ( vw::int32 i = 0; i < CompoundNumChannels<ValT>::value; i++ ) {
        result[i] = destructive_median(m_values[i]);
      }
      return result;
    }
  };

  template <class AccumT>
  class PixelAccumulator : public AccumT {
  public:
    template <class ArgT>
    void operator()( ArgT const& pix ) {
      if ( ::vw::is_valid(pix) )
        AccumT::operator()( remove_mask( pix ) );
    }
  };

  // Functions

  template <class ViewT>
  typename boost::enable_if< IsCompound<typename UnmaskedPixelType<typename ViewT::pixel_type>::type>,typename UnmaskedPixelType<typename ViewT::pixel_type>::type>::type
  min_pixel_value( ImageViewBase<ViewT> const& view ) {
    typedef typename UnmaskedPixelType<typename ViewT::pixel_type>::type accum_type;
    PixelAccumulator<EWMinMaxAccumulator<accum_type> > accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.minimum();
  }

  template <class ViewT>
    typename boost::disable_if< IsCompound<typename UnmaskedPixelType<typename ViewT::pixel_type>::type>,typename PixelChannelType<typename ViewT::pixel_type>::type>::type
  min_pixel_value( ImageViewBase<ViewT> const& view ) {
    return min_channel_value( view );
  }

  template <class ViewT>
    typename boost::enable_if< IsCompound<typename UnmaskedPixelType<typename ViewT::pixel_type>::type>,typename UnmaskedPixelType<typename ViewT::pixel_type>::type>::type
  max_pixel_value( ImageViewBase<ViewT> const& view ) {
    typedef typename UnmaskedPixelType<typename ViewT::pixel_type>::type accum_type;
    PixelAccumulator<EWMinMaxAccumulator<accum_type> > accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.maximum();
    return max_channel_value( view );
  }

  template <class ViewT>
    typename boost::disable_if< IsCompound<typename UnmaskedPixelType<typename ViewT::pixel_type>::type>,typename PixelChannelType<typename ViewT::pixel_type>::type>::type
  max_pixel_value( ImageViewBase<ViewT> const& view ) {
    return max_channel_value( view );
  }

  template <class ViewT>
  void min_max_pixel_values( ImageViewBase<ViewT> const& view,
                             typename boost::enable_if< IsCompound<typename UnmaskedPixelType<typename ViewT::pixel_type>::type>, typename UnmaskedPixelType<typename ViewT::pixel_type>::type>::type &min,
                             typename UnmaskedPixelType<typename ViewT::pixel_type>::type &max ) {
    typedef typename UnmaskedPixelType<typename ViewT::pixel_type>::type accum_type;
    PixelAccumulator<EWMinMaxAccumulator<accum_type> > accumulator;
    for_each_pixel( view, accumulator );
    min = accumulator.minimum();
    max = accumulator.maximum();
  }

  template <class ViewT>
  void min_max_pixel_values( ImageViewBase<ViewT> const& view,
                             typename boost::disable_if< IsCompound<typename UnmaskedPixelType<typename ViewT::pixel_type>::type>, typename PixelChannelType<typename ViewT::pixel_type>::type>::type &min,
                             typename PixelChannelType<typename ViewT::pixel_type>::type &max ) {
    min_max_channel_values( view, min, max );
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

  template <class ViewT>
  typename PixelChannelCast<typename UnmaskedPixelType<typename ViewT::pixel_type>::type,double>::type
  mean_pixel_value( ImageViewBase<ViewT> const& view ) {
    typedef typename PixelChannelCast<typename UnmaskedPixelType<typename ViewT::pixel_type>::type,double>::type accum_type;
    PixelAccumulator<MeanAccumulator<accum_type> > accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.value();
  }

  template <class ViewT>
  typename boost::enable_if< IsCompound<typename UnmaskedPixelType<typename ViewT::pixel_type>::type>,typename UnmaskedPixelType<typename ViewT::pixel_type>::type>::type
  stddev_pixel_value( ImageViewBase<ViewT> const& view ) {
    typedef typename UnmaskedPixelType<typename ViewT::pixel_type>::type accum_type;
    PixelAccumulator<EWStdDevAccumulator<accum_type> > accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.value();
  }

  template <class ViewT>
  typename boost::disable_if< IsCompound<typename UnmaskedPixelType<typename ViewT::pixel_type>::type>,double>::type
  stddev_pixel_value( ImageViewBase<ViewT> const& view ) {
    return stddev_channel_value( view );
  }

  template <class ViewT>
    typename boost::enable_if< IsCompound<typename UnmaskedPixelType<typename ViewT::pixel_type>::type>,typename UnmaskedPixelType<typename ViewT::pixel_type>::type>::type
  median_pixel_value( ImageViewBase<ViewT> const& view ) {
    typedef typename UnmaskedPixelType<typename ViewT::pixel_type>::type accum_type;
    PixelAccumulator<EWMedianAccumulator<accum_type> > accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.value();
  }

  template <class ViewT>
    typename boost::disable_if< IsCompound<typename UnmaskedPixelType<typename ViewT::pixel_type>::type>,typename PixelChannelType<typename ViewT::pixel_type>::type>::type
  median_pixel_value( ImageViewBase<ViewT> const& view ) {
    return median_channel_value( view );
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

  // Find the histogram of an image. 
  template <class ViewT>
  void histogram( const ImageViewBase<ViewT> &view, int num_bins, std::vector<double> &hist){
    
    VW_ASSERT(num_bins > 0, ArgumentErr() << "histogram: number of input bins must be positive");
    
    // Find the maximum and minimum
    double max_val = -std::numeric_limits<double>::max(), min_val = -max_val;
    for (int row = 0; row < view.impl().rows(); row++){
      for (int col = 0; col < view.impl().cols(); col++){
        if ( !is_valid(view.impl()(col, row)) ) continue;
        double val = view.impl()(col, row);
        if (val < min_val) min_val = val;
        if (val > max_val) max_val = val;
      }
    }
    if (max_val == min_val) max_val = min_val + 1.0;
      
    hist.assign(num_bins, 0.0);
    for (int row = 0; row < view.impl().rows(); row++){
      for (int col = 0; col < view.impl().cols(); col++){
        if ( !is_valid(view.impl()(col, row)) ) continue;
        double val = view.impl()(col, row);
        int bin = (int)round( (num_bins - 1) * ( (val - min_val)/(max_val - min_val) ) );
        hist[bin]++;
      }
    }

    return;
  }


  /// Find the min and max values in an image
  /// - TODO: Why are there two methods for doing this?
  template <class ViewT>
  void find_image_min_max( const ImageViewBase<ViewT> &view, double &min_val, double &max_val){

    max_val = -std::numeric_limits<double>::max();
    min_val = -max_val;
    for (int row = 0; row < view.impl().rows(); row++){
      for (int col = 0; col < view.impl().cols(); col++){
        if ( !is_valid(view.impl()(col, row)) ) 
          continue;
        double val = view.impl()(col, row);
        if (val < min_val) min_val = val;
        if (val > max_val) max_val = val;
      }
    }
  }

  /// Overload that takes precomputed min and max values
  template <class ViewT>
  void histogram( const ImageViewBase<ViewT> &view, int num_bins, double min_val, double max_val,
                  std::vector<double> &hist){
    
    VW_ASSERT(num_bins > 0, ArgumentErr() << "histogram: number of input bins must be positive");
    
    if (max_val == min_val) max_val = min_val + 1.0;
      
    hist.assign(num_bins, 0.0);
    for (int row = 0; row < view.impl().rows(); row++){
      for (int col = 0; col < view.impl().cols(); col++){
        if ( !is_valid(view.impl()(col, row)) ) continue;
        double val = view.impl()(col, row);
        int bin = (int)round( (num_bins - 1) * ( (val - min_val)/(max_val - min_val) ) );
        hist[bin]++;
      }
    }

    return;
  }

  /// Returns the histogram bin index corresponding to the specified percentile
  inline size_t get_histogram_percentile(std::vector<double> const& hist, double percentile) {
  
    // Verify the input percentile is in the legal range
    if ((percentile < 0) || (percentile > 1.0)) {
      vw_throw(ArgumentErr() << "get_histogram_percentile: illegal percentile request: " << percentile << "\n");
    }

    // Get the total pixel count
    const size_t num_bins = hist.size();
    double num_pixels = 0;
    for (size_t i=0; i<num_bins; ++i)
      num_pixels += hist[i];

    // Now go through and find the requested percentile
    double running_percentile = 0;
    for (size_t i=0; i<num_bins; ++i) {
      double this_percent = hist[i] / num_pixels;
      running_percentile += this_percent;
      //std::cout << "p - " << running_percentile << std::endl;
      if (running_percentile >= percentile)
        return i;
    }
    vw_throw(LogicErr() << "get_histogram_percentile: Illegal histogram encountered!");
  }

  // Find the optimal Otsu threshold for splitting a gray scale image
  // into black and white pixels.
  // Reference: http://www.labbookpages.co.uk/software/imgProc/otsuThreshold.html
  // It agrees with Matlab for uint8 images.
  // Note: The optimal threshold is normalized, so it is between 0 and 1.
  // To find its value in pixel space, you may want to compute:
  // minImageVal + threshold*(maxImageVal - minImageVal).
  template <class ViewT>
  double optimal_threshold( const ImageViewBase<ViewT> &view){

    std::vector<double> hist;
    int num_bins = 256;
    histogram(view, num_bins, hist);

    double sum = 0.0;
    for (size_t i = 0; i < hist.size(); i++) sum += hist[i];
    if (sum == 0.0) return 0.0;

    // Normalize the histogram
    for (size_t i = 0; i < hist.size(); i++) hist[i] /= sum;
    sum = 1.0;
    
    double totalAccum = 0.0;
    for (size_t i = 0; i < hist.size(); i++) totalAccum += i*hist[i];

    // Find the variance between classes
    std::vector<double> V;
    V.assign(num_bins, 0.0);

    double leftProb = 0.0, leftAccum = 0.0, rightProb = 0.0, rightAccum = 0.0;
    for (size_t i = 0; i < hist.size(); i++){

      leftProb  += hist[i];
      leftAccum += i*hist[i];

      rightProb  = sum - leftProb;
      rightAccum = totalAccum - leftAccum;

      if (leftProb == 0 || rightProb == 0.0) continue;
      
      double leftMean = leftAccum/leftProb;
      double rightMean = rightAccum/rightProb;
      V[i] = leftProb*rightProb*(leftMean-rightMean)*(leftMean-rightMean);
      
    }

    double maxV = *std::max_element(V.begin(), V.end());

    // If the maximum is reached at more than one index, find the average index
    double indexSum = 0, numIndices = 0;
    for (size_t i = 0; i < V.size(); i++){
      if (V[i] == maxV){
        indexSum += i;
        numIndices++;
      }
    }
    double meanIndex = indexSum/numIndices;

    // Normalize the value
    return meanIndex/(num_bins - 1.0);
  }
  
}  // namespace vw

#endif // __VW_IMAGE_STATISTICS_H__
