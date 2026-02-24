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

#include <cfloat>
#include <boost/type_traits.hpp>
#include <vw/Math/Statistics.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/BlockImageOperator.h>
#include <vw/Image/EdgeExtension.h>

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
  double
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
      } else {
        for ( size_t i = 0; i < CompoundNumChannels<ValT>::value; i++ )
          if ( arg[i] < m_min[i] )
            m_min[i] = arg[i];
          else
            if ( arg[i] > m_max[i] )
              m_max[i] = arg[i];
      }
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

  /// This class wraps another accumulator functor to add pixel mask handling.
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
  typename boost::enable_if< IsCompound<typename UnmaskedPixelType<typename ViewT::pixel_type>::type>,
                                        typename UnmaskedPixelType<typename ViewT::pixel_type>::type>::type
  min_pixel_value( ImageViewBase<ViewT> const& view ) {
    typedef typename UnmaskedPixelType<typename ViewT::pixel_type>::type accum_type;
    PixelAccumulator<EWMinMaxAccumulator<accum_type> > accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.minimum();
  }

  template <class ViewT>
    typename boost::disable_if< IsCompound<typename UnmaskedPixelType<typename ViewT::pixel_type>::type>,
                                           typename PixelChannelType<typename ViewT::pixel_type>::type>::type
  min_pixel_value( ImageViewBase<ViewT> const& view ) {
    return min_channel_value( view );
  }

  template <class ViewT>
    typename boost::enable_if< IsCompound<typename UnmaskedPixelType<typename ViewT::pixel_type>::type>,
                                          typename UnmaskedPixelType<typename ViewT::pixel_type>::type>::type
  max_pixel_value( ImageViewBase<ViewT> const& view ) {
    typedef typename UnmaskedPixelType<typename ViewT::pixel_type>::type accum_type;
    PixelAccumulator<EWMinMaxAccumulator<accum_type> > accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.maximum();
    return max_channel_value( view );
  }

  template <class ViewT>
    typename boost::disable_if< IsCompound<typename UnmaskedPixelType<typename ViewT::pixel_type>::type>,
                                           typename PixelChannelType<typename ViewT::pixel_type>::type>::type
  max_pixel_value( ImageViewBase<ViewT> const& view ) {
    return max_channel_value( view );
  }

  template <class ViewT>
  void min_max_pixel_values( ImageViewBase<ViewT> const& view,
                             typename boost::enable_if< IsCompound<typename UnmaskedPixelType<typename ViewT::pixel_type>::type>, 
                             typename UnmaskedPixelType<typename ViewT::pixel_type>::type>::type &min,
                             typename UnmaskedPixelType<typename ViewT::pixel_type>::type        &max ) {
    typedef typename UnmaskedPixelType<typename ViewT::pixel_type>::type accum_type;
    PixelAccumulator<EWMinMaxAccumulator<accum_type> > accumulator;
    for_each_pixel( view, accumulator );
    min = accumulator.minimum();
    max = accumulator.maximum();
  }

  template <class ViewT>
  void min_max_pixel_values( ImageViewBase<ViewT> const& view,
                             typename boost::disable_if< IsCompound<typename UnmaskedPixelType<typename ViewT::pixel_type>::type>, 
                             typename PixelChannelType<typename ViewT::pixel_type>::type>::type &min,
                             typename PixelChannelType<typename ViewT::pixel_type>::type        &max ) {
    min_max_channel_values( view, min, max );
  }

  /// Compute the sum of all valid pixels in the image.
  template <class ViewT>
  typename PixelChannelCast<typename UnmaskedPixelType<typename ViewT::pixel_type>::type,
                            typename AccumulatorType<typename PixelChannelType<typename ViewT::pixel_type>::type>::type>::type
  sum_of_pixel_values( const ImageViewBase<ViewT>& view ) {
    typedef typename PixelChannelCast<typename UnmaskedPixelType<typename ViewT::pixel_type>::type,
                          typename AccumulatorType<typename PixelChannelType<typename ViewT::pixel_type>::type>::type>::type accum_type;
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
  typename boost::enable_if< IsCompound<typename UnmaskedPixelType<typename ViewT::pixel_type>::type>,
                                        typename UnmaskedPixelType<typename ViewT::pixel_type>::type>::type
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
    typename boost::enable_if< IsCompound<typename UnmaskedPixelType<typename ViewT::pixel_type>::type>,
                                          typename UnmaskedPixelType<typename ViewT::pixel_type>::type>::type
  median_pixel_value( ImageViewBase<ViewT> const& view ) {
    typedef typename UnmaskedPixelType<typename ViewT::pixel_type>::type accum_type;
    PixelAccumulator<EWMedianAccumulator<accum_type> > accumulator;
    for_each_pixel( view, accumulator );
    return accumulator.value();
  }

  template <class ViewT>
    typename boost::disable_if< IsCompound<typename UnmaskedPixelType<typename ViewT::pixel_type>::type>,
                                           typename PixelChannelType<typename ViewT::pixel_type>::type>::type
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

  /// Find the min and max values in an image
  template <class ViewT>
  void find_image_min_max( const ImageViewBase<ViewT> &view, double &min_val, double &max_val);

  /// Find histogram
  template <class ViewT>
  void histogram( const ImageViewBase<ViewT> &view, int num_bins, double min_val, double max_val,
                  math::Histogram &hist);

  // Find the optimal Otsu threshold for splitting a gray scale image
  // into black and white pixels.
  // Reference: http://www.labbookpages.co.uk/software/imgProc/otsuThreshold.html
  // This returns the scaled threshold, so 
  // minImageVal + normalized_threshold*(maxImageVal - minImageVal).
  template <class ViewT>
  double otsu_threshold(const ImageViewBase<ViewT> &view);

  // Another implementation of the Otsu threshold, from:
  // https://github.com/opencv/opencv/blob/master/modules/imgproc/src/thresh.cpp
  // This gives the same result as the one above if the number of bins
  // is 256, and num_sample_cols is same as number of cols, and
  // num_sample_rows is same as number of rows.
  template <class ViewT>
  double otsu_threshold(const ImageViewBase<ViewT> &view,
                        int num_sample_rows, int num_sample_cols,
                        int num_bins);
  
  /// Converts a single channel image into a uint8 image with percentile based intensity scaling.
  template <class ViewT>
  void percentile_scale_convert(ImageViewBase<ViewT> const& input_image,
                                ImageView<vw::uint8> &output_image,
                                double low_percentile=0.02, double high_percentile=0.98,
                                int num_bins=256);

  /// Converts a single channel image into a uint8 image using a
  /// simple stretch between the min and max values.
  template <class ViewT>
  void u8_convert(ImageViewBase<ViewT> const& input_image,
                  ImageView<PixelGray<vw::uint8> > &output_image) {
    // First get the min and max values
    double min_val, max_val;
    find_image_min_max(input_image, min_val, max_val);

    // Scale the image using the computed values and convert to uint8
    if (max_val == min_val)
      max_val = min_val + 1.0;
    output_image = pixel_cast<vw::uint8>(normalize( clamp(input_image, min_val, max_val),
                                                    min_val, max_val, 0.0, 255.0 ));
  }


  // TODO: Does this already exist somewhere in the code?
  /// An adapter to let a functor handle a single channel image mask.
  template <class AccumT>
  class SingleChannelAccumulator {
    AccumT * m_functor;
  public:
    
    SingleChannelAccumulator(AccumT* ptr) : m_functor(ptr) {}
    
    template <class ArgT>
    void operator()( ArgT const& pix ) {
      if ( is_valid(pix) )
        m_functor->operator()(remove_mask(pix));
    }
  };

  // TODO: Replace!
  template <typename T>
  struct PixelCollector {
    std::vector<T> m_vec;
    
    void operator()(T p) {
      m_vec.push_back(p);
    }
  };
  
  /// Thread safe functor to accumulate CDF results on multiple single channel images.
  /// - A CDF is computed for each image and they are then merged together.
  template<typename T>
  class ParallelCdfFunctor {

    typedef vw::math::CDFAccumulator<T> CdfType;

    CdfType * m_cdf_ptr;
    int m_subsample_amt;
    Mutex  m_mutex;

  public:

    /// Constructor takes a pointer to the CDF object that will be populated.
    ParallelCdfFunctor(CdfType* ptr, int subsample_amt=1)
      : m_cdf_ptr(ptr), m_subsample_amt(subsample_amt) {}

    /// Process an image and incorporate it into the input CDF object.
    template <class ImageT>
    void operator()(ImageView<ImageT> const& image, BBox2i const& bbox) {

      // Compute a CDF on just this input image.
      //CdfType local_cdf(1000, 20);
      //SingleChannelAccumulator<CdfType> accumulator(&local_cdf);

      // TODO: The CDF class cannot merge properly, so use this alternative
      //       method until it is fixed!
      typedef PixelCollector<float> PC; // TODO: Fix type!
      PC pixel_accum;
      float est_num_pixels = (image.rows()/m_subsample_amt)*(image.cols()/m_subsample_amt);
      pixel_accum.m_vec.reserve(int(est_num_pixels * 1.2));
      SingleChannelAccumulator<PC> accumulator(&pixel_accum);
      for_each_pixel( subsample( edge_extend(image, ConstantEdgeExtension()),
                                 m_subsample_amt ),
                      accumulator);

      // Merge the local CDF with the main CDF.
      m_mutex.lock();
      /*
      std::cout << "\nBLOCK val =";
      std::cout << "\nval = " << local_cdf.quantile(0);
      std::cout << "\nval = " << local_cdf.quantile(1); // Max
      std::cout << "\nval = " << local_cdf.approximate_mean();
      std::cout << "\nval = " << local_cdf.approximate_stddev();
      std::cout << "\nval = " << local_cdf.quantile(0.02); // Percentile values
      std::cout << "\nval = " << local_cdf.quantile(0.98) << std::endl;
      */
      //m_cdf_ptr->operator()(local_cdf);  // TODO: Fix this!
      for (size_t i=0; i<pixel_accum.m_vec.size(); ++i)
        m_cdf_ptr->operator()(pixel_accum.m_vec[i]);
      /*
      std::cout << "\nNEW val =";
      std::cout << "\nval = " << m_cdf_ptr->quantile(0);
      std::cout << "\nval = " << m_cdf_ptr->quantile(1); // Max
      std::cout << "\nval = " << m_cdf_ptr->approximate_mean();
      std::cout << "\nval = " << m_cdf_ptr->approximate_stddev();
      std::cout << "\nval = " << m_cdf_ptr->quantile(0.02); // Percentile values
      std::cout << "\nval = " << m_cdf_ptr->quantile(0.98) << std::endl;*/
      m_mutex.unlock();
    }
  }; // End class ParallelCdfFunctor
  

  /// Compute the CDF of an image using multiple threads.
  /// - The CDF object must be "fresh" when passed to this function.
  /// - Consider initializing the CDFAccumulator object with a large buffer.
  template <class ViewT>
  void block_cdf_computation(ImageViewBase<ViewT> const& image,
                             math::CDFAccumulator<float> &cdf,
                             int      subsample_amt = 1,
                             Vector2i block_size    = Vector2i(256,256)) {
    // Set up the functor, then execute it in parallel.
    ParallelCdfFunctor<float> cdf_functor(&cdf, subsample_amt);

    // No need for a cache since each tile will be visited only once.
    block_op(image, cdf_functor, block_size);
  }


#include <vw/Image/Statistics.tcc>

}  // namespace vw

#endif // __VW_IMAGE_STATISTICS_H__
