// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file CameraCurve.h
///
/// Functions for deducing the camera response curve by comparing
/// relative brightness values from several images of the same scene.
///
#ifndef __VW_HDR_CAMERACURVE_H__
#define __VW_HDR_CAMERACURVE_H__

#include <vw/Image/PerPixelViews.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Image/Statistics.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>

// Number of LDR intensity pairs to sample
const int VW_HDR_DEFAULT_NUM_PIXEL_SAMPLES = 300;

namespace vw {
namespace hdr {
  namespace detail {
    inline uint32 dice(uint32 max_) {
      return static_cast<uint32>(max_ * (double(rand()) / double(RAND_MAX)));
    }
  }


  /// Generate a set of brightness values based on a known ratio of
  /// exposures between images.  Note that this brightness value will
  /// only be correct in a relative sense.
  std::vector<double> brightness_values_from_exposure_ratio(double exposure_ratio, int size);

  /// Generate a list of brightness values from EXIF information
  /// stored in a list of files.
  std::vector<double> brightness_values_from_exif(std::vector<std::string> const& filenames);


  // -----------------------------------------------------------------------
  // Pixel Pair Extraction
  // -----------------------------------------------------------------------
  //
  // These routines are used internally by the HDR module to sample
  // LDR images for pairs of corresponding pixels with different
  // exposures.  These pixel pairs are used by the camera curve
  // estimation functions (later in this file) to estimate the
  // camera's response function.
  //

  /// Samples an image channel at the specified indices using a
  /// sample region of size kernel_size. kernel_size should be odd.
  template <class ViewT>
  typename PixelChannelType<typename ViewT::pixel_type>::type sample_image(ImageViewBase<ViewT> const& image,
                                                                           int x, int y, int channel, int kernel_size) {

    typedef typename PixelChannelType<typename ViewT::pixel_type>::type channel_type;
    EdgeExtensionView<ViewT, ConstantEdgeExtension> edge_extended_image(image.impl(), ConstantEdgeExtension());

    int halfsize = kernel_size / 2;
    double average = 0;

    for ( int col = x - halfsize; col <= x + halfsize; ++col ) {
      for ( int row = y - halfsize; row <= y + halfsize; ++row ) {
        average += select_channel(edge_extended_image, channel)(col, row);
      }
    }
    average /= (kernel_size * kernel_size);
    return channel_type(average);
  }

  /// Generates an Nx3 matrix where each row contains a channel value
  /// from one LDR image, the corresponding pixel value from a second
  /// LDR image, and the ratio of exposure between these two images.
  template <class ViewT>
  Matrix<typename PixelChannelType<typename ViewT::pixel_type>::type> generate_ldr_intensity_pairs(std::vector<ViewT> const &images,
                                                                                                   std::vector<double> const &brightness_values,
                                                                                                   int num_pairs, uint32 channel,
                                                                                                   int kernel_size = 1) {

    typedef typename PixelChannelType<typename ViewT::pixel_type>::type channel_type;
    uint32 n_channels = PixelNumChannels<typename ViewT::pixel_type>::value;

    // Error checking
    VW_ASSERT(images.size() > 1, ArgumentErr() << "Need at least two images.");
    VW_ASSERT(channel < n_channels, ArgumentErr() << "No such channel.");

    Matrix<channel_type> pair_list(num_pairs, 3);
    int height = images[0].impl().rows();
    int width = images[0].impl().cols();

    srand(time(0)); // Initialize random number generator
    int i = 0;
    while (i < num_pairs) {
      // Generate random indices for two images
      int rand_x = detail::dice(width);
      int rand_y = detail::dice(height);

      // Pick two distinct images to sample from
      int id1 = detail::dice(images.size());
      int id2;
      while (true) {
        id2 = detail::dice(images.size());
        if (id1 != id2) break;
      }

      // Sample both images at those indices
      channel_type I_1 = sample_image(images[id1].impl(), rand_x, rand_y, channel, kernel_size);
      channel_type I_2 = sample_image(images[id2].impl(), rand_x, rand_y, channel, kernel_size);

      // We would prefer to avoid points that are very close to
      // saturated.  If the values are within 90% of the dynamic
      // range, add the sample pair to the list along with the
      // exposure ratio
      if ((I_1 > 0.01) && (I_1 < 0.99) && (I_2 > 0.01) && (I_2 < 0.99)) {

        // FOR DEBUGGING:
        //                vw_out() << "Adding ["<< id1 << " " << id2 <<"] " << I_1 << "   " << I_2 << "   " << brightness_values[id1] << "    " << brightness_values[id2] << "\n";
        //                vw_out() << "      " << pow(2.0, (brightness_values[id2] - brightness_values[id1]) * 0.5) << "  " << brightness_values[id2]/brightness_values[id1] << "\n";

        pair_list(i, 0) = I_1;
        pair_list(i, 1) = I_2;
        pair_list(i, 2) = brightness_values[id2]/brightness_values[id1];
        ++i;
      }
    }

    return pair_list;
  }

  /// Generates an Nx3 matrix where each row contains a channel value
  /// from one LDR image, the corresponding pixel value from a second
  /// LDR image, and the ratio of exposure between these two images.
  template <class ViewT>
  Matrix<typename PixelChannelType<typename ViewT::pixel_type>::type> sample_ldr_images(std::vector<ViewT> const &images,
                                                                                        std::vector<double> const &/*brightness_values*/,
                                                                                        int num_pairs, int channel,
                                                                                        int kernel_size = 1) {

    typedef typename PixelChannelType<typename ViewT::pixel_type>::type channel_type;
    uint32 n_channels = PixelNumChannels<typename ViewT::pixel_type>::value;

    // Error checking
    VW_ASSERT(images.size() > 1, ArgumentErr() << "Need at least two images.");
    VW_ASSERT((channel >= 0) && (channel < int(n_channels)), ArgumentErr() << "No such channel.");

    Matrix<channel_type> pair_list(num_pairs, images.size());
    int height = images[0].impl().rows();
    int width = images[0].impl().cols();

    srand(time(0)); // Initialize random number generator
    int i = 0;
    while (i < num_pairs) {
      // Generate random indices for two images
      int rand_x = detail::dice(width);
      int rand_y = detail::dice(height);

      for (unsigned j = 0; j < images.size(); ++j)
        pair_list(i,j) = sample_image(images[j].impl(), rand_x, rand_y, channel, kernel_size);
      ++i;
    }

    return pair_list;
  }

  /// This is useful for debugging if you want to save out the
  /// collection of corresponding pixels.  You may want to us this to
  /// develop new camera response curve estimators in some other
  /// program, such as MATLAB.
  template <class ElemT>
  void write_camera_pairs(std::string const& pairs_file,
                          vw::Matrix<ElemT> const &pairs) {
    FILE* output_file = fopen(pairs_file.c_str(), "w");
    if ( !output_file ) vw_throw( IOErr() << "write_camera_pairs: failed to open file for writing." );
    for (double i = 0; i < pairs.rows(); ++i) {
      for ( unsigned j = 0; j < pairs.cols(); ++j ) {
        fprintf(output_file, "%f ", pairs(i,j));
      }
      fprintf(output_file, "\n");
    }
    fclose(output_file);
  }

  // ---------------------------------------------------------------------
  // Camera Curve Estimation
  // ---------------------------------------------------------------------

  class CameraCurveFn {

    // The lookup table stores the log luminance values for each pixel
    // value.
    std::vector<Vector<double> > m_lookup_tables;

  public:
    CameraCurveFn(std::vector<Vector<double> > const& lookup_tables) :
      m_lookup_tables(lookup_tables) {}

    CameraCurveFn() {}

    // Returns the luminance value for a given pixel value.  Linearly
    // interpolates between points in the lookup table.  Pixel value
    // is expected to be a value in the range [0.0 1.0]
    double operator() (double pixel_val, unsigned channel) const {
      if (channel >= m_lookup_tables.size())
        vw_throw(ArgumentErr() << "CameraCurveFn: unknown lookup table.");

      double scaled_pixel_val = pixel_val*(m_lookup_tables[channel].size()-1);
      int idx1 = int(floor(scaled_pixel_val));
      int idx2 = int(ceil(scaled_pixel_val));
      double val1 = m_lookup_tables[channel][idx1];
      double val2 = m_lookup_tables[channel][idx2];

      double frac = scaled_pixel_val - idx1;
      return exp(val1 + (val2-val1) * frac);
    }

    template <class PixelT>
    typename CompoundChannelCast<PixelT, double>::type operator() (PixelT pixel_val) const {
      typedef typename CompoundChannelCast<PixelT, double>::type pixel_type;

      if (CompoundNumChannels<PixelT>::value != int32(this->num_channels()))
        vw_throw(ArgumentErr() << "CameraCurveFn: pixel does not have the same number of channels as there are curves.");

      pixel_type result;
      for (int c = 0; c < CompoundNumChannels<PixelT>::value; ++c) {
        result[c] = this->operator()(double(pixel_val[c]), c);
      }
      return result;
    }

    unsigned num_channels() const { return m_lookup_tables.size(); }

    Vector<double> const& lookup_table(int channel) const {
      if (channel < 0 || (size_t)channel >= m_lookup_tables.size())
        vw_throw(ArgumentErr() << "CameraCurveFn: unknown lookup table.");

      return m_lookup_tables[channel];
    }

  };

  /// Returns a lookup table, indexed by pixel values, that contains
  /// their corresponding log luminance values.
  Vector<double> estimate_camera_curve(vw::Matrix<double> const& pixels,
                                      std::vector<double> const& brightness_values);

  /// Computes the camera curve for LDR images of the same scene.
  ///
  /// The input to this function, 'images', is a std::vector of images
  /// sorted from darkest to brightest. Each image should be
  /// single-plane with channels stored as a floating point channel
  /// type with values in the range [0.0 1.0].
  ///
  /// sample_region_size is given in units of pixels, and it
  /// determines the size of the neighborhood that is averaged when
  /// picking corresponding points samples between LDR images.
  template <class ViewT>
  CameraCurveFn camera_curves(std::vector<ViewT> const &images,
                              std::vector<double> brightness_values,
                              int sample_region_size = 1) {

    int32 n_channels = PixelNumChannels<typename ViewT::pixel_type>::value;

    // Sample each image channel
    std::vector<vw::Matrix<double> > pixels(n_channels);
    for ( int32 i = 0; i < n_channels; ++i ) {
      pixels[i] = sample_ldr_images(images, brightness_values, VW_HDR_DEFAULT_NUM_PIXEL_SAMPLES, i, sample_region_size);
    }

    // Compute camera response curve for each channel.
    std::vector<Vector<double> > lookup_tables(n_channels);
    for ( int32 i = 0; i < n_channels; ++i ) {
      lookup_tables[i] = estimate_camera_curve(pixels[i], brightness_values);
    }

    return CameraCurveFn(lookup_tables);
  }

  /// Write the camera curve values in a tabulated format on disk.
  void write_curves(std::string const& curves_file, CameraCurveFn const &curves);

  /// Read the camera curve values from a tabulated format on disk.
  CameraCurveFn read_curves(std::string const& curves_file);

  /// A pixel casting functor, used by \ref pixel_cast().
  template <class PixelT>
  class LuminanceFunc : public ReturnFixedType<typename CompoundChannelCast<PixelT,double>::type> {
    CameraCurveFn m_curves;
    double m_brightness_val;

  public:
    LuminanceFunc(CameraCurveFn const& curves, double brightness_val) :
      m_curves(curves), m_brightness_val(brightness_val) {}

    inline typename CompoundChannelCast<PixelT,double>::type operator()( PixelT pixel ) const {
      return m_brightness_val*m_curves(pixel);
    }
  };

  /// Create a new image view that contains luminance values rather than pixel values.
  template <class ImageT>
  inline UnaryPerPixelView<ImageT,LuminanceFunc<typename ImageT::pixel_type> > luminance_image( ImageViewBase<ImageT> const& image,
                                                                                                CameraCurveFn const& curves,
                                                                                                double brightness_val) {
    return UnaryPerPixelView<ImageT,LuminanceFunc<typename ImageT::pixel_type> >( image.impl(), LuminanceFunc<typename ImageT::pixel_type>(curves, brightness_val) );
  }

}} // namespace vw::HDR


#endif  // __VW_HDR_CAMERACURVE_H__
