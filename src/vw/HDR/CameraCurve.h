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

/// \file CameraCurve.h
/// 
/// Functions for deducing the camera response curve by comparing 
/// relative brightness values from several images of the same scene.
///
#ifndef __VW_HDR_CAMERACURVE_H__
#define __VW_HDR_CAMERACURVE_H__

#include <vw/Core/Exception.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Optimization.h>

#define VW_HDR_GET_RAND(max) int((max) * double(rand()) / (RAND_MAX + 1.0))

// Normal ratio of one f-stop is sqrt(2)
const double VW_HDR_DEFAULT_FSTOP_RATIO = 1.414213562;

// Approximate response curves with polynomials of order 8
const int VW_HDR_RESPONSE_POLYNOMIAL_ORDER = 8;

// Number of LDR intensity pairs to sample
const int VW_HDR_NUM_PAIRS = 3000;

namespace vw { 
namespace hdr {

  void write_curves_file(std::string const& generated_curves_file, 
                         std::vector<vw::Vector<double> > const &curves);

  void read_curves_file(std::vector<vw::Vector<double> > &curves, 
                        std::string const& curves_input_file);

  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   * Approximate estimation of the camera response curve as a 
   * polynomial.  This appoarch uses a random sampling of pixel 
   * values from mutliple exposures of the same image.  The 
   * solution can be computed in closed form, and a future version
   * of this function should probably rely on the closed for sol'n,
   * which can be computed much faster.
   * 
   * The input to this function, 'pairs', is an Nx3 matrix where each
   * row contains a pixel value from image 1, the corresponding pixel 
   * value from image 2, and the ratio of exposure between these two
   * images.
   * 
   * So far, we have found that 1 photographic stop (a factor of 2 
   * in shutter speed or 1 f-stop) equates to a difference of sqrt(2) 
   * in the amount of light that falls on the sensor.  Therefore, for 
   * photos that are seperated by 1 f-stop, the ratios in the third 
   * column of pixel_pairs should be powers of sqrt(2).
   * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  void estimate_camera_curve(vw::Matrix<double> &pixel_pairs,
                             vw::Vector<double> &poly,
                             unsigned int polynomial_order);
    
  void invert_curve(vw::Vector<double> &curve, vw::Vector<double> &inverse,
                      unsigned int polynomial_order,
                      double lower_bound = 0.0,
                      double upper_bound = 1.0);


    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     * Samples an image channel at the specified indices using a
     * sample region of size kernel_size. kernel_size should be odd.
     * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    template <class PixelT>
    typename PixelChannelType<PixelT>::type sample_image(ImageView<PixelT> image, int x, int y, int channel, int kernel_size) {
      typedef typename PixelChannelType<PixelT>::type channel_type;
      
      int halfsize = kernel_size / 2;
      if ((x - halfsize < 0) || (x + halfsize >= image.cols()) ||
          (y - halfsize < 0) || (y + halfsize >= image.rows())) {
        return 0;
      }
      channel_type average = 0;
      for (int col = x - halfsize; col <= x + halfsize; col++) {
        for (int row = y - halfsize; row <= y + halfsize; row++) {
          average += select_channel(image, channel)(col, row);
        }
      }
      average /= kernel_size * kernel_size;
      return average;
    }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   * Generates an Nx3 matrix where each row contains a channel value
   * from one LDR image, the corresponding pixel value from a second
   * LDR image, and the ratio of exposure between these two images.
   * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  template <class PixelT>
  Matrix<typename PixelChannelType<PixelT>::type> generate_ldr_intensity_pairs(std::vector<ImageView<PixelT> > const &images, int num_pairs, int channel, std::vector<double> const &brightness_values) {
    typedef typename PixelChannelType<PixelT>::type channel_type;

    Matrix<channel_type> pair_list(num_pairs, 3);
    int n_images = images.size();
    VW_ASSERT(n_images > 1, ArgumentErr() << "Need at least two images.");
    int height = images[0].rows();
    int width = images[0].cols();
    int n_channels = PixelNumChannels<PixelT>::value;
    VW_ASSERT((channel >= 0) && (channel < n_channels), ArgumentErr() << "No such channel.");
    srand(time(0)); // Initialize random number generator
    int i = 0;
    while (i < num_pairs) {
      // Generate random indices for two images
      int rand_x = VW_HDR_GET_RAND(width);
      int rand_y = VW_HDR_GET_RAND(height);

      // Pick two distinct images to sample from
      int id1 = VW_HDR_GET_RAND(n_images);
      int id2;
      while (true) {
        id2 = VW_HDR_GET_RAND(n_images);
        if (id1 != id2) break;
      }

      // Sample both images at those indices
      int kernel_size = 1;
      channel_type I_1 = sample_image(images[id1], rand_x, rand_y, channel, kernel_size);
      channel_type I_2 = sample_image(images[id2], rand_x, rand_y, channel, kernel_size);

      // Add the sample pair to the list along with the exposure ratio
      if ((I_1 > 0.01) && (I_1 < 0.99) && (I_2 > 0.1) && (I_2 < 0.9)) {
        pair_list(i, 0) = I_1;
        pair_list(i, 1) = I_2;
        pair_list(i, 2) = pow(2.0, (brightness_values[id2] - brightness_values[id1]) * 0.5);
        i++;
      }
    }

    return pair_list;
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   * Generates an Nx3 matrix where each row contains a channel value
   * from one LDR image, the corresponding pixel value from a second
   * LDR image, and the ratio of exposure between these two images.
   * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  template <class PixelT>
  Matrix<typename PixelChannelType<PixelT>::type> generate_ldr_intensity_pairs(std::vector<ImageView<PixelT> > const &images, int num_pairs, int channel, double ev_ratio) {
    std::vector<double> brightness_values(images.size());
    for (int i=0; i < brightness_values.size();++i) {
      brightness_values[i] = pow(ev_ratio, i);
    }
    return generate_ldr_intensity_pairs(images, num_pairs, channel, brightness_values);
  }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     * Computes the camera curve for LDR images of the same scene. This
     * approach uses a polynomial approximation of the camera response
     * curve for each channel (interactions between channels are ignored)
     * developed by Mitsunaga and Nayar.
     *
     * The input to this function, 'images', is a std::vector of
     * images sorted from darkest to brightest. Consecutive images
     * should have a constant exposure ratio (if separated by one
     * f-stop, the default value sqrt(2) may be used). Each
     * image should be single-plane with channels stored as a
     * floating point type.
     * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    template <class PixelT>
    std::vector<Vector<double> > camera_curves(std::vector<ImageView<PixelT> > const &images, double ratio = VW_HDR_DEFAULT_FSTOP_RATIO) {
      typedef typename PixelChannelType<PixelT>::type channel_type;
      
      int n_channels = PixelNumChannels<PixelT>::value;
      std::vector<Matrix<double> > pairs(n_channels);
      std::vector<Vector<double> > curves(n_channels);
      
      // Sample each image channel
      for (int i = 0; i < n_channels; i++) {
        pairs[i] = generate_ldr_intensity_pairs(images, VW_HDR_NUM_PAIRS, i, ratio);
      }
      
      // Compute camera response curve for each channel. 
      for (int i = 0; i < n_channels; i++) {
        estimate_camera_curve(pairs[i], curves[i], VW_HDR_RESPONSE_POLYNOMIAL_ORDER);
      }
      
      return curves;
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     * Converts each pixel value in the image to scaled radiance
     * values based on a set of polynomial response curves, one
     * curve for each image channel.
     * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    template <class PixelT>
    void psi(ImageView<PixelT> &image, std::vector<Vector<double> > const &curves) {
      int width = image.cols();
      int height = image.rows();
      int n_channels = PixelNumChannels<PixelT>::value;
  
      for (int channel = 0; channel < n_channels; channel++) {
        Vector<double> theta = curves[channel];
        for (int col = 0; col < width; col++) {
          for (int row = 0; row < height; row++) {
            // Evaluate the polynomial
            double x = image(col,row)[channel];
            double f_x = 0.0;
            double pow_x = 1.0;
            for (int e = 0; e < theta.size(); e++) {
              f_x += theta(e) * pow_x;
              pow_x *= x;
            }
            image(col, row)[channel] = f_x;
          }
        }
      }
    } 


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   * Cost function object for the nonlinear least squares optimizer
   *
   * The VW least squares optimizer will optimize a function that 
   * is a CRTP subclass of LeastSquaresModelBase.  Here we specialize
   * for our problem of solving for the camera curve by implementing
   * the methods operator() and jacobian().
   * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/  
  class CameraCurveCostFn : public vw::math::LeastSquaresModelBase<CameraCurveCostFn> {

    Matrix<double> m_pairs;
    
  public:
    typedef Vector<double> domain_type;
    typedef Vector<double> result_type;
    typedef Matrix<double> jacobian_type;
    
    // Constructor 
    CameraCurveCostFn(Matrix<double> pixel_pairs): m_pairs(pixel_pairs) {}
    
    result_type operator()( domain_type const& x ) const;

    // Compute the actual jacobian.  This method overrides the
    // numerical differentiation in the base class
    Matrix<double> jacobian( domain_type const& x ) const;
  };
    
  class InverseCostFn : public vw::math::LeastSquaresModelBase<InverseCostFn> {
    Matrix<double> m_pairs;
  public:
    typedef Vector<double> domain_type;
    typedef Vector<double> result_type;
    typedef Matrix<double> jacobian_type;

    InverseCostFn(Matrix<double> pairs) : m_pairs(pairs) {}

    result_type operator()( domain_type const& x ) const;
  };
  
}} // namespace vw::HDR 


#endif  // __VW_HDR_CAMERACURVE_H__
