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
#ifndef __VW_HDR_POLYNOMIALCAMERACURVE_H__
#define __VW_HDR_POLYNOMIALCAMERACURVE_H__

#include <vw/Core/Exception.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>
#include <vw/Math/LevenbergMarquardt.h>

// Approximate response curves with polynomials of order 8
const int VW_HDR_DEFAULT_POLYNOMIAL_ORDER = 13;

// Number of LDR intensity pairs to sample
const int VW_HDR_NUM_PAIRS = 3000;

namespace vw { 
namespace hdr {

  // -----------------------------------------------------------------------
  // Camera Curve Generation
  // -----------------------------------------------------------------------
  //
  // Given a long list of corresponding pixels and their relative
  // exposure values, these functions can be used to estimate the
  // Camera's response function in either the forward (luminance to
  // pixel value) or inverse (pixel value to luminance) directions.


  /// Approximate estimation of the camera response curve as a 
  /// polynomial.  This appoarch uses a random sampling of pixel 
  /// values from mutliple exposures of the same image.  The 
  /// solution can be computed in closed form, and a future version
  /// of this function should probably rely on the closed form sol'n,
  /// which can be computed much faster.
  /// 
  /// The input to this function, 'pairs', is an Nx3 matrix where each
  /// row contains a pixel value from image 1, the corresponding pixel 
  /// value from image 2, and the ratio of exposure between these two
  /// images.
  /// 
  /// So far, we have found that 1 photographic stop (a factor of 2 
  /// in shutter speed or 1 f-stop) equates to a difference of sqrt(2) 
  /// in the amount of light that falls on the sensor.  Therefore, for 
  /// photos that are seperated by 1 f-stop, the ratios in the third 
  /// column of pixel_pairs should be powers of sqrt(2).
  vw::Vector<double> estimate_polynomial_camera_curve(vw::Matrix<double> &pixel_pairs,
                                                      unsigned int polynomial_order);
      

  /// Computes the camera curve for LDR images of the same scene. This
  /// approach uses a polynomial approximation of the camera response
  /// curve for each channel (interactions between channels are ignored)
  /// developed by Mitsunaga and Nayar.
  ///
  /// The input to this function, 'images', is a std::vector of
  /// images sorted from darkest to brightest. Consecutive images
  /// should have a constant exposure ratio (if separated by one
  /// f-stop, the default value sqrt(2) may be used). Each
  /// image should be single-plane with channels stored as a
  /// floating point type.
  ///
  /// sample_region_size is given in units of pixels, and it
  /// determines the size of the neighborhood that is averaged when
  /// picking corresponding points samples between LDR images.
  template <class ViewT>
  std::vector<Vector<double> > polynomial_camera_curves(std::vector<ViewT> const &images, 
                                                        std::vector<double> brightness_values,
                                                        int sample_region_size = 1,
                                                        double polynomial_order = VW_HDR_DEFAULT_POLYNOMIAL_ORDER) {
    
    int32 n_channels = PixelNumChannels<typename ViewT::pixel_type>::value;
    
    // Sample each image channel
    std::vector<Matrix<double> > pairs(n_channels);
    for ( int32 i = 0; i < n_channels; ++i ) {
      pairs[i] = generate_ldr_intensity_pairs(images, brightness_values, VW_HDR_NUM_PAIRS, i, sample_region_size);
    }

    // Compute camera response curve for each channel. 
    std::vector<Vector<double> > curves(n_channels);
    for ( int32 i = 0; i < n_channels; ++i ) {
      curves[i] = estimate_polynomial_camera_curve(pairs[i], polynomial_order);
    }

    return curves;
  }

  
  /// Evaluate a polynomial given the set of coefficients provided in
  /// theta as the N coefficients, and the first coefficient is the
  /// constant term.
  inline double psi(double x, Vector<double> const &theta) {    

    // Compute the main term:    (sum(k = 0, N-1, theta_k * x^k)) 
    double L = 0;
    for (unsigned int i = 0; i < theta.size(); i++) {
      L += theta(i) * pow(x, i); 
    }
    return L;
  }

  /// Evaluate a polynomial given the set of coefficients provided in
  /// theta as the N coefficients, and the first coefficient is the
  /// constant term.
  template <class PixelT>
  inline typename CompoundChannelCast<PixelT, double>::type 
  psi(PixelT x, std::vector<Vector<double> > const &curves) {    
    typedef typename CompoundChannelCast<PixelT, double>::type pixel_type;
    
    pixel_type return_val;

    for (int c = 0; c < PixelNumChannels<pixel_type>::value; ++c) {
      // Compute the main term:    (sum(k = 0, N-1, theta_k * x^k)) 
      double L = 0;
      for (unsigned int i = 0; i < curves[c].size(); i++) {
        L += curves[c](i) * pow(x[c], i); 
      }
      return_val[c] = L;
    }
    return return_val;
  }

  /// Converts each pixel value in the image to scaled illuminance
  /// values based on a set of polynomial response curves, one
  /// curve for each image channel.
  template <class ViewT>
  class PolynomialLuminanceView : public ImageViewBase<PolynomialLuminanceView<ViewT> > {
    
    ViewT m_view;
    std::vector<vw::Vector<double> > m_curves;
    float m_brightness_val;
    
  public:

    typedef typename ViewT::pixel_type src_pixel_type;
    typedef typename CompoundChannelCast<src_pixel_type, double>::type pixel_type;
    typedef pixel_type result_type;
    typedef ProceduralPixelAccessor<PolynomialLuminanceView> pixel_accessor;

    PolynomialLuminanceView( ImageViewBase<ViewT> const& view, 
                   std::vector<vw::Vector<double> > const& curves,
                   double brightness_val) : 
      m_view(view.impl()), m_curves(curves), m_brightness_val(brightness_val) {}

    inline int32 cols() const { return m_view.cols(); }
    inline int32 rows() const { return m_view.rows(); }
    inline int32 planes() const { return m_view.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(*this); }

    inline pixel_type operator()( int i, int j, int p=0 ) const { 
      pixel_type pix;

      // The psi() function returns a value between 0.0 and 1.0.  The
      // supplied brightness value is the luminance corresponding to a
      // gray pixel ( psi(x) = 0.5 ), so we must scale the brightness
      // by a factor of 2 to achieve the correct correspondence.
      return 2 * m_brightness_val * psi(m_view(i,j,p), m_curves);
    }

    /// \cond INTERNAL
    typedef PolynomialLuminanceView<typename ViewT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const { return prerasterize_type( m_view.prerasterize(bbox), m_curves, m_brightness_val ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };
  
  /// Given a set polynomial coefficients representing of forward or
  /// inverse camera curves, produce the inverted version of the
  /// curves.
  std::vector<vw::Vector<double> > invert_polynomial(std::vector<Vector<double> > &curves,
                                                     double lower_bound = 0.0, double upper_bound = 1.0);


  /// Read and write the coefficients for the Nth order polynomial
  /// used to approximate the camera curves.
  void write_polynomial_curve_coefficients(std::string const& generated_curves_file, 
                                std::vector<vw::Vector<double> > const &curves);
  void read_polynomial_curve_coefficients(std::vector<vw::Vector<double> > &curves, 
                               std::string const& curves_input_file);

  /// Read and write the actual camera curves themselves in a
  /// tabulated data format.
  void write_polynomial_curves(std::string const& curves_file, 
                    std::vector<vw::Vector<double> > const &curves);

}} // namespace vw::HDR 


#endif  // __VW_HDR_POLYNOMIALCAMERACURVE_H__
