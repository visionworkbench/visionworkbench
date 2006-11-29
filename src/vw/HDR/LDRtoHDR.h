/* LDRtoHDR.h
 * 
 * Functions for stitching multiple LDR images of the same
 * scene into one HDR image.
 * 
 * Copyright 2006 NASA Ames Research Center.
 * Created 22 June 2006 by mihelich
 */ 

#ifndef __LDRTOHDR_H__
#define __LDRTOHDR_H__

#include <vw/Image/Filter.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/Manipulation.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/Algorithms.h>
#include <vw/FileIO.h>
#include <vw/Core/Exception.h>
#include <vw/HDR/CameraCurve.h>

#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <ctime>
#include <math.h>


namespace vw { 
namespace hdr {

  /// Stitches the set of LDR images into one HDR image.  The camera
  /// response curve is estimated from the provided images, however
  /// you must supply a vector of brightness values corresponding to
  /// the sensitivity of the exposure of each image.  (only the ratio
  /// of the brightness values is important here, although you are
  /// best off working in physical units of lumens, etc.)  Brightness
  /// is a function of camera ISO (sensitivity or sensor gain),
  /// shutter speed, and aperture.
  ///
  /// This code roughly implements Mitsunaga nod Nayar technique for
  /// reconstructing HDR images from a bracketed set of LDR
  /// exposures. See:
  ///
  /// T. Mitsunaga and S.K.Nayar. "Radiometric Self Calibration,"
  /// in Proceedings of the IEEE Conference on Computer Vision and
  /// Pattern Recoognition, Fort Collins, CO. June 1999.
  ///
  template <class PixelT>
  ImageView<PixelT> process_ldr_images(std::vector<ImageView<PixelT> > const &images, std::vector<Vector<double> > &ret_curves, std::vector<double> brightness_values) {
    typedef typename PixelChannelType<PixelT>::type channel_type;
    
    int n_channels = PixelNumChannels<PixelT>::value;
    std::vector<Matrix<double> > pairs(n_channels);
    std::vector<Vector<double> > curves(n_channels);

    // Sample each image channel
    for (int i = 0; i < n_channels; i++) {
      pairs[i] = generate_ldr_intensity_pairs(images, VW_HDR_NUM_PAIRS, i, brightness_values);
    }

    // Compute camera response curve for each channel. See CameraCurve.h.
    for (int i = 0; i < n_channels; i++) {
      estimate_camera_curve(pairs[i], curves[i], VW_HDR_RESPONSE_POLYNOMIAL_ORDER);
    }

    ret_curves = curves;

    // Stitch LDR image set into HDR image
    ImageView<PixelT> hdr_image = ldr_to_hdr(images, curves, brightness_values);
    return hdr_image;
  }

  template <class PixelT>
  ImageView<PixelT> ldr_to_hdr(std::vector<ImageView<PixelT> > images, std::vector<Vector<double> > const &curves, std::vector<double> const &brightness_values) {
    int width = images[0].cols();
    int height = images[0].rows();
    int n_channels = PixelNumChannels<PixelT>::value;
    ImageView<PixelT> hdr_image(width, height);
    ImageView<PixelT> image_sum(width, height);
    ImageView<double> weight_sum(width, height);
    
    // Convert all pixel values to scaled radiance values
    for (int i = 0; i < images.size(); i++) {
      psi(images[i], curves);
    }
    
    // Bring all images into same domain and average pixels across images using
    // a weighting function that favors pixels in middle of dynamic range
    double min_bv = *(min_element(brightness_values.begin(), brightness_values.end()));
    for (int i = 0; i < images.size(); i++) {
      ImageView<PixelGray<double> > gray = images[i];
      ImageView<double> weight = 1.0 - 2.0 * abs(0.5 - channels_to_planes(gray));
      double ratio = pow(2.0, (brightness_values[i] - min_bv) * 0.5);
      //      std::cout << "Ratio: " << ratio << "\n";
      hdr_image += weight * images[i] * ratio;
      weight_sum += weight;
    }
    
    // Divide by sum of weights, compensating for very small weightings
    std::vector<Vector<int,2> > sat;
    for (int col = 0; col < width; col++) {
      for (int row = 0; row < height; row++) {
        if (weight_sum(col, row) < 0.2) {
          sat.push_back(Vector<int,2>(col, row));
          weight_sum(col, row) = 1.0;
        }
      }
    }
    hdr_image /= weight_sum;
    ImageView<PixelGray<double> > gray = hdr_image;
    ImageView<PixelGray<double> > bright_or_dark = images[images.size() / 2];
    double min, max;
    min_max_channel_values(bright_or_dark, min, max);
    double threshold = (min + max) / 2;
    min_max_channel_values(gray, min, max);
    
    // Assign very bright/dark pixels
    for (int i = 0; i < sat.size(); i++) {
      int col = sat[i].x();
      int row = sat[i].y();
      double value;
      if (bright_or_dark(col, row)[0] > threshold) {
        value = max;
      } else {
        value = min;
      }
      for (int channel = 0; channel < n_channels; channel++) {
        hdr_image(col, row)[channel] = value;
      }
    }
    
    return hdr_image;
  }


  /// Stitches the set of LDR images into one HDR image.  Use this
  /// variant of the function if you only know the constant exposure
  /// rations (spacing between consecutive exposures).  
  template <class PixelT>
  ImageView<PixelT> process_ldr_images(std::vector<ImageView<PixelT> > const &images, std::vector<Vector<double> > &ret_curves, double ratio = VW_HDR_DEFAULT_FSTOP_RATIO) {
    std::vector<double> brightness_values(images.size());
    for (int i=0; i < brightness_values.size();++i) {
      brightness_values[i] = pow(ratio, i);
    }
    return process_ldr_images(images, ret_curves, brightness_values);      
  }
  
}} // vw::HDR 
        
#endif  // __LDRtoHDR_H__
