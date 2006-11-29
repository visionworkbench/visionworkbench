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

 #include <vw/HDR/GlobalToneMap.h>

 #include <vw/Image/Statistics.h>
 #include <vw/Image/ImageMath.h>
 #include <vw/Image/Filter.h>
 #include <vw/Image/Algorithms.h>

 const int DRAGO_MIN_BASE = 2;
 const int DRAGO_MAX_BASE = 10;
 const float DRAGO_DEFAULT_BIAS = 0.85;

class DragoFunctor: public vw::UnaryReturnSameType {
 private:
   double L_wmax, power, c;
   int b_min, b_diff;
 public:
   DragoFunctor(double L_dmax, double L_wmax, double bias, int b_min, int b_max) : L_wmax(L_wmax), b_min(b_min) {
     b_diff = b_max - b_min;
     c = L_dmax * 0.01 * log(b_max) / log(L_wmax + 1.0);
     power = log(bias) / log(0.5);
   }

   double operator() (double L_w) const {
     return c * log(L_w + 1.0) / log(b_min + b_diff * pow(L_w / L_wmax, power));
   }
 };

vw::ImageView<vw::PixelRGB<double> > vw::hdr::drago_tone_map(ImageView<PixelRGB<double> > hdr_image, 
                                                             double bias, 
                                                             double exposure_factor, 
                                                             double L_dmax) {

  ImageView<double> L_map = pixel_cast<PixelGray<double> >(hdr_image);
  ImageView<PixelRGB<double> > out_image = copy(hdr_image) / L_map;
  
  // Compute world adaptation luminance
  double L_wa = exp(mean_channel_value(log(L_map + 0.00001)));

  // Perform overall brightness correction
  L_wa /= pow(1 + bias - DRAGO_DEFAULT_BIAS, 5.0);
  
  // Scale luminances
  L_map *= exposure_factor / L_wa;
  
  // Compute display luminances
  double L_wmax = max_channel_value(L_map);
  DragoFunctor drago_tm(L_dmax, L_wmax, bias, DRAGO_MIN_BASE, DRAGO_MAX_BASE);
  ImageView<double> L_d_map = per_pixel_filter(L_map, drago_tm);
  
  // Recombine luminance values into color image
  return normalize(out_image * L_d_map);
}
