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


#include <vw/HDR/LocalToneMap.h>

#include <vw/Image/Statistics.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/Filter.h>
#include <vw/Image/Algorithms.h>
#include <vw/Core/Functors.h>

#include <vector>
#include <iostream>

using namespace vw;
using namespace vw::hdr;

// ********************************************************************
//  Ashikhmin operator
// ********************************************************************
const unsigned ASH_MAX_KERNEL = 10;

ImageView<double> ashikhmin_world_adaptation_luminance(ImageView<double> L_w, double threshold) {
  typedef ImageView<double> Map;

  std::vector<Map> L_w_blur(ASH_MAX_KERNEL * 2);
  std::vector<Map> V(ASH_MAX_KERNEL);

  vw_out() << "Computing L_w_blur\n";
  for ( unsigned s = 1; s <= ASH_MAX_KERNEL * 2; ++s ) {
    if ((s < ASH_MAX_KERNEL) || (s % 2 == 0)) {
      vw_out() << "L_w_blur[" << (s) << "]\n";
      L_w_blur[s-1] = gaussian_filter(L_w, 1.0, 1.0, s, s);
    }
  }

  vw_out() << "Computing Vs\n";
  for ( unsigned s = 1; s <= ASH_MAX_KERNEL; ++s ) {
    V[s-1] = abs((L_w_blur[s-1] - L_w_blur[2*s - 1]) / (L_w_blur[s-1] + 0.0001));
  }

  vw_out() << "Computing L_wa\n";
  Map L_wa(L_w.cols(), L_w.rows());
  for ( int32 x = 0; x < L_wa.cols(); ++x ) {
    for ( int32 y = 0; y < L_wa.rows(); ++y ) {
      unsigned s_t = 1;
      while ((s_t < ASH_MAX_KERNEL) && (V[s_t - 1](x,y) <= threshold)) {
        ++s_t;
      }
      L_wa(x,y) = L_w_blur[s_t - 1](x,y);
    }
  }

  return L_wa;
}

struct AshikhminCompressiveFunctor : ReturnFixedType<double> {
private:
  double C_L_wmin, k;

public:
  AshikhminCompressiveFunctor(double L_wmin, double L_wmax, double L_dmax = 1.0) {
    C_L_wmin = C(L_wmin);
    k = L_dmax / (C(L_wmax) - C_L_wmin);
    vw_out() << "C(L_wmin) = " << C_L_wmin << "\n";
    vw_out() << "k = " << k << "\n";
  }

  double C(double L) const {
    if (L < 0.0034) return (L / 0.0014);
    if (L < 1.0) return (2.4483 + log10(L/0.0034) / 0.4027);
    if (L < 7.2444) return (16.5630 + (L-1) / 0.4027);
    return (32.0693 + log10(L/7.2444) / 0.0556);
  }

  double operator() (double L_wa) const {
    return k * (C(L_wa) - C_L_wmin);
  }
};

ImageView<PixelRGB<double> > vw::hdr::ashikhmin_tone_map(ImageView<PixelRGB<double> > hdr_image, double threshold) {
  ImageView<PixelGray<double> > gray = pixel_cast<PixelGray<double> >(hdr_image);
  ImageView<double> L_w = channels_to_planes(gray);
  //  write_image("lw.jpg", L_w);
  ImageView<PixelRGB<double> > out_image = copy(hdr_image) / L_w;

  // Compute world adaptation luminance
  ImageView<double> L_wa = ashikhmin_world_adaptation_luminance(L_w, threshold);
  //  write_image("lwa.jpg", L_wa);

  // Compute display luminances
  vw_out() << "Computing L_wmin and L_wmax\n";
  double L_wmin, L_wmax;
  min_max_channel_values(L_w, L_wmin, L_wmax);
  vw_out() << "Computing F(L_wa)\n";
  AshikhminCompressiveFunctor F(L_wmin, L_wmax);
  ImageView<double> F_L_wa = per_pixel_filter(L_wa, F);
  //  write_image("flwa.jpg", F_L_wa);
  vw_out() << "Computing display luminances\n";
  ImageView<double> L_d = F_L_wa * L_w / L_wa;
  //  write_image("ld.jpg", L_d);

  // Recombine luminance values into color image
  return normalize(out_image * L_d);
}
