// __BEGIN_LICENSE__
//  Copyright (c) 2006-2026, United States Government as represented by the
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

/// \file Detector.cc
///
/// Interest point detection functions
///

#include <vw/InterestPoint/Detector.h>
#include <vw/Image/Filter.h>
#include <vw/Image/Statistics.h>

namespace vw {
namespace ip {

float get_orientation(vw::ImageViewRef<float> const& x_grad,
                     vw::ImageViewRef<float> const& y_grad,
                     float i0, float j0, float sigma_ratio) {

  // The size, in pixels, of the image patch used to compute the orientation.
  static const int IP_ORIENTATION_WIDTH = 10;

  // Nominal feature support patch is WxW at the base scale, with
  // W = IP_ORIENTATION_HALF_WIDTH * 2 + 1, and
  // we multiply by sigma[k]/sigma[1] for other planes.
  //
  // Get bounds for scaled WxW window centered at (i,j) in plane k
  int halfwidth = (int)(IP_ORIENTATION_WIDTH/2*sigma_ratio + 0.5);
  int left  = int(roundf(i0 - halfwidth));
  int top   = int(roundf(j0 - halfwidth));

  // Compute (gaussian weight)*(edge magnitude) kernel
  ImageView<float> weight(IP_ORIENTATION_WIDTH,IP_ORIENTATION_WIDTH);
  make_gaussian_kernel_2d(weight, 6 * sigma_ratio, IP_ORIENTATION_WIDTH);

  // We must compute the average orientation in quadrature.
  double weight_sum = sum_of_pixel_values(weight);

  // Compute the gaussian weighted average x_gradient
  ImageView<float> weighted_grad = weight * crop(edge_extend(x_grad.impl()),left,top,IP_ORIENTATION_WIDTH,IP_ORIENTATION_WIDTH);
  double avg_x_grad = sum_of_pixel_values(weighted_grad) / weight_sum;

  // Compute the gaussian weighted average y_gradient
  weighted_grad = weight * crop(edge_extend(y_grad.impl()),left,top,IP_ORIENTATION_WIDTH,IP_ORIENTATION_WIDTH);
  double avg_y_grad = sum_of_pixel_values(weighted_grad) / weight_sum;

  return atan2(avg_y_grad,avg_x_grad);
}

}} // namespace vw::ip
