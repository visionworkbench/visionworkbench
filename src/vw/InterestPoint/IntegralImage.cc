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

#include <vw/InterestPoint/IntegralImage.h>
#include <cmath>

namespace vw {
namespace ip {

/// Using an integral image, compute the summed value of a region
/// in the original image.
float
IntegralBlock(ImageView<float> const& integral,
              Vector2i         const& top_left,
              Vector2i         const& bottom_right) {
  VW_DEBUG_ASSERT(top_left.x() < integral.cols(),
                  vw::ArgumentErr() << "x0 out of bounds. "
                  << integral.cols() << " : "
                  << top_left << bottom_right << "\n");
  VW_DEBUG_ASSERT(bottom_right.x() < integral.cols(),
                  vw::ArgumentErr() << "x1 out of bounds. "
                  << integral.cols() << " : "
                  << top_left << bottom_right << "\n");
  VW_DEBUG_ASSERT(top_left.y() < integral.rows(),
                  vw::ArgumentErr() << "y0 out of bounds. "
                  << integral.rows() << " : "
                  << top_left << bottom_right << "\n");
  VW_DEBUG_ASSERT(bottom_right.y() < integral.rows(),
                  vw::ArgumentErr() << "y1 out of bounds. "
                  << integral.rows() << " : "
                  << top_left << bottom_right << "\n");

  float result = 0;
  result = integral(top_left.x(), top_left.y());
  result += integral(bottom_right.x(), bottom_right.y());
  result -= integral(top_left.x(), bottom_right.y());
  result -= integral(bottom_right.x(), top_left.y());

  return result;
}

/// X First Derivative (not implemented)
float
XFirstDerivative(ImageView<float> const& /*integral*/,
                 int const& /*x*/, int const& /*y*/,
                 unsigned const& /*filter_size*/) {
  vw_throw(vw::NoImplErr()
           << "First derivative filter has not been implemented yet\n");
  float derivative = 0;
  return derivative;
}

/// X Second Derivative
float
XSecondDerivative(ImageView<float> const& integral,
                  int const& x, int const& y,
                  unsigned const& filter_size) {
  unsigned lobe      = filter_size / 3;
  unsigned half_lobe = (unsigned)floor(float(lobe) / 2.0);
  float derivative = 0;

  // Adding positive left
  derivative = IntegralBlock(integral,
                             Vector2i(x - lobe - half_lobe,
                                      y - lobe + 1),
                             Vector2i(x - half_lobe, y + lobe));

  // Adding negative middle
  derivative -= 2.0f * IntegralBlock(integral,
                                     Vector2i(x - half_lobe,
                                              y - lobe + 1),
                                     Vector2i(x + half_lobe + 1,
                                              y + lobe));

  // Adding positive right
  derivative += IntegralBlock(integral,
                              Vector2i(x + half_lobe + 1,
                                       y - lobe + 1),
                              Vector2i(x + half_lobe + lobe + 1,
                                       y + lobe));

  derivative /= filter_size * filter_size;

  return derivative;
}

/// Y First Derivative (not implemented)
float
YFirstDerivative(ImageView<float> const& /*integral*/,
                 int const& /*x*/, int const& /*y*/,
                 unsigned const& /*filter_size*/) {
  vw_throw(vw::NoImplErr()
           << "First derivative filter has not been implemented yet\n");
  float derivative = 0;
  return derivative;
}

/// Y Second Derivative
float
YSecondDerivative(ImageView<float> const& integral,
                  int const& x, int const& y,
                  unsigned const& filter_size) {
  unsigned lobe = filter_size / 3;
  unsigned half_lobe = (unsigned)floor(float(lobe) / 2.0);
  float derivative = 0;

  // Adding positive top
  derivative = IntegralBlock(integral,
                             Vector2i(x - lobe + 1,
                                      y - lobe - half_lobe),
                             Vector2i(x + lobe, y - half_lobe));

  // Adding negative middle
  derivative -= 2.0f * IntegralBlock(integral,
                                     Vector2i(x - lobe + 1,
                                              y - half_lobe),
                                     Vector2i(x + lobe,
                                              y + half_lobe + 1));

  // Adding positive bottom
  derivative += IntegralBlock(integral,
                              Vector2i(x - lobe + 1,
                                       y + half_lobe + 1),
                              Vector2i(x + lobe,
                                       y + half_lobe + lobe + 1));

  derivative /= filter_size * filter_size;

  return derivative;
}

/// XY Derivative
float
XYDerivative(ImageView<float> const& integral,
             int const& x, int const& y,
             unsigned const& filter_size) {

  unsigned lobe = filter_size / 3;
  float derivative = 0;

  // Adding positive top left
  derivative = IntegralBlock(integral,
                             Vector2i(x - lobe, y - lobe),
                             Vector2i(x, y));

  // Adding negative top right
  derivative -= IntegralBlock(integral,
                              Vector2i(x + 1, y - lobe),
                              Vector2i(x + lobe + 1, y));

  // Adding negative bottom left
  derivative -= IntegralBlock(integral,
                              Vector2i(x - lobe, y + 1),
                              Vector2i(x, y + lobe + 1));

  // Adding positive bottom right
  derivative += IntegralBlock(integral,
                              Vector2i(x + 1, y + 1),
                              Vector2i(x + 1 + lobe,
                                       y + 1 + lobe));

  derivative /= filter_size * filter_size;

  return derivative;
}

}} // end namespace vw::ip
