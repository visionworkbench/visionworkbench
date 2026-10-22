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

}} // end namespace vw::ip
