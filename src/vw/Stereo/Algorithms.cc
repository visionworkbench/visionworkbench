// __BEGIN_LICENSE__
//  Copyright (c) 2006-2025, United States Government as represented by the
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

#include <vw/Stereo/Algorithms.h>
#include <vw/Math/Functors.h>

#include <vector>

namespace vw {
namespace stereo {

void disparity_median_filter(ImageView<PixelMask<Vector2f>> const& disparity_in,
                             ImageView<PixelMask<Vector2f>>      & disparity_out,
                             int kernel_size) {

  int half_kernel = (kernel_size - 1) / 2;
  int num_vals = kernel_size * kernel_size;
  disparity_out = disparity_in;

  if (kernel_size < 3) // No smoothing called for
    return;

  // Output pixel loop
  for (int row = half_kernel; row < disparity_in.rows() - half_kernel; ++row) {
    for (int col = half_kernel; col < disparity_in.cols() - half_kernel; ++col) {

      if (!is_valid(disparity_in(col, row)))
        continue;

      // Loop through the kernel
      std::vector<double> dx(num_vals), dy(num_vals);
      int index = 0;
      for (int r = row - half_kernel; r <= row + half_kernel; ++r) {
        for (int c = col - half_kernel; c <= col + half_kernel; ++c) {
          if (is_valid(disparity_in(c, r))) {
            dx[index] = disparity_in(c, r)[0];
            dy[index] = disparity_in(c, r)[1];
            ++index;
          }
        }
      }
      if (index == 0)
        continue;

      dx.resize(index);
      dy.resize(index);
      double median_x = math::destructive_median(dx);
      double median_y = math::destructive_median(dy);

      disparity_out(col, row) = PixelMask<Vector2f>(median_x, median_y);
    }
  } // End loop through pixels
}

void disparity_neighbor_filter(ImageView<PixelMask<Vector2i>> const& disparity_in,
                               ImageView<PixelMask<Vector2i>>      & disparity_out) {

  const int COPY_COUNT = 5;

  disparity_out = disparity_in;

  std::vector<int> counts(8);
  std::vector<PixelMask<Vector2i>> vals(8);
  for (int row = 1; row < disparity_in.rows() - 1; ++row) {
    for (int col = 1; col < disparity_in.cols() - 1; ++col) {

      vals[0] = disparity_in(col-1, row-1);
      vals[1] = disparity_in(col,   row-1);
      vals[2] = disparity_in(col+1, row-1);
      vals[3] = disparity_in(col-1, row  );
      vals[4] = disparity_in(col+1, row  );
      vals[5] = disparity_in(col-1, row+1);
      vals[6] = disparity_in(col,   row+1);
      vals[7] = disparity_in(col+1, row+1);

      int max_count = 0;
      int max_index = 0;
      for (int i = 0; i < 8; ++i) {
        counts[i] = 0;
        if (!is_valid(vals[i]))
          continue;
        for (int j = 0; j < 8; ++j) {
          if (vals[i] == vals[j])
            ++counts[i];
        }
        if (counts[i] > max_count) {
          max_count = counts[i];
          max_index = i;
        }
      } // End max finder

      if (max_count >= COPY_COUNT)
        disparity_out(col, row) = vals[max_index];
    }
  } // End loop through pixels
}

}} // namespace vw::stereo
