// __BEGIN_LICENSE__
//  Copyright (c) 2009-2026, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NGT platform is licensed under the Apache License, Version 2.0 (the
//  "License"); you may not use this file except in compliance with the
//  License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__

/// \file SGMCosts.cc
/// Disparity cost computation routines for Semi-Global Matching

#include <vw/Stereo/SGM.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/CensusTransform.h>
#include <vw/Image/PixelMask.h>

#include <cmath>

namespace vw {

// Helper function
template <typename T>
void get_hamming_distance_costs(ImageView<T> const& left_binary_image,
                                ImageView<T> const& right_binary_image,
                                int min_row, int max_row,
                                int min_col, int max_col,
                                int kernel_size,
                                ImageView<Vector4i> const& disp_bound_image,
                                boost::shared_array<SemiGlobalMatcher::CostType> cost_buf) {

  const int half_kernel = (kernel_size - 1) / 2;

  // Now compute the disparity costs for each pixel.
  // Make sure we don't go out of bounds here due to the disparity shift and kernel.
  size_t cost_index = 0;
  for (int r = min_row; r <= max_row; r++) { // For each row in left
    int output_row = r - min_row;
    int binary_row = r - half_kernel;
    for (int c = min_col; c <= max_col; c++) { // For each column in left
      int output_col = c - min_col;
      int binary_col = c - half_kernel;

      Vector4i pixel_disp_bounds = disp_bound_image(output_col, output_row);

      for (int dy = pixel_disp_bounds[1]; dy <= pixel_disp_bounds[3]; dy++) { // For each disparity
        for (int dx = pixel_disp_bounds[0]; dx <= pixel_disp_bounds[2]; dx++) {

          SemiGlobalMatcher::CostType cost
            = hamming_distance(left_binary_image(binary_col, binary_row),
                               right_binary_image(binary_col+dx, binary_row+dy));
          cost_buf[cost_index] = cost;
          ++cost_index;
        }
      } // End disparity loops
    } // End x loop
  } // End y loop
}

void SemiGlobalMatcher::compute_patch_mean_std(ImageView<uint8> const& image, int x, int y,
                            double &mean, double &std) const {

  if (m_kernel_size == 1) {
    mean = static_cast<double>(image(x, y));
    std  = 1;
    return;
  }

  const int half_kernel_size = (m_kernel_size-1) / 2;

  // Compute the mean
  mean = 0;
  double count=0;
  for (int j=-half_kernel_size; j<=half_kernel_size; ++j) {
    for (int i=-half_kernel_size; i<=half_kernel_size; i++) {
      mean += static_cast<double>(image(x +i, y +j));
      count += 1.0;
    }
  }
  mean /= count;

  // Compute the STD
  std = 0;
  for (int j=-half_kernel_size; j<=half_kernel_size; ++j) {
    for (int i=-half_kernel_size; i<=half_kernel_size; i++) {
      double val = static_cast<double>(image(x +i, y +j)) - mean;
      std += val*val;
    }
  }
  std = sqrt(std / (count - 1.0));
} // End function compute_patch_mean

// TODO: Replace with ASP implementation?
SemiGlobalMatcher::CostType
SemiGlobalMatcher::get_cost_block(ImageView<uint8> const& left_image,
                                  ImageView<uint8> const& right_image,
                                  double mean1, double std1,
                                  int left_x, int left_y, int right_x, int right_y,
                                  bool debug) const {

  if (m_kernel_size == 1) { // Special handling for single pixel case
    int diff = static_cast<int>(left_image (left_x,  left_y)) -
               static_cast<int>(right_image(right_x, right_y));
    return static_cast<CostType>(std::abs(diff));
  }

  const CostType MAX_RESULT = 255; // For uint8
  const int      half_kernel_size = (m_kernel_size-1) / 2;
  const int      count = static_cast<int>(m_kernel_size*m_kernel_size);

  // Normalized cross correlation
  if (m_cost_type == 2) {

    // Already have stats for the left image patch, compute them for the right image patch.
    double mean2, std2;
    compute_patch_mean_std(right_image, right_x, right_y, mean2, std2);

    double sum=0;
    double coeff = 1.0 / (std1*std2);
    for (int j=-half_kernel_size; j<=half_kernel_size; ++j) {
      for (int i=-half_kernel_size; i<=half_kernel_size; i++) {
        double a = static_cast<double>(left_image (left_x +i, left_y +j)) - mean1;
        double b = static_cast<double>(right_image(right_x+i, right_y+j)) - mean2;
        sum += a*b*coeff;
      }
    }
    sum /= static_cast<double>(count); // Result is in the range: -1 to +1

    // Saturate the result
    if (sum > MAX_RESULT)
      return MAX_RESULT;
    const double OFFSET = 128;
    return static_cast<CostType>(sum*OFFSET + OFFSET); // Convert to 0-255 range.
  } // End normalized cross correlation case

  // Mean of abs differences
  int sum=0, diff=0;
  for (int j=-half_kernel_size; j<=half_kernel_size; ++j) {
    for (int i=-half_kernel_size; i<=half_kernel_size; i++) {
      diff = static_cast<int>(left_image (left_x +i, left_y +j)) -
             static_cast<int>(right_image(right_x+i, right_y+j));
      sum += std::abs(diff);
    }
  }
  sum /= count;
  if (sum > MAX_RESULT)
    return MAX_RESULT;

  return static_cast<CostType>(sum);
}

void SemiGlobalMatcher::fill_costs_block(ImageView<uint8> const& left_image,
                                         ImageView<uint8> const& right_image) {
  // Make sure we don't go out of bounds here due to the disparity shift and kernel.
  size_t cost_index = 0;
  for (int r = m_min_row; r <= m_max_row; r++) { // For each row in left
    int output_row = r - m_min_row;
    for (int c = m_min_col; c <= m_max_col; c++) { // For each column in left
      int output_col = c - m_min_col;

      Vector4i pixel_disp_bounds = m_disp_bound_image(output_col, output_row);

      // These statistics are used repeatedly for this particular cost type.
      double mean_left=0, std_left=1;
      if (m_cost_type == 2)
        compute_patch_mean_std(left_image, c, r, mean_left, std_left);

      // Only compute costs in the search radius for this pixel
      for (int dy = pixel_disp_bounds[1]; dy <= pixel_disp_bounds[3]; dy++) { // For each disparity
        for (int dx = pixel_disp_bounds[0]; dx <= pixel_disp_bounds[2]; dx++) {

          CostType cost = get_cost_block(left_image, right_image, mean_left, std_left,
                                         c, r, c+dx, r+dy, false);
          m_cost_buffer[cost_index] = cost;
          ++cost_index;
        }
      } // End disparity loops
    } // End x loop
  } // End y loop
}

void SemiGlobalMatcher::fill_costs_census3x3(ImageView<uint8> const& left_image,
                                             ImageView<uint8> const& right_image) {
  const int half_kernel = (m_kernel_size - 1) / 2;
  const int padding     = 2*half_kernel;

  if (m_cost_type == CENSUS_TRANSFORM) {
    ImageView<uint8> left_census (left_image.cols()-padding,  left_image.rows()-padding),
                     right_census(right_image.cols()-padding, right_image.rows()-padding);

    for (int r = 0; r < left_census.rows(); r++)
      for (int c = 0; c < left_census.cols(); c++)
        left_census(c,r) = get_census_value_3x3(left_image, c+half_kernel, r+half_kernel);
    for (int r = 0; r < right_census.rows(); r++)
      for (int c = 0; c < right_census.cols(); c++)
        right_census(c,r) = get_census_value_3x3(right_image, c+half_kernel, r+half_kernel);
    get_hamming_distance_costs(left_census, right_census, m_min_row, m_max_row, m_min_col,
                               m_max_col, m_kernel_size, m_disp_bound_image, m_cost_buffer);
  } else {
    ImageView<uint16> left_census (left_image.cols()-padding,  left_image.rows()-padding),
                      right_census(right_image.cols()-padding, right_image.rows()-padding);

    for (int r = 0; r < left_census.rows(); r++)
      for (int c = 0; c < left_census.cols(); c++)
        left_census(c,r) = get_census_value_ternary_3x3(left_image, c+half_kernel, r+half_kernel, m_ternary_census_threshold);
    for (int r = 0; r < right_census.rows(); r++)
      for (int c = 0; c < right_census.cols(); c++)
        right_census(c,r) = get_census_value_ternary_3x3(right_image, c+half_kernel, r+half_kernel, m_ternary_census_threshold);
    get_hamming_distance_costs(left_census, right_census, m_min_row, m_max_row, m_min_col,
                               m_max_col, m_kernel_size, m_disp_bound_image, m_cost_buffer);
  }
}

void SemiGlobalMatcher::fill_costs_census5x5(ImageView<uint8> const& left_image,
                                             ImageView<uint8> const& right_image) {
  const int half_kernel = (m_kernel_size - 1) / 2;
  const int padding     = 2*half_kernel;

  ImageView<uint32> left_census (left_image.cols()-padding,  left_image.rows()-padding),
                    right_census(right_image.cols()-padding, right_image.rows()-padding);

  if (m_cost_type == CENSUS_TRANSFORM) {
    for (int r = 0; r < left_census.rows(); r++)
      for (int c = 0; c < left_census.cols(); c++)
        left_census(c,r) = get_census_value_5x5(left_image, c+half_kernel, r+half_kernel);
    for (int r = 0; r < right_census.rows(); r++)
      for (int c = 0; c < right_census.cols(); c++)
        right_census(c,r) = get_census_value_5x5(right_image, c+half_kernel, r+half_kernel);
  } else { // TERNARY_CENSUS_TRANSFORM
    for (int r = 0; r < left_census.rows(); r++)
      for (int c = 0; c < left_census.cols(); c++)
        left_census(c,r) = get_census_value_ternary_5x5(left_image, c+half_kernel, r+half_kernel, m_ternary_census_threshold);
    for (int r = 0; r < right_census.rows(); r++)
      for (int c = 0; c < right_census.cols(); c++)
        right_census(c,r) = get_census_value_ternary_5x5(right_image, c+half_kernel, r+half_kernel, m_ternary_census_threshold);
  }
  get_hamming_distance_costs(left_census, right_census, m_min_row, m_max_row, m_min_col,
                             m_max_col, m_kernel_size, m_disp_bound_image, m_cost_buffer);
}

void SemiGlobalMatcher::fill_costs_census7x7(ImageView<uint8> const& left_image,
                                             ImageView<uint8> const& right_image) {
  const int half_kernel = (m_kernel_size - 1) / 2;
  const int padding     = 2*half_kernel;

  ImageView<uint64> left_census (left_image.cols()-padding,  left_image.rows()-padding),
                    right_census(right_image.cols()-padding, right_image.rows()-padding);

  if (m_cost_type == CENSUS_TRANSFORM) {
    for (int r = 0; r < left_census.rows(); r++)
      for (int c = 0; c < left_census.cols(); c++)
        left_census(c,r) = get_census_value_7x7(left_image, c+half_kernel, r+half_kernel);
    for (int r = 0; r < right_census.rows(); r++)
      for (int c = 0; c < right_census.cols(); c++)
        right_census(c,r) = get_census_value_7x7(right_image, c+half_kernel, r+half_kernel);
  } else { // TERNARY_CENSUS_TRANSFORM
    for (int r = 0; r < left_census.rows(); r++)
      for (int c = 0; c < left_census.cols(); c++)
        left_census(c,r) = get_census_value_ternary_7x7(left_image, c+half_kernel, r+half_kernel, m_ternary_census_threshold);
    for (int r = 0; r < right_census.rows(); r++)
      for (int c = 0; c < right_census.cols(); c++)
        right_census(c,r) = get_census_value_ternary_7x7(right_image, c+half_kernel, r+half_kernel, m_ternary_census_threshold);
  }
  get_hamming_distance_costs(left_census, right_census, m_min_row, m_max_row, m_min_col,
                             m_max_col, m_kernel_size, m_disp_bound_image, m_cost_buffer);
}

void SemiGlobalMatcher::fill_costs_census9x9(ImageView<uint8> const& left_image,
                                             ImageView<uint8> const& right_image) {
  const int half_kernel = (m_kernel_size - 1) / 2;
  const int padding     = 2*half_kernel;

  if (m_cost_type == CENSUS_TRANSFORM) {
    ImageView<uint32> left_census (left_image.cols()-padding,  left_image.rows()-padding),
                      right_census(right_image.cols()-padding, right_image.rows()-padding);

    for (int r = 0; r < left_census.rows(); r++)
      for (int c = 0; c < left_census.cols(); c++)
        left_census(c,r) = get_census_value_9x9(left_image, c+half_kernel, r+half_kernel);
    for (int r = 0; r < right_census.rows(); r++)
      for (int c = 0; c < right_census.cols(); c++)
        right_census(c,r) = get_census_value_9x9(right_image, c+half_kernel, r+half_kernel);
    get_hamming_distance_costs(left_census, right_census, m_min_row, m_max_row, m_min_col,
                               m_max_col, m_kernel_size, m_disp_bound_image, m_cost_buffer);
  } else { // TERNARY_CENSUS_TRANSFORM
    ImageView<uint64> left_census (left_image.cols()-padding,  left_image.rows()-padding),
                      right_census(right_image.cols()-padding, right_image.rows()-padding);

    for (int r = 0; r < left_census.rows(); r++)
      for (int c = 0; c < left_census.cols(); c++)
        left_census(c,r) = get_census_value_ternary_9x9(left_image, c+half_kernel, r+half_kernel, m_ternary_census_threshold);
    for (int r = 0; r < right_census.rows(); r++)
      for (int c = 0; c < right_census.cols(); c++)
        right_census(c,r) = get_census_value_ternary_9x9(right_image, c+half_kernel, r+half_kernel, m_ternary_census_threshold);
    get_hamming_distance_costs(left_census, right_census, m_min_row, m_max_row, m_min_col,
                               m_max_col, m_kernel_size, m_disp_bound_image, m_cost_buffer);
  }
}

// TODO: Add multithreading capability to this function!
void SemiGlobalMatcher::compute_disparity_costs(ImageView<uint8> const& left_image,
                                                ImageView<uint8> const& right_image) {

  // Call the appropriate cost function
  switch(m_cost_type) {
  case CENSUS_TRANSFORM:
  case TERNARY_CENSUS_TRANSFORM:
    switch(m_kernel_size) {
    case 3:  fill_costs_census3x3(left_image, right_image); break;
    case 5:  fill_costs_census5x5(left_image, right_image); break;
    case 7:  fill_costs_census7x7(left_image, right_image); break;
    case 9:  fill_costs_census9x9(left_image, right_image); break;
    default:
      vw_throw(ArgumentErr() << "Unsupported kernel size for census transform: "
                             << m_kernel_size);
    };
    break;
  default: // All other cost functions
    fill_costs_block(left_image, right_image);
  };

} // end compute_disparity_costs()

} // end namespace vw
