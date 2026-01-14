// __BEGIN_LICENSE__
//  Copyright (c) 2009-2025, United States Government as represented by the
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

#include <queue>
#include <math.h>
#include <vw/Stereo/SGM.h>
#include <vw/Stereo/SGMAssist.h>
#include <vw/Core/Debugging.h>
#include <vw/Core/ThreadPool.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/Algorithms.h>

#if defined(VW_ENABLE_SSE) && (VW_ENABLE_SSE==1)
  #include <emmintrin.h>
  #include <smmintrin.h> // SSE4.1
#endif

namespace vw {

namespace stereo {

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

void SemiGlobalMatcher::set_parameters(CostFunctionType cost_type,
                                       bool use_mgm,
                                       int min_disp_x, int min_disp_y,
                                       int max_disp_x, int max_disp_y,
                                       int kernel_size,
                                       SgmSubpixelMode subpixel,
                                       Vector2i search_buffer,
                                       size_t memory_limit_mb,
                                       uint16 p1, uint16 p2,
                                       int ternary_census_threshold) {
  m_cost_type   = cost_type;
  m_use_mgm     = use_mgm;
  m_min_disp_x  = min_disp_x;
  m_min_disp_y  = min_disp_y;
  m_max_disp_x  = max_disp_x;
  m_max_disp_y  = max_disp_y;
  m_kernel_size = kernel_size;
  m_ternary_census_threshold = ternary_census_threshold;
  m_subpixel_type   = subpixel;
  m_search_buffer   = search_buffer;
  m_memory_limit_mb = memory_limit_mb;

  m_num_disp_x = m_max_disp_x - m_min_disp_x + 1;
  m_num_disp_y = m_max_disp_y - m_min_disp_y + 1;
  size_t size_check = m_num_disp_x * m_num_disp_y;;
  if (size_check > (size_t)std::numeric_limits<DisparityType>::max())
    vw_throw(NoImplErr() << "Number of disparities is too large for data type.\n");
  m_num_disp = m_num_disp_x * m_num_disp_y;

  if (p1 > 0) // User provided
    m_p1 = p1;
  else { // Choose based on the kernel size and cost type
    switch(cost_type) {
    case CENSUS_TRANSFORM:
      switch(kernel_size) {
      case 3:  m_p1 = 3;  break; // 0 - 8
      case 5:  m_p1 = 15; break; // 0 - 24
      case 7:  m_p1 = 30; break; // 0 - 48
      case 9:  m_p1 = 20; break; // 0 - 32
      default: m_p1 = 3; break; // Unsupported!
      };
      break;
    case TERNARY_CENSUS_TRANSFORM:
      switch(kernel_size) {
      case 3:  m_p1 = 12; break; // 0 - 16
      case 5:  m_p1 = 30; break; // 0 - 48
      case 7:  m_p1 = 40; break; // 0 - 64
      case 9:  m_p1 = 40; break; // 0 - 64
      default: m_p1 = 30; break; // Unsupported!
      };
      break;
    default: // MAD
      m_p1 = 3; // 0 - 255 range
      break;
    };
  } // End p1 cases

  if (p2 > 0) // User provided
    m_p2 = p2;
  else {
    switch(cost_type) {
    case CENSUS_TRANSFORM:
      switch(kernel_size) {
      case 3:  m_p2 = 70;   break; // 0 - 8
      case 5:  m_p2 = 750;  break; // 0 - 24
      case 7:  m_p2 = 1500; break; // 0 - 48
      case 9:  m_p2 = 1000; break; // 0 - 32
      default: m_p2 = 22;   break; // Unsupported!
      };
      break;
    case TERNARY_CENSUS_TRANSFORM:
      switch(kernel_size) {
      case 3:  m_p2 = 600;  break; // 0 - 16
      case 5:  m_p2 = 1500; break; // 0 - 48
      case 7:  m_p2 = 2000; break; // 0 - 64
      case 9:  m_p2 = 2000; break; // 0 - 64
      default: m_p2 = 30;   break; // Unsupported!
      };
      break;
    default: // MAD
      m_p2 = 250; // 0 - 255 range
      break;
    };
  } // End p2 cases

  vw_out(DebugMessage, "stereo") << "SGM m_p1 = " << m_p1 << "\n";
  vw_out(DebugMessage, "stereo") << "SGM m_p2 = " << m_p2 << "\n";

} // End set_parameters

ImageView<PixelMask<Vector2i>> calc_disparity_sgm(
  CostFunctionType cost_type,
  ImageView<PixelGray<float>> const& left_in,
  ImageView<PixelGray<float>> const& right_in,
  BBox2i                 const& left_region,   // Valid region in the left image
  Vector2i               const& search_volume, // Max disparity to search in right image
  Vector2i               const& kernel_size,   // The kernel dimensions are always equal
  bool                   const  use_mgm,
  SemiGlobalMatcher::SgmSubpixelMode const& subpixel_mode,
  Vector2i               const  search_buffer, // Search buffer applied around prev_disparity
  size_t                 const  memory_limit_mb,
  boost::shared_ptr<SemiGlobalMatcher> &matcher_ptr,
  ImageView<uint8>       const* left_mask_ptr,
  ImageView<uint8>       const* right_mask_ptr,
  SemiGlobalMatcher::DisparityImage const* prev_disparity) {

  // Sanity checks
  VW_DEBUG_ASSERT(kernel_size[0] % 2 == 1 && kernel_size[1] % 2 == 1,
                    ArgumentErr() << "calc_disparity_sgm: Kernel input not sized with odd values.");
  VW_DEBUG_ASSERT(kernel_size[0] <= left_region.width() &&
                    kernel_size[1] <= left_region.height(),
                    ArgumentErr() << "calc_disparity_sgm: Kernel size too large of active region.");
  VW_DEBUG_ASSERT(left_region.min().x() >= 0 &&  left_region.min().y() >= 0 &&
                    left_region.max().x() <= left_in.impl().cols() &&
                    left_region.max().y() <= left_in.impl().rows(),
                    ArgumentErr() << "calc_disparity_sgm: Region not inside left image.");

  Vector2i search_volume_inclusive = search_volume;

  // Rasterize input so that we can do a lot of processing on it.
  BBox2i right_region = left_region;
  right_region.max() += search_volume_inclusive;

  vw_out(VerboseDebugMessage, "stereo")
    << "calc_disparity_sgm: left  region = " << left_region << "\n";
  vw_out(VerboseDebugMessage, "stereo")
    << "calc_disparity_sgm: right region = " << right_region << "\n";
  vw_out(VerboseDebugMessage, "stereo")
    << "calc_disparity_sgm: search_volume_inclusive = " << search_volume_inclusive << "\n";

  // TODO: Ignore masked values when computing this!
  // Convert the input image to uint8
  // - Any "smart" stretching here can cause problems.
  ImageView<PixelGray<vw::uint8>> left, right;
  u8_convert(crop(left_in.impl(), left_region),  left);
  u8_convert(crop(right_in.impl(), right_region), right);

  try {
    matcher_ptr.reset(new SemiGlobalMatcher(cost_type, use_mgm, 0, 0,
                      search_volume_inclusive[0], search_volume_inclusive[1],
                      kernel_size[0], subpixel_mode, search_buffer, memory_limit_mb));
    return matcher_ptr->semi_global_matching_func(left, right, left_mask_ptr,
                                                  right_mask_ptr, prev_disparity);

  } catch (const std::exception& e) {
    vw::vw_throw(vw::ArgumentErr()
      << "Failed to compute the correlation. See the online documentation "
      << "(next_steps.html) for how to handle failures.\n"
      << "Detailed error message: " << e.what() << "\n");
  }

  return ImageView<PixelMask<Vector2i>>();
} // End function calc_disparity

void SemiGlobalMatcher::populate_constant_disp_bound_image() {
  // Allocate the image
  m_disp_bound_image.set_size(m_num_output_cols, m_num_output_rows);
  // Fill it up with an identical vector
  Vector4i bounds_vector(m_min_disp_x, m_min_disp_y, m_max_disp_x, m_max_disp_y);
  size_t buffer_size = m_num_output_cols*m_num_output_rows;

  std::fill(m_disp_bound_image.data(), m_disp_bound_image.data()+buffer_size, bounds_vector);
}

bool SemiGlobalMatcher::populate_disp_bound_image(ImageView<uint8> const* left_image_mask,
                                                  ImageView<uint8> const* right_image_mask,
                                                  DisparityImage const* prev_disparity) {

  // The masks are assumed to be the same size as the output image.
  // TODO: Check or automatically compute the left valid mask size!
  bool left_mask_valid = false, right_mask_valid = false;
  if (left_image_mask) {
    vw_out(VerboseDebugMessage, "stereo") << "Left  mask image size:" << bounding_box(*left_image_mask) << "\n";
    // Currently the left mask is required to EXACTLY match the output size...
    if ((left_image_mask->cols() == m_disp_bound_image.cols()) &&
         (left_image_mask->rows() == m_disp_bound_image.rows()))
      left_mask_valid = true;
    else
      vw_throw(LogicErr() << "Left mask size does not match the output size.\n");
  }
  if (right_image_mask) {
    vw_out(VerboseDebugMessage, "stereo") << "Right mask image size:" << bounding_box(*right_image_mask) << "\n";

    // Currently this class assumes a positive search range and requires the right mask to
    //  be big enough to contain the output size PLUS the search range.
    if ((right_image_mask->cols() >= m_disp_bound_image.cols()+m_num_disp_x-1) &&
         (right_image_mask->rows() >= m_disp_bound_image.rows()+m_num_disp_y-1))
      right_mask_valid = true;
    else
      vw_throw(LogicErr()
               << "Right mask size is not large enough to support search range.\n");
  }

  // The low-res disparity image must be half-resolution.
  const int SCALE_UP = 2;

  // Require that the right image mask is valid for this percentage of the
  // search range for each pixel in the left image.
  // - If this value is too low, many border pixels will be assigned the full
  //   search range and significantly slow down SGM!
  //const double MIN_MASK_OVERLAP = 0.95;

  // There needs to be some "room" in the disparity search space for
  // us to discard prior results on the edge as not trustworthy predictors.
  // In other words, don't mark any pixels as edge if the search space is too small!
  // - The edges come from the previous resolution level, so don't check edges unless
  //   the previous level was at least 5 pixels in that axis.
  const int disp_range_width  = (m_max_disp_x - m_min_disp_x + 1);
  const int disp_range_height = (m_max_disp_y - m_min_disp_y + 1);
  const bool check_x_edge = (disp_range_width  >= 10);
  const bool check_y_edge = (disp_range_height >= 10);

  double area = 0, percent_trusted = 0, percent_masked = 0;

  // NOTE: We don't use this class anymore because the right mask checking it was doing
  //       was too conservative in checking for right image area and was causing valid
  //       stereo area to go unprocessed.
  // This class will check the right image mask in an efficient manner.
  // - This implementation relies on >= 0 disparity ranges!
  //IterativeMaskBoxCounter right_mask_checker(right_image_mask, Vector2i(m_num_disp_x, m_num_disp_y));

  const Vector4i ZERO_SEARCH_AREA(0, 0, -1, -1);

  ImageView<uint8> full_search_image(m_disp_bound_image.cols(), m_disp_bound_image.rows());

  // Find the upper and lower bounds of valid pixels in the right mask image.
  int min_valid_right_row = 0, max_valid_right_row = 0;

  if (right_mask_valid) {

    min_valid_right_row = right_image_mask->rows()-1;

    int num_cols = m_disp_bound_image.cols();

    // For each column, find the maximum and minimum valid rows.
    for (int c=0; c<num_cols; ++c) {

      for (int i=right_image_mask->rows()-1; i>0; --i) { // Top to bottom
        if (right_image_mask->operator()(c, i) > 0) {
          if (i > max_valid_right_row)
            max_valid_right_row = i;
          break; // Stop at first pixel
        }
      } // End first row loop
      for (int i = 0; i < right_image_mask->rows(); i++) { // Bottom to top
        if (right_image_mask->operator()(c, i) > 0) {
          if (i < min_valid_right_row)
            min_valid_right_row = i;
          break; // Stop at first pixel
        }
      } // End second row loop
    } // End col loop

  } // End mask valid case

  // Loop through the output disparity image and compute a search range for each pixel
  int r_in, c_in;
  int dx_scaled, dy_scaled;
  PixelMask<Vector2i> input_disp;
  Vector4i bounds;
  double missing_prev_disp_count = 0;
  for (int r=0; r<m_disp_bound_image.rows(); ++r) { // Start loop through rows
    r_in = r / SCALE_UP;

    int min_valid_right_column = -1, max_valid_right_column = -2;
    if (right_mask_valid) {
      // For this row, find the maximum and minimum valid columns.
      for (int i=right_image_mask->cols()-1; i>0; --i) {
        if (right_image_mask->operator()(i, r) > 0) {
          max_valid_right_column = i;
          break;
        }
      }
      if (max_valid_right_column > 0) { // Only do if there are valid pixels
        for (int i = 0; i < right_image_mask->cols(); i++) {
          if (right_image_mask->operator()(i, r) > 0) {
            min_valid_right_column = i;
            break;
          }
        }
      }
    } // End mask column checking

    for (int c=0; c<m_disp_bound_image.cols(); ++c) {
/*
      // TODO: This will fail if not used with positive search ranges!
      // Verify that there is sufficient overlap with the right image mask
      if (right_mask_valid) {
        double right_percent  = right_mask_checker.next_pixel();
        bool   right_check_ok = (right_percent >= MIN_MASK_OVERLAP);
        // If none of the right mask pixels were valid, flag this pixel as invalid.
        if (!right_check_ok) {
          m_disp_bound_image(c,r) = ZERO_SEARCH_AREA;
          ++percent_masked;
          continue;
        }
      } // End right mask handling
*/
      // If the left mask is invalid here, flag the pixel as invalid.
      //// - Do this second so that our right image pixel tracker stays up to date.
      if (left_mask_valid && (left_image_mask->operator()(c,r) == 0)) {
        m_disp_bound_image(c,r) = ZERO_SEARCH_AREA;
        ++percent_masked;
        continue;
      }

      // If a previous disparity was provided, see if we have a valid disparity
      // estimate for this pixel.
      bool good_disparity = false;
      c_in = c / SCALE_UP; // Pixel location in prior disparity map if it exists

      if (prev_disparity) {

        // Verify that the pixel we want exists
        if ((c_in >= prev_disparity->cols()) || (r_in >= prev_disparity->rows())) {
          missing_prev_disp_count += 1;
        } else {
          input_disp = prev_disparity->operator()(c_in,r_in);
          // Disparity values on the edge of our 2D search range are not considered trustworthy!
          dx_scaled = input_disp[0] * SCALE_UP;
          dy_scaled = input_disp[1] * SCALE_UP;

          bool on_edge = ((check_x_edge && ((dx_scaled <= m_min_disp_x)   ||
                                            (dx_scaled >= m_max_disp_x))) ||
                          (check_y_edge && ((dy_scaled <= m_min_disp_y)   ||
                                            (dy_scaled >= m_max_disp_y))));
          good_disparity = (is_valid(input_disp) && !on_edge);
        }

      } // End prev disparity check

      if (good_disparity) {

        // We are more confident in the prior disparity, search nearby.
        bounds[0]  = dx_scaled - m_search_buffer[0]; // Min x
        bounds[2]  = dx_scaled + m_search_buffer[0]; // Max X
        bounds[1]  = dy_scaled - m_search_buffer[1]; // Min y
        bounds[3]  = dy_scaled + m_search_buffer[1]; // Max y

        // Constrain to global limits
        if (bounds[0] < m_min_disp_x) bounds[0] = m_min_disp_x;
        if (bounds[1] < m_min_disp_y) bounds[1] = m_min_disp_y;
        if (bounds[2] > m_max_disp_x) bounds[2] = m_max_disp_x;
        if (bounds[3] > m_max_disp_y) bounds[3] = m_max_disp_y;

        percent_trusted += 1.0;
      } else {
        // Not a trusted prior disparity, search the entire range!
        // - This takes a long time.
        bounds = Vector4i(m_min_disp_x, m_min_disp_y, m_max_disp_x, m_max_disp_y);
        full_search_image(c,r) = 255;
      }

      // Restrict search range to the right image mask
      // - This could be improved and more efficient!
      if (right_mask_valid) {

        // Get the intersection of the valid pixels and the offsets
        BBox2i valid_region(Vector2i(min_valid_right_column, min_valid_right_row),
                            Vector2i(max_valid_right_column, max_valid_right_row));
        BBox2i valid_offsets = valid_region - Vector2i(c,r);
        BBox2i full_offsets(Vector2i(bounds[0], bounds[1]),
                            Vector2i(bounds[2], bounds[3]));
        valid_offsets.crop(full_offsets);

        // Case with no valid offsets available
        if ((valid_offsets.min()[0] > valid_offsets.max()[0]) ||
            (valid_offsets.min()[1] > valid_offsets.max()[1])) {
          m_disp_bound_image(c,r) = ZERO_SEARCH_AREA;
          ++percent_masked;
          full_search_image(c,r) = 0;
          continue;
        }

        bounds[0] = valid_offsets.min()[0];
        bounds[1] = valid_offsets.min()[1];
        bounds[2] = valid_offsets.max()[0];
        bounds[3] = valid_offsets.max()[1];
      } // End mask checking case

      m_disp_bound_image(c,r) = bounds;
      int this_area = (bounds[3]-bounds[1]+1)*(bounds[2]-bounds[0]+1);
      area += this_area;
    } // End col loop
    //right_mask_checker.advance_row();
  } // End row loop

  double num_pixels = m_disp_bound_image.rows()*m_disp_bound_image.cols();
  percent_masked  /= num_pixels;
  percent_trusted /= num_pixels;

  vw_out(DebugMessage, "stereo")
    << "Percent masked pixels  = " << percent_masked  << "\n";
  vw_out(DebugMessage, "stereo")
    << "Percent trusted prior disparities = " << percent_trusted << "\n";
  vw_out(DebugMessage, "stereo")
    << "Number missing prev disp = " << missing_prev_disp_count << "\n";

  // Next we attempt to refine the search ranges to stay within the memory usage limits.
  const int NO_CONSERVE  = 0;
  const int MAX_CONSERVE = 3;
  bool result = false;
  for (int conserve_level=NO_CONSERVE; conserve_level<=MAX_CONSERVE; ++conserve_level) {
    // Refine the disparity search ranges with the current conservation level.
    vw_out(DebugMessage, "stereo")
      << "Refining search ranges with conservation level "
      << conserve_level << "\n";
    constrain_disp_bound_image(full_search_image, prev_disparity, percent_trusted,
                               percent_masked, area, conserve_level);
    try { // Check if the computed boundaries fall within the user specified memory limit
      result = calc_main_buf_size();
      break; // Memory usage is ok!
    } catch(...) { // Used too much memory!
      vw_out(InfoMessage, "stereo")
        << "Uses too much memory with conservation level " << conserve_level << ".\n"
        << "Consider increasing --sgm-memory-limit-mb and/or decreasing the "
        << "number of threads.\n";
      result = false;
    }
  } // End loop through search range conservation attempts

  return result;
}


bool SemiGlobalMatcher::constrain_disp_bound_image(ImageView<uint8> const &full_search_image,
                                                   DisparityImage const* prev_disparity,
                                                   double percent_trusted, double percent_masked, double area,
                                                   int conserve_memory) {

  const Vector4i ZERO_SEARCH_AREA(0, 0, -1, -1);

  const BBox2i max_range_bbox(Vector2i(m_min_disp_x, m_min_disp_y), Vector2i(m_max_disp_x, m_max_disp_y));
  const double max_search_area = max_range_bbox.area();

  // Shrink the search range of full range pixels based on neighbors
  // TODO: Make this more efficient!!!
  int NEARBY_DISP_SEARCH_RANGE = 10; // Look this many pixels in each direction
  int NEARBY_DISP_EXPANSION    = 2; // Grow search range from what nearby pixels have
  if (conserve_memory == 1) // Look further, but failing pixels are discarded.
    NEARBY_DISP_SEARCH_RANGE = 25;
  if (conserve_memory == 2) // Look less and failing pixels are discarded.
    NEARBY_DISP_SEARCH_RANGE = 3;
  if (conserve_memory == 3) // Most conservative, don't look for ranges at all.
    NEARBY_DISP_SEARCH_RANGE = 0;
  double percent_shrunk = 0, shrunk_area = area;
  float conserved = 0;
  if (prev_disparity) { // No point doing this if a previous disparity image was not provided

    //// DEBUG
    //std::stringstream s;
    //s << m_disp_bound_image.cols();
    //if (conserve_memory)
    //  s << "_conserve_";
    //write_image("full_search_image"+s.str()+".tif", full_search_image);

    // Debug image to record the search size for each pixel
    //ImageView<int> search_size_image(full_search_image.cols(), full_search_image.rows());

    for (int r=0; r<m_disp_bound_image.rows(); ++r) {

      // Get vertical search range
      int min_search_r = r - NEARBY_DISP_SEARCH_RANGE;
      int max_search_r = r + NEARBY_DISP_SEARCH_RANGE;
      if (min_search_r <  0)
        min_search_r = 0;
      if (max_search_r >= m_disp_bound_image.rows())
        max_search_r = m_disp_bound_image.rows()-1;

      for (int c=0; c<m_disp_bound_image.cols(); ++c) {
        // Skip pixels without a full search range
        //Vector4i bounds = m_disp_bound_image(c,r);
        //int this_area = (bounds[3]-bounds[1]+1)*(bounds[2]-bounds[0]+1);
        if (!full_search_image(c,r)) {
          //search_size_image(c,r) = this_area;
          continue;
        }

        /* We stopped using this star pattern search because it could lead to
         * some weird looking results that were undesirable.
        // Look in a star pattern to cheaply search a larger space
        // - Keep expanding outward and quit when we find valid pixels.
        BBox2i new_range;
        for (int i=1; i<=NEARBY_DISP_SEARCH_RANGE; i++) {
          std::vector<Vector2i> coords;
          coords.push_back(Vector2i(c+i,r)); // Next set of pixels
          coords.push_back(Vector2i(c,r+i));
          coords.push_back(Vector2i(c-i,r));
          coords.push_back(Vector2i(c,r-i));
          coords.push_back(Vector2i(c-i,r-i));
          coords.push_back(Vector2i(c+i,r-i));
          coords.push_back(Vector2i(c-i,r+i));
          coords.push_back(Vector2i(c+i,r+i));

          for (size_t j=0; j<coords.size(); ++j) { // Loop through the four points.
            int cs = coords[j][0];
            int rs = coords[j][1];
            if ((cs < 0) || (rs < 0) ||
                (cs >= m_disp_bound_image.cols()) || (rs >= m_disp_bound_image.rows()))
              continue; // Skip out of bounds pixels
            if (full_search_image(cs,rs)) // Don't look at other uncertain pixels
              continue;
            Vector4i vec = m_disp_bound_image(cs,rs);
            if (vec == ZERO_SEARCH_AREA) // Skip zeroed out pixels too
              continue;
            new_range.grow(Vector2i(vec[0], vec[1])); // Expand to contain this pixel
            new_range.grow(Vector2i(vec[2], vec[3]));
            if ((c == DEBUG_X) && (r == DEBUG_Y))
              vw_out() << "Expand to contain: " << vec << "\n";
          }
          if (!new_range.empty())
            break; // Quit when we hit valid pixels
        }
        */


        // Get horizontal search range
        int min_search_c = c - NEARBY_DISP_SEARCH_RANGE;
        int max_search_c = c + NEARBY_DISP_SEARCH_RANGE;
        if (min_search_c <  0)
          min_search_c = 0;
        if (max_search_c >= m_disp_bound_image.cols())
          max_search_c = m_disp_bound_image.cols()-1;

        // Look through the search range
        BBox2i new_range;
        for (int rs=min_search_r; rs<=max_search_r; ++rs) {
          for (int cs=min_search_c; cs<=max_search_c; ++cs) {
            if (full_search_image(cs,rs)) // Don't look at other uncertain pixels
              continue;
            // Expand bounding box
            Vector4i vec = m_disp_bound_image(cs,rs);
            if (vec == ZERO_SEARCH_AREA) // Skip zeroed out pixels too
              continue;
            new_range.grow(Vector2i(vec[0], vec[1]));
            new_range.grow(Vector2i(vec[2], vec[3]));
          }
        }

        if (new_range.empty()) { // If we did not find a new estimate
          // If worried about memory, don't try to solve pixels with no estimate.
          if (conserve_memory > 0) {
            m_disp_bound_image(c,r) = ZERO_SEARCH_AREA;
            percent_shrunk += 1.0;
            shrunk_area -= max_search_area;
            conserved += 1.0;
            //search_size_image(c,r) = 0;
          }
          //else
          //  search_size_image(c,r) = this_area;

          // Otherwise use the full search range for them.
          continue;
        }
        // Grow the bounding box a bit and then record it
        new_range.expand(NEARBY_DISP_EXPANSION);
        new_range.crop(max_range_bbox); // Constrain to global limits
        m_disp_bound_image(c,r) = Vector4i(new_range.min().x(),   new_range.min().y(),
                                           new_range.max().x(),   new_range.max().y());
        percent_shrunk += 1.0;
        shrunk_area -= (max_search_area - new_range.area());
        //search_size_image(c,r) = new_range.area();

      } // End col loop
    } // End row loop
    //write_image("search_size_image"+s.str()+".tif", search_size_image);
  } // End of search range shrinking code

  // Compute some statistics for help improving the speed
  double num_pixels = m_disp_bound_image.rows()*m_disp_bound_image.cols();
  area            /= num_pixels;
  shrunk_area     /= num_pixels;
  percent_shrunk  /= num_pixels;
  conserved       /= num_pixels;

  double initial_percent_full_range = 1.0 - (percent_trusted+percent_masked);
  double final_percent_full_range   = initial_percent_full_range - percent_shrunk;
  vw_out(DebugMessage, "stereo") << "Max pixel search area  = " << max_search_area << "\n";
  vw_out(DebugMessage, "stereo") << "Mean pixel search area (initial) = " << area << "\n";
  vw_out(DebugMessage, "stereo") << "Mean pixel search area (final) = "
    << shrunk_area << "\n";
  vw_out(DebugMessage, "stereo") << "Percent shrunk pixels = " << percent_shrunk << "\n";
  vw_out(DebugMessage, "stereo") << "Percent full search range pixels (initial) = "
    << initial_percent_full_range << "\n";
  vw_out(DebugMessage, "stereo") << "Percent full search range pixels (final) = "
    << final_percent_full_range   << "\n";
  if (conserve_memory)
    vw_out(DebugMessage, "stereo") << "Percent pixels conserved for memory = "
      << conserved << "\n";

  // Return false if the image cannot be processed
  if ((area <= 0) || (percent_masked >= 100))
    return false;
  return true;

} // End populate_disp_bound_image

// Calculate the size of the main buffer. This will be used in
// populate_disp_bound_image() to shrink the search ranges if the total
// buffer size is more than what the user expects.
size_t SemiGlobalMatcher::calc_main_buf_size() {

  // Init the starts data storage
  m_buffer_starts.set_size(m_num_output_cols, m_num_output_rows);

  vw_out(DebugMessage, "stereo") << "SGM: Num pixels = "      << m_num_output_rows * m_num_output_cols << "\n";
  vw_out(DebugMessage, "stereo") << "SGM: Num disparities = " << m_num_disp << "\n";

  // For each pixel, record the starting offset and add the disparity search area
  //  of this pixel to the running offset total.
  size_t main_buf_size = 0;
  for (int r=0; r<m_num_output_rows; ++r) {
    for (int c=0; c<m_num_output_cols; ++c) {
      m_buffer_starts(c,r) = main_buf_size;
      main_buf_size += get_num_disparities(c, r);
    }
  }
  // Finished computing the pixel offsets.

  // It is not clear why this is needed
  if (main_buf_size < 6)
    main_buf_size = 6;

  // Record this value early, as it is used right below when calling one_buf_size()
  // and/or multi_buf_size().
  m_main_buf_size = main_buf_size;

  // Verify that allocating the "small" buffers won't put us over the limit.
  // - m_main_buf_size must be set for this to work.
  const int PATHS_PER_PASS  = 1;
  const int NUM_MGM_BUFFERS = 4;
  size_t small_buf_size = 0;
  if (m_use_mgm) {
    // Four vertical buffers and four horizontal buffers.
    small_buf_size
      = MultiAccumRowBuffer::multi_buf_size(this, PATHS_PER_PASS, true) * NUM_MGM_BUFFERS
      + MultiAccumRowBuffer::multi_buf_size(this, PATHS_PER_PASS, false) * NUM_MGM_BUFFERS;
  } else {
    int num_threads   = vw_settings().default_num_threads();
    small_buf_size = OneLineBuffer::one_buf_size(this) * num_threads;
  }

  double main_buf_MB  = double(main_buf_size) * MainBufToMB;
  double small_buf_MB = double(small_buf_size) * SmallBufToMB;
  double total_buf_MB = main_buf_MB + small_buf_MB;

  if (total_buf_MB > m_memory_limit_mb)
    vw_throw(ArgumentErr()
             << "SGM: Required memory usage is " << total_buf_MB
             << " MB which is greater than the cap of --corr-memory-mb "
             << m_memory_limit_mb << " MB. Consider increasing this or "
             << "reducing the number of threads.\n");

  return main_buf_size;
}

void SemiGlobalMatcher::allocate_large_buffers() {

  //Timer timer_total("Memory allocation");

  const size_t BYTES_PER_MB = 1024*1024;
  const size_t main_buf_size          = calc_main_buf_size();
  const size_t cost_buffer_num_bytes  = main_buf_size * sizeof(CostType);
  const size_t accum_buffer_num_bytes = main_buf_size * sizeof(AccumCostType);

  vw_out(DebugMessage, "stereo") << "SGM: Allocating buffer of size: " << cost_buffer_num_bytes/BYTES_PER_MB << " MB\n";

  m_cost_buffer.reset(new CostType[main_buf_size]);

  vw_out(DebugMessage, "stereo") << "SGM: Allocating buffer of size: " << accum_buffer_num_bytes/BYTES_PER_MB << " MB\n";

  // Allocate the requested memory and init all to zero
  m_accum_buffer.reset(new AccumCostType[main_buf_size]);
  memset(m_accum_buffer.get(), 0, accum_buffer_num_bytes);
}



void SemiGlobalMatcher::populate_adjacent_disp_lookup_table() {

  const int TABLE_WIDTH = 8;
  m_adjacent_disp_lookup.resize(m_num_disp*TABLE_WIDTH);

  // Loop through the disparities
  int d = 0;
  for (int dy=m_min_disp_y; dy<=m_max_disp_y; ++dy) {
    // Figure out above and below disparities with bounds checking
    int y_less = dy - 1;
    int y_more = dy + 1;
    if (y_less < m_min_disp_y) y_less = dy;
    if (y_more > m_max_disp_y) y_more = dy;

    int y_less_o = y_less - m_min_disp_y; // Offset from min y disparity
    int y_o      = dy     - m_min_disp_y;
    int y_more_o = y_more - m_min_disp_y;

    for (int dx=m_min_disp_x; dx<=m_max_disp_x; ++dx) {

      // Figure out left and right disparities with bounds checking
      int x_less = dx - 1;
      int x_more = dx + 1;
      if (x_less < m_min_disp_x) x_less = dx;
      if (x_more > m_max_disp_x) x_more = dx;

      int x_less_o = x_less - m_min_disp_x; // Offset from min x disparity
      int x_o      = dx     - m_min_disp_x;
      int x_more_o = x_more - m_min_disp_x;

      // Record the disparity indices of each of the adjacent pixels
      int table_pos = d*TABLE_WIDTH;
      m_adjacent_disp_lookup[table_pos+0] = y_less_o*m_num_disp_x + x_o;      // The four adjacent pixels
      m_adjacent_disp_lookup[table_pos+1] = y_o     *m_num_disp_x + x_less_o;
      m_adjacent_disp_lookup[table_pos+2] = y_o     *m_num_disp_x + x_more_o;
      m_adjacent_disp_lookup[table_pos+3] = y_more_o*m_num_disp_x + x_o;
      m_adjacent_disp_lookup[table_pos+4] = y_less_o*m_num_disp_x + x_less_o; // The four diagonal pixels
      m_adjacent_disp_lookup[table_pos+5] = y_less_o*m_num_disp_x + x_more_o;
      m_adjacent_disp_lookup[table_pos+6] = y_more_o*m_num_disp_x + x_less_o;
      m_adjacent_disp_lookup[table_pos+7] = y_more_o*m_num_disp_x + x_more_o;

      ++d;
    }
  }
} // End function populate_adjacent_disp_lookup_table


#if not defined(VW_ENABLE_SSE) || (VW_ENABLE_SSE==0)
// Note: local and output are the same size.
// full_prior_buffer is always length m_num_disps and comes in initialized to a
//  large flag value.  When the function quits the buffer must be returned to this state.
void SemiGlobalMatcher::evaluate_path(int col, int row, int col_p, int row_p,
                       AccumCostType* const prior,
                       AccumCostType*       full_prior_buffer,
                       CostType     * const local,
                       AccumCostType*       output,
                       int path_intensity_gradient, bool debug) {

  // Decrease p2 (jump cost) with increasing disparity along the path
  AccumCostType p2_mod = m_p2;
  if (path_intensity_gradient > 0)
    p2_mod /= path_intensity_gradient;
  if (p2_mod < m_p1)
    p2_mod = m_p1;

  //int num_disparities   = get_num_disparities(col,   row); // Can be input arg
  //int num_disparities_p = get_num_disparities(col_p, row_p);

  Vector4i pixel_disp_bounds   = m_disp_bound_image(col, row);
  Vector4i pixel_disp_bounds_p = m_disp_bound_image(col_p, row_p);

  // Init the min prior in case the previous pixel is invalid.
  AccumCostType BAD_VAL = get_bad_accum_val();
  AccumCostType min_prior = BAD_VAL;

  // Insert the valid disparity scores into full_prior buffer so they are
  //  easy to access quickly within the pixel loop below.
  int d = 0;
  for (int dy=pixel_disp_bounds_p[1]; dy<=pixel_disp_bounds_p[3]; ++dy) {

    // Get initial fill linear storage index for this dy row
    int full_index = xy_to_disp(pixel_disp_bounds_p[0], dy);

    for (int dx=pixel_disp_bounds_p[0]; dx<=pixel_disp_bounds_p[2]; ++dx) {

      // Get the min prior while we are at it.
      if (prior[d] < min_prior) {
        min_prior  = prior[d];
      }

      full_prior_buffer[full_index] = prior[d];
      ++full_index;
      ++d;
    }
  }
  AccumCostType min_prev_disparity_cost = min_prior + p2_mod;
  const int LOOKUP_TABLE_WIDTH = 8;

  // Loop through disparities for this pixel
  int packed_d = 0; // Index for cost and output vectors
  for (int dy=pixel_disp_bounds[1]; dy<=pixel_disp_bounds[3]; ++dy) {

    // Need the disparity index from all of m_num_disp for proper indexing into full_prior_buffer
    int full_d = xy_to_disp(pixel_disp_bounds[0], dy);

    for (int dx=pixel_disp_bounds[0]; dx<=pixel_disp_bounds[2]; ++dx) {

      // Start with the cost for the same disparity in the previous pixel
      AccumCostType lowest_combined_cost = full_prior_buffer[full_d];

      // Compare to the eight adjacent disparities using the lookup table
      const int lookup_index = full_d*LOOKUP_TABLE_WIDTH;

      // TODO: This is the slowest part of the algorithm!
      // Note that the lookup table indexes into a full size buffer of disparities, not the compressed
      //  buffers that are stored for each pixel.  This allows us to use a single lookup table for every pixel
      //  and avoid any bounds checking logic inside this loop.
      AccumCostType lowest_adjacent_cost =                  full_prior_buffer[m_adjacent_disp_lookup[lookup_index  ]];
      lowest_adjacent_cost = std::min(lowest_adjacent_cost, full_prior_buffer[m_adjacent_disp_lookup[lookup_index+1]]);
      lowest_adjacent_cost = std::min(lowest_adjacent_cost, full_prior_buffer[m_adjacent_disp_lookup[lookup_index+2]]);
      lowest_adjacent_cost = std::min(lowest_adjacent_cost, full_prior_buffer[m_adjacent_disp_lookup[lookup_index+3]]);
      lowest_adjacent_cost = std::min(lowest_adjacent_cost, full_prior_buffer[m_adjacent_disp_lookup[lookup_index+4]]);
      lowest_adjacent_cost = std::min(lowest_adjacent_cost, full_prior_buffer[m_adjacent_disp_lookup[lookup_index+5]]);
      lowest_adjacent_cost = std::min(lowest_adjacent_cost, full_prior_buffer[m_adjacent_disp_lookup[lookup_index+6]]);
      lowest_adjacent_cost = std::min(lowest_adjacent_cost, full_prior_buffer[m_adjacent_disp_lookup[lookup_index+7]]);

      // Now add the adjacent penalty cost and compare to the local cost
      lowest_adjacent_cost += m_p1;
      lowest_combined_cost = std::min(lowest_combined_cost, lowest_adjacent_cost);

      // Compare to the lowest prev disparity cost regardless of location
      lowest_combined_cost = std::min(lowest_combined_cost, min_prev_disparity_cost);

      // The output cost = local cost + lowest combined cost - min_prior
      // - Subtracting out min_prior avoids overflow.
      output[packed_d] = local[packed_d] + lowest_combined_cost - min_prior;

      //if (debug) {
      //  printf("Details %d: local = %d, lowest_adjacent_cost = %d, prev_cost = %d, lowest_combined_cost = %d, output = %d\n",
      //        packed_d, local[packed_d], lowest_adjacent_cost, full_prior_buffer[full_d], lowest_combined_cost, output[packed_d]);
      //}

      ++packed_d;
      ++full_d;
    }
  } // End loop through this disparity

  if (debug) {
    int min_val   = 99999;
    int min_index = 0;

    int i=0;
    for (int dy=pixel_disp_bounds[1]; dy<=pixel_disp_bounds[3]; ++dy) {
      for (int dx=pixel_disp_bounds[0]; dx<=pixel_disp_bounds[2]; ++dx) {
        if (output[i] < min_val) {
          min_val   = output[i];
          min_index = i;
        }
        ++i;
      }
    }

    DisparityType dx, dy;
    disp_to_xy(min_index, dx, dy);
  }
  // Remove the valid disparity scores from full_prior buffer.
  for (int dy=pixel_disp_bounds_p[1]; dy<=pixel_disp_bounds_p[3]; ++dy) {

    // Get initial fill linear storage index for this dy row
    int full_index = xy_to_disp(pixel_disp_bounds_p[0], dy);

    for (int dx=pixel_disp_bounds_p[0]; dx<=pixel_disp_bounds_p[2]; ++dx) {

      full_prior_buffer[full_index] = BAD_VAL;
      ++full_index;
    }
  }

}
#endif

#if defined(VW_ENABLE_SSE) && (VW_ENABLE_SSE==1)
void compute_path_internals_sse(uint16* dL, uint16* d0, uint16* d1,
                                uint16* d2, uint16* d3,
                                uint16* d4, uint16* d5, uint16* d6,
                                uint16* d7, uint16* d8,
                                __m128i& _dJ, __m128i& _dP, __m128i& _dp1, uint16* dRes,
                                int sse_index, int &output_index,
                                SemiGlobalMatcher::AccumCostType* output) {
  // Load data from arrays into SSE registers
  __m128i _dL = _mm_load_si128((__m128i*) dL);
  __m128i _d0 = _mm_load_si128((__m128i*) d0);
  __m128i _d1 = _mm_load_si128((__m128i*) d1);
  __m128i _d2 = _mm_load_si128((__m128i*) d2);
  __m128i _d3 = _mm_load_si128((__m128i*) d3);
  __m128i _d4 = _mm_load_si128((__m128i*) d4);
  __m128i _d5 = _mm_load_si128((__m128i*) d5);
  __m128i _d6 = _mm_load_si128((__m128i*) d6);
  __m128i _d7 = _mm_load_si128((__m128i*) d7);
  __m128i _d8 = _mm_load_si128((__m128i*) d8);

  // Operation = min(min(d1...d8)+dp1, d0, dJ) + dL - dP

  // Start computing the min
  __m128i _min12   = _mm_min_epu16(_d1, _d2);
  __m128i _min34   = _mm_min_epu16(_d3, _d4);
  __m128i _min56   = _mm_min_epu16(_d5, _d6);
  __m128i _min78   = _mm_min_epu16(_d7, _d8);
  // Keep computing the min
  __m128i _min1234 = _mm_min_epu16(_min12, _min34);
  __m128i _min5678 = _mm_min_epu16(_min56, _min78);
  // Finish computing the min
  __m128i _minAdj = _mm_min_epu16(_min1234, _min5678);
  __m128i _minO   = _mm_min_epu16(_d0, _dJ);

  // Perform the required computations
  __m128i _result = _mm_adds_epu16(_minAdj, _dp1);
  _result = _mm_min_epu16(_result, _minO);
  _result = _mm_adds_epu16(_result, _dL);
  _result = _mm_subs_epu16(_result, _dP);

  // Fetch results from the output register
  _mm_store_si128((__m128i*) dRes, _result);

  // Copy the valid results from the register.
  for (int i = 0; i < sse_index; i++) {
    output[output_index++] = dRes[i];
  }
} // end function compute_path_internals_sse
#endif

/// Non-sse backup for compute_path_internals_sse
void compute_path_internals(uint16* dL, uint16* d0, uint16* d1, uint16* d2, uint16* d3,
                            uint16* d4, uint16* d5, uint16* d6, uint16* d7, uint16* d8,
                            SemiGlobalMatcher::AccumCostType dJ,
                            SemiGlobalMatcher::AccumCostType dP,
                            SemiGlobalMatcher::AccumCostType dp1,
                            uint16* dRes, int sse_index, int &output_index,
                            SemiGlobalMatcher::AccumCostType* output) {
  // Operation = min(min(d1...d8)+dp1, d0, dJ) + dL - dP

  for (int i = 0; i < sse_index; i++) {

    uint16 minAdj = std::min(d1[i], d2[i]);
    minAdj = std::min(minAdj, d3[i]);
    minAdj = std::min(minAdj, d4[i]);
    minAdj = std::min(minAdj, d5[i]);
    minAdj = std::min(minAdj, d6[i]);
    minAdj = std::min(minAdj, d7[i]);
    minAdj = std::min(minAdj, d8[i]);
    minAdj += dp1;

    uint16 minVal = std::min(minAdj, d0[i]);
           minVal = std::min(minVal, dJ);
    output[output_index++] = minVal + (dL[i] - dP);
  }
} // end function compute_path_internals

#if defined(VW_ENABLE_SSE) && (VW_ENABLE_SSE==1)
void SemiGlobalMatcher::evaluate_path(int col, int row, int col_p, int row_p,
                       AccumCostType* const prior,
                       AccumCostType*       full_prior_buffer,
                       CostType     * const local,
                       AccumCostType*       output,
                       int path_intensity_gradient, bool debug) {

  // Decrease p2 (jump cost) with increasing disparity along the path
  AccumCostType p2_mod = m_p2;
  if (path_intensity_gradient > 0)
    p2_mod /= path_intensity_gradient;
  if (p2_mod < m_p1)
    p2_mod = m_p1;

  Vector4i pixel_disp_bounds   = m_disp_bound_image(col, row);
  Vector4i pixel_disp_bounds_p = m_disp_bound_image(col_p, row_p);

  // Init the min prior in case the previous pixel is invalid.
  AccumCostType BAD_VAL = get_bad_accum_val();
  AccumCostType min_prior = BAD_VAL;

  // Insert the valid disparity scores into full_prior buffer so they are
  //  easy to access quickly within the pixel loop below.
  // - If we don't use a full sized buffer, our adjacent disparity lookup
  //   table could not be used!
  int d = 0;
  for (int dy=pixel_disp_bounds_p[1]; dy<=pixel_disp_bounds_p[3]; ++dy) {

    // Get initial fill linear storage index for this dy row
    int full_index = xy_to_disp(pixel_disp_bounds_p[0], dy);

    for (int dx=pixel_disp_bounds_p[0]; dx<=pixel_disp_bounds_p[2]; ++dx) {

      // Get the min prior while we are at it.
      if (prior[d] < min_prior) {
        min_prior  = prior[d];
      }

      full_prior_buffer[full_index] = prior[d];
      ++full_index;
      ++d;
    }
  }
  AccumCostType min_prev_disparity_cost = min_prior + p2_mod;

  const int LOOKUP_TABLE_WIDTH = 8;

  // Allocate linear storage for data to pass to SSE instructions
  const int SSE_BUFF_LEN = 8;
  uint16 d_packed[SSE_BUFF_LEN*11] __attribute__ ((aligned (16))); // TODO: Could be passed in!
  uint16* dL   = &(d_packed[0*SSE_BUFF_LEN]);
  uint16* d0   = &(d_packed[1*SSE_BUFF_LEN]);
  uint16* d1   = &(d_packed[2*SSE_BUFF_LEN]);
  uint16* d2   = &(d_packed[3*SSE_BUFF_LEN]);
  uint16* d3   = &(d_packed[4*SSE_BUFF_LEN]);
  uint16* d4   = &(d_packed[5*SSE_BUFF_LEN]);
  uint16* d5   = &(d_packed[6*SSE_BUFF_LEN]);
  uint16* d6   = &(d_packed[7*SSE_BUFF_LEN]);
  uint16* d7   = &(d_packed[8*SSE_BUFF_LEN]);
  uint16* d8   = &(d_packed[9*SSE_BUFF_LEN]);
  uint16* dRes = &(d_packed[10*SSE_BUFF_LEN]); // The results

  // Set up constant SSE registers that never change
  __m128i _dJ  = _mm_set1_epi16(static_cast<int16>(min_prev_disparity_cost));
  __m128i _dP  = _mm_set1_epi16(static_cast<int16>(min_prior));
  __m128i _dp1 = _mm_set1_epi16(static_cast<int16>(m_p1));

  //printf("dJ = %d, dP = %d, dp1 = %d\n", min_prev_disparity_cost, min_prior, m_p1);

  // Loop through disparities for this pixel
  int sse_index = 0, output_index = 0;
  int packed_d = 0; // Index for cost and output vectors
  for (int dy=pixel_disp_bounds[1]; dy<=pixel_disp_bounds[3]; ++dy) {

    // Need the disparity index from all of m_num_disp for proper indexing into full_prior_buffer
    int full_d = xy_to_disp(pixel_disp_bounds[0], dy);

    for (int dx=pixel_disp_bounds[0]; dx<=pixel_disp_bounds[2]; ++dx) {

      // Get local value and matching disparity value
      dL[sse_index] = local[packed_d];
      d0[sse_index] = full_prior_buffer[full_d];

      // Get the 8 surrounding values.
      // - Is there any way to speed this up?
      const int lookup_index = full_d*LOOKUP_TABLE_WIDTH;
      d1[sse_index] = full_prior_buffer[m_adjacent_disp_lookup[lookup_index  ]];
      d2[sse_index] = full_prior_buffer[m_adjacent_disp_lookup[lookup_index+1]];
      d3[sse_index] = full_prior_buffer[m_adjacent_disp_lookup[lookup_index+2]];
      d4[sse_index] = full_prior_buffer[m_adjacent_disp_lookup[lookup_index+3]];
      d5[sse_index] = full_prior_buffer[m_adjacent_disp_lookup[lookup_index+4]];
      d6[sse_index] = full_prior_buffer[m_adjacent_disp_lookup[lookup_index+5]];
      d7[sse_index] = full_prior_buffer[m_adjacent_disp_lookup[lookup_index+6]];
      d8[sse_index] = full_prior_buffer[m_adjacent_disp_lookup[lookup_index+7]];

      ++packed_d;
      ++full_d;
      ++sse_index;

      // Keep packing the SSE buffers until they are filled up, then use SSE to operate on
      // all of the data at once.
      if (sse_index == SSE_BUFF_LEN) {

        compute_path_internals_sse(dL, d0, d1, d2, d3, d4, d5, d6, d7, d8,
                                   _dJ, _dP, _dp1, dRes, sse_index, output_index, output);
        //compute_path_internals(dL, d0, d1, d2, d3, d4, d5, d6, d7, d8,
        //                   min_prev_disparity_cost, min_prior, m_p1, dRes, sse_index, output_index, output);

        sse_index = 0;
      } // End SSE operations

    }
  } // End loop through this disparity

  // If there is data left over in the buffer, process it now.
  if (sse_index > 0) {
    compute_path_internals_sse(dL, d0, d1, d2, d3, d4, d5, d6, d7, d8,
                           _dJ, _dP, _dp1, dRes, sse_index, output_index, output);
    //compute_path_internals(dL, d0, d1, d2, d3, d4, d5, d6, d7, d8,
    //                       min_prev_disparity_cost, min_prior, m_p1, dRes, sse_index, output_index, output);
  }

  // Remove the valid disparity scores from full_prior buffer.
  for (int dy=pixel_disp_bounds_p[1]; dy<=pixel_disp_bounds_p[3]; ++dy) {

    // Get initial fill linear storage index for this dy row
    int full_index = xy_to_disp(pixel_disp_bounds_p[0], dy);

    for (int dx=pixel_disp_bounds_p[0]; dx<=pixel_disp_bounds_p[2]; ++dx) {

      full_prior_buffer[full_index] = BAD_VAL;
      ++full_index;
    }
  }

} // End evaluate_path SSE
#endif


/* This function is not 100% successful at removing "multiple minimums"
   but it usually works.  Multiple minimums are most commonly
   encountered over large regions of low quality image where most
   costs are zero.  The result is very poor subpixel performance,
   oscillating disparity values, and other artifacts.
*/
int SemiGlobalMatcher::select_best_disparity(AccumCostType * accum_vec,
                                             Vector4i const& bounds,
                                             int &min_index,
                                             std::vector<AccumCostType> & buffer,
                                             bool debug) {
  // The input cost buffer is always a rectangle
  int height   = (bounds[3] - bounds[1] + 1);
  int width    = (bounds[2] - bounds[0] + 1);
  int num_vals = height*width;

  // Find and count the minimum value and create a copy of the input cost buffer
  buffer.resize(num_vals);
  int index     = 0;
  int min_count = 0;
  min_index = 0;
  AccumCostType min_val = std::numeric_limits<AccumCostType>::max();
  AccumCostType value;
  for (int row=0; row<height; ++row) {
    for (int col=0; col<width; ++col) {
      value = accum_vec[index]; // Copy value
      buffer[index] = value;

      // Track the minimum value
      if (value == min_val)
        ++min_count;
      if (value < min_val) {
        min_index = index;
        min_val = value;
        min_count = 1;
      }
      ++index;
    }
  } // End double loop

  // Filtering steps below will go back and forth between buffers.
  AccumCostType * input_array  = accum_vec;
  AccumCostType * output_array = &(buffer[0]);
  AccumCostType * swap = 0;

  // If there is more than one value tied for minimum cost,
  //  apply smoothing kernels to the accumulation values until
  //  we are reduced to a single minimum or we give up.

  //TODO: Function calls to do this!  The VW convolve functions are messy.
  const int FILTER_WIDTH = 3;
  std::vector<double> filter(FILTER_WIDTH);
  filter[0] = 1.0/3.0; filter[1] = 1.0/3.0; filter[2] = 1.0/3.0;
  int filter_radius = (FILTER_WIDTH -1) / 2;

  const int MAX_ITERATIONS = 6; // This many iterations total
  const int VERT_ITERATION = 5; // After this many iterations switch to vertical filtering
  int iter_count = 0;
  while (min_count > 1) {

    // Swap buffers
    // - The first swap means we go from buffer to accum_vec
    swap         = input_array;
    input_array  = output_array;
    output_array = swap;

    index     = 0;
    min_count = 0;
    min_val = std::numeric_limits<AccumCostType>::max();
    min_index = 0;
    for (int row=0; row<height; ++row) {
      // Convolve this row
      for (int col=0; col<width; ++col) {
        int min = -1*filter_radius;
        int max =    filter_radius;
        double result      = 0;
        double weight_total = 0;

        // Usually horizontal filtering is good enough but sometimes vertical
        //  filtering is needed to get down to one minimum.
        if (iter_count < VERT_ITERATION) {
          // Horizontal filtering
          if (min + col < 0) min = 0; // TODO: Hard coded to radius 1!!
          if (max + col >= width) max = 0;

          for (int k=min; k<=max; ++k) {
            double weight = filter[k+filter_radius];
            result       += static_cast<double>(input_array[index + k]) * weight;
            weight_total += weight;
          }
        } else {
          // Vertical filtering
          if (min + row < 0) min = 0; // TODO: Hard coded to radius 1!!
          if (max + row >= height) max = 0;

          for (int k=min; k<=max; ++k) {
            double weight = filter[k+filter_radius];
            result       += static_cast<double>(input_array[index + k*width]) * weight;
            weight_total += weight;
          }
        }
        value = static_cast<AccumCostType>(round(result / weight_total));

        // Track the minimum value and the minimum count
        if (value == min_val)
          ++min_count;
        if (value < min_val) {
          min_index = index;
          min_val = value;
          min_count = 1;
        }
        output_array[index] = value;

        ++index;
      } // End col loop
    } // End row loop

    iter_count++;

    // Make sure this never becomes an infinite loop
    if (iter_count >= MAX_ITERATIONS)
      break;
  } // End smoothing iteration loop

  // Copy result back to the accum vector if needed
  if ((iter_count > 0) && (iter_count % 2 == 0)) {
    for (int i = 0; i < index; i++)
      input_array[i] = output_array[i];
  }

  return min_count;
} // End function select_best_disparity

SemiGlobalMatcher::DisparityImage
SemiGlobalMatcher::create_disparity_view() {
  // Init output vector
  DisparityImage disparity(m_num_output_cols, m_num_output_rows);

  // For each element in the accumulated costs matrix,
  //  select the disparity with the lowest accumulated cost.
  //Timer timer("Calculate Disparity Minimum");

  /* These are for experimenting with ways to evaluate the
  ImageView<int> wto_intervals   (m_num_output_cols, m_num_output_rows);
  ImageView<int> second_intervals(m_num_output_cols, m_num_output_rows);
  ImageView<int> max_intervals   (m_num_output_cols, m_num_output_rows);
  ImageView<int> min_counts      (m_num_output_cols, m_num_output_rows);
  ImageView<int> wto_ratio       (m_num_output_cols, m_num_output_rows);
  ImageView<int> raw             (m_num_output_cols, m_num_output_rows);
  ImageView<int> cost_second     (m_num_output_cols, m_num_output_rows);
  ImageView<int> cost_ratio      (m_num_output_cols, m_num_output_rows);
  ImageView<int> best_costs      (m_num_output_cols, m_num_output_rows);
  ImageView<int> worst_costs     (m_num_output_cols, m_num_output_rows);
  double mean;
  AccumCostType best, second;
  int count = 0, worst=0;
  */

  DisparityType dx, dy;
  int min_index=0;
  std::vector<AccumCostType> accum_buffer;
  for (int j = 0; j < m_num_output_rows; j++) {
    for (int i = 0; i < m_num_output_cols; i++) {

      int num_disp = get_num_disparities(i, j);

      if (num_disp == 0) {
        // Pixels with no search area were never valid.
        disparity(i,j) = DisparityImage::pixel_type();
        invalidate(disparity(i,j));

        // DEBUG
        /*wto_intervals(i,j) = 0;
        max_intervals(i,j) = 0;
        second_intervals(i,j) = 0;
        min_counts(i,j) = 0;
        wto_ratio(i,j) = 0;
        raw(i,j) = 0;
        cost_second(i,j) = 0;
        cost_ratio(i,j) = 0;
        best_costs(i,j) = 0;
        worst_costs(i,j) = 0;*/
        continue;
      }

      // Valid pixel, choose the best disparity.

      const Vector4i bounds = m_disp_bound_image(i, j);
      AccumCostType  *accum_vec = get_accum_vector(i, j);
      bool debug = false;
      select_best_disparity(accum_vec, bounds, min_index, accum_buffer, debug);
      disp_index_to_xy(min_index, i, j, dx, dy);
      disparity(i,j) = DisparityImage::pixel_type(dx, dy);

      // DEBUG
      /*
      // Cost measurements
      CostType min_cost = 255, second_cost=0, worst_cost = 0;
      double mean_cost = 0;
      CostType *cost_vec = get_cost_vector(i,j);
      for (size_t n=0; n<num_disp; ++n) {
        mean_cost += cost_vec[n];
        if (cost_vec[n] < min_cost) {
          second_cost = min_cost;
          min_cost = cost_vec[n];
        } else
          if (cost_vec[n] < second_cost)
            second_cost = cost_vec[n];
        if (cost_vec[n] > worst_cost)
          worst_cost = cost_vec[n];
      }
      mean_cost /= static_cast<double>(num_disp);
      cost_second(i,j) = second_cost - min_cost;
      if (min_cost > 0)
        cost_ratio(i,j) = floor(mean_cost /  (double)min_cost);
      else
        cost_ratio(i,j) = 0;
      best_costs (i,j) = min_cost;
      worst_costs(i,j) = worst_cost;

      // Accumulation measurements
      wto_intervals(i,j) = floor(mean - (double)best);
      second_intervals(i,j) = second - best;
      max_intervals(i,j) = worst - best;
      min_counts(i,j) = count;
      if (best > 0)
        wto_ratio(i,j) = floor(mean /  (double)best);
      else
        wto_ratio(i,j) = 0;
      raw(i,j) = best;
      */

      /*
      if ((i >= 15) && (i <= 15) && (j==40)) { // DEBUG!
        printf("ACC costs (%d,%d): %d, %d\n", i, j, dx, dy);

        AccumCostType* vec = get_accum_vector(i, j);
        const int num_disp = get_num_disparities(i, j);
      }
      */
    }
  }
  // write_image("wta_intervals.tif",    wto_intervals);
  // write_image("second_intervals.tif", second_intervals);
  // write_image("max_intervals.tif",    max_intervals);
  // write_image("min_counts.tif",       min_counts);
  // write_image("wto_ratios.tif",       wto_ratio);
  // write_image("raw_scores.tif",       raw);
  // write_image("cost_second.tif",      cost_second);
  // write_image("cost_ratio.tif",       cost_ratio);
  // write_image("best_costs.tif",       best_costs);
  // write_image("worst_costs.tif",      worst_costs);

  vw_out(DebugMessage, "stereo") << "Finished creating integer disparity image.\n";
  return disparity;
}

// A number of proposed subpixel algorithms
double linearFit(double x) {
  return x/2.0;
}
double parabolaFit(double x) { // Use our 2d implementation instead of this one
  return x/(x+1.0);
}
double poly4Fit(double x) {
  return (x*x*x*x + x)/4.0;
}
double sinFit(double x) {
  const double PI = 3.14159265359; // TODO: MOVE THIS
  return 0.5 * (sin(x*PI/2.0 - PI/2.0) + 1.0);
}
double cosFit(double x) {
  const double PI = 3.14159265359; // TODO: MOVE THIS
  return (1 - cos(x*PI/3.0));
}
double alg16(double x) {
  return std::max(cosFit(x), poly4Fit(x));
}
double lcBlendFit(double x) {
  const double PI = 3.14159265359; // TODO: MOVE THIS

  double factor = 1.195 - cos(x*(PI/2.3));
  return cosFit(x)*factor + linearFit(x)*(1.0-factor);
}

/// Crude way of performing subpixel interpolation when only two values are available
/// - Returns fraction of distance from the primary to the other value.
double SemiGlobalMatcher::two_value_subpixel(AccumCostType primary, AccumCostType other) {
  return 0.5*(static_cast<double>(primary) / static_cast<double>(other));
}

/// Add this offset to the integer disparity to get the final result.
double SemiGlobalMatcher::compute_subpixel_offset(AccumCostType prev, AccumCostType center,
                                                  AccumCostType next,
                                                  bool left_bound, bool right_bound,
                                                  bool debug) {
  double ld = prev - center; // Left  diff
  double rd = next - center; // Right diff
  if ((rd == 0) && (ld == 0)) // Handle case where all values are equal or both bounds are set
    return 0;

  // If only two values are available, use a lower quality two value interpolation.
  if (left_bound)
    return two_value_subpixel(center, next); // Shift right
  if (right_bound)
    return -1.0*two_value_subpixel(center, prev); // Shift left

  // Pick which direction to interpolate in.
  // - Never interpolate in a bounded direction.
  double x = rd/ld;
  double mult = -1.0;
  if (ld < rd) {
    x    = ld/rd;
    mult = 1.0;
  }

  // Use the selected subpixel function
  double value = 0;
  switch(m_subpixel_type) {
    case SUBPIXEL_POLY4:    value = poly4Fit  (x); break;
    case SUBPIXEL_COSINE:   value = cosFit    (x); break;
    case SUBPIXEL_LC_BLEND: value = lcBlendFit(x); break;
    default:                value = linearFit (x); break;
  };
  // Complete computation
  return (value - 0.5)*mult;
}
/*
double SemiGlobalMatcher::compute_subpixel_ratio(AccumCostType prev, AccumCostType center, AccumCostType next) {
  // Just return the ratio that would be used in compute_subpixel_offset
  // - If this value is written to disk, subpixel functions can be experimented with Matlab or something.
  double ld = prev - center;
  double rd = next - center;
  if ((rd == 0) && (ld == 0))
    return 0;
  if (ld < rd)
    return ld/rd;
  else
    return rd/ld;
}
*/


// Test out alternate subpixel methods
ImageView<PixelMask<Vector2f> > SemiGlobalMatcher::
create_disparity_view_subpixel(DisparityImage const& integer_disparity) {

  //Timer timer("Calculate Subpixel Disparity");

  typedef  PixelMask<Vector2f> p_type;
  ImageView<p_type> disparity(m_num_output_cols, m_num_output_rows);

  ParabolaFit2d fitter; // Only used with parabola2d

  vw_out(DebugMessage, "stereo") << "Creating subpixel disparity image...\n";

  // DEBUG
  //math::Histogram hist_dx(201, -1.0, 1.0), hist_dy(201, -1.0, 1.0);
  //std::ofstream rawFile("raw.csv");

  // For each element in the accumulated costs matrix,
  //  select the disparity with the lowest accumulated cost.
  double percent_bad = 0;
  double delta_x, delta_y;
  for (int j = 0; j < m_num_output_rows; j++) {
    for (int i = 0; i < m_num_output_cols; i++) {

      const Vector4i bounds = m_disp_bound_image(i,j);
      //int height   = (bounds[3] - bounds[1] + 1);
      int width    = (bounds[2] - bounds[0] + 1);

      // Check the input image to find masked pixels
      PixelMask<Vector2i> integer_pixel = integer_disparity(i, j);
      if (!is_valid(integer_pixel)) {
        disparity(i,j) = integer_pixel;
        invalidate(disparity(i,j));
        continue; // No need for subpixel here
      }

      int dx = integer_pixel[0]; // The input integer disparity
      int dy = integer_pixel[1];

      // Stop here if not doing subpixel processing
      if (m_subpixel_type == SUBPIXEL_NONE) {
        disparity(i,j) = p_type(dx, dy);
        continue;
      }

      // Linear index of the min offset that will be checked
      int min_index = (dy-bounds[1])*width + (dx-bounds[0]);

      // Get the indices of where to find the four adjacent pixels
      int x_left  = -1;
      int x_right =  1;
      int y_up    = -width;
      int y_down  =  width;

      // Don't interpolate out of bounds
      bool top_bound, bottom_bound, left_bound, right_bound;
      top_bound = bottom_bound = left_bound = right_bound = false;
      if (dx == bounds[0]) { x_left  = 0; left_bound   = true; }
      if (dx == bounds[2]) { x_right = 0; right_bound  = true; }
      if (dy == bounds[1]) { y_up    = 0; top_bound    = true; }
      if (dy == bounds[3]) { y_down  = 0; bottom_bound = true; }

      // Apply subpixel correction and apply
      AccumCostType const* accum_vec = get_accum_vector(i, j);

      bool debug = false;//((j==2937) && (i >= 4635) && (i <= 4645));

      bool valid = true;
      if (m_subpixel_type == SUBPIXEL_PARABOLA) {
        valid = fitter.find_peak(accum_vec[min_index+x_left+y_up],
                                 accum_vec[min_index+y_up],
                                 accum_vec[min_index+x_right+y_up],
                                 accum_vec[min_index+x_left],
                                 accum_vec[min_index],
                                 accum_vec[min_index+x_right],
                                 accum_vec[min_index+x_left+y_down],
                                 accum_vec[min_index+y_down],
                                 accum_vec[min_index+x_right+y_down],
                                 delta_x, delta_y);
      } else {
        // This branch handles all 1D interpolation methods.
        // - These methods are always considered valid.
        delta_x = compute_subpixel_offset(accum_vec[min_index+x_left], accum_vec[min_index],
                                          accum_vec[min_index+x_right],
                                          left_bound, right_bound, debug);
        delta_y = compute_subpixel_offset(accum_vec[min_index+y_up], accum_vec[min_index],
                                          accum_vec[min_index+y_down],
                                          top_bound, bottom_bound, false);
      }

      // To assist development, write the internal subpixel input ratio to a file.
      //double temp = compute_subpixel_ratio(accum_vec[min_index+x_left], accum_vec[min_index], accum_vec[min_index+x_right]);
      //rawFile << temp << "\n";

      if (valid) {
        disparity(i,j) = p_type(dx+delta_x, dy+delta_y);
        //hist_dx.add_value(delta_x);
        //hist_dy.add_value(delta_y);
      } else {
        disparity(i,j) = p_type(dx, dy);
        percent_bad += 1.0;
      }
    } // End col loop
  } // End row loop

  if (m_subpixel_type == SUBPIXEL_PARABOLA) { // Only ever failures with this mode
    percent_bad /= (double)(m_num_output_rows*m_num_output_cols);
    vw_out(DebugMessage, "stereo") << "Subpixel interpolation failure percentage: " << percent_bad << "\n";
  }

  // Write these out for debugging/development
  //write_image("subpixel_disp.tif", disparity);
  //hist_dx.write_to_disk("delta_x.csv");
  //hist_dy.write_to_disk("delta_y.csv");
  //rawFile.close();

  return disparity;

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
    return static_cast<CostType>(abs(diff));
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
  } // End normalized cross correlaton case

  // Mean of abs differences
  int sum=0, diff=0;
  for (int j=-half_kernel_size; j<=half_kernel_size; ++j) {
    for (int i=-half_kernel_size; i<=half_kernel_size; i++) {
      diff = static_cast<int>(left_image (left_x +i, left_y +j)) -
             static_cast<int>(right_image(right_x+i, right_y+j));
      sum += abs(diff);
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
                                         c, r, c+dx,r+dy, false);
          m_cost_buffer[cost_index] = cost;
          ++cost_index;
        }
      } // End disparity loops
    } // End x loop
  }// End y loop

} // End function fill_costs_block

void SemiGlobalMatcher::fill_costs_census3x3(ImageView<uint8> const& left_image,
                                             ImageView<uint8> const& right_image) {
  const int half_kernel = (m_kernel_size - 1) / 2;
  const int padding     = 2*half_kernel;

  // Compute the census value for each pixel.
  // - ROI handling could be fancier but this is simple and works.
  // - The 0,0 pixels in the left and right images are assumed to be aligned.

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

  // Compute the census value for each pixel.
  // - ROI handling could be fancier but this is simple and works.
  // - The 0,0 pixels in the left and right images are assumed to be aligned.
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

  // Compute the census value for each pixel.
  // - ROI handling could be fancier but this is simple and works.
  // - The 0,0 pixels in the left and right images are assumed to be aligned.
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

  // Compute the census value for each pixel.
  // - ROI handling could be fancier but this is simple and works.
  // - The 0,0 pixels in the left and right images are assumed to be aligned.

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
  //Timer timer("\tSGM Cost Calculation");
  if ((m_cost_type == CENSUS_TRANSFORM) || (m_cost_type == TERNARY_CENSUS_TRANSFORM)) {
    switch(m_kernel_size) {
    case 3:  fill_costs_census3x3(left_image, right_image); break;
    case 5:  fill_costs_census5x5(left_image, right_image); break;
    case 7:  fill_costs_census7x7(left_image, right_image); break;
    case 9:  fill_costs_census9x9(left_image, right_image); break;
    default: vw_throw(NoImplErr()
                      << "Census transforms are only available in size 3, 5, 7, and 9.\n"
                      << "Other cost mode options do not have this restriction.");
    };
  } else {
    // Use the default mean of diff cost function.
    // Throw an error here, as this approach produces bad results.
    vw::vw_throw(NoImplErr() << "With SGM/MGM, only the census transform cost mode "
                 << "gives good results.\n");
    fill_costs_block(left_image, right_image);
  }

  /*
  // DEBUG: Record the min and max cost for each pixel
  ImageView<CostType> min_cost_image(m_num_output_cols, m_num_output_rows),
                      max_cost_image(m_num_output_cols, m_num_output_rows),
                      mean_cost_image(m_num_output_cols, m_num_output_rows);
  size_t cost_index = 0;
  for (int r = m_min_row; r <= m_max_row; r++) { // For each row in left
    int output_row = r - m_min_row;
    for (int c = m_min_col; c <= m_max_col; c++) { // For each column in left
      int output_col = c - m_min_col;

      Vector4i pixel_disp_bounds = m_disp_bound_image(output_col, output_row);

      CostType min_cost = 255; // Max value of uint8
      CostType max_cost = 0;
      double mean_cost = 0;
      double num_disp =0;
      for (int dy = pixel_disp_bounds[1]; dy <= pixel_disp_bounds[3]; dy++) { // For each disparity
        for (int dx = pixel_disp_bounds[0]; dx <= pixel_disp_bounds[2]; dx++) {

          CostType value = m_cost_buffer[cost_index];
          if (value < min_cost)
            min_cost = value;
          if (value > max_cost)
            max_cost = value;
          mean_cost += static_cast<double>(value);
          num_disp += 1.0;

          ++cost_index;
        }
      } // End disparity loops
      min_cost_image(output_col, output_row) = min_cost;
      max_cost_image(output_col, output_row) = max_cost;
      mean_cost_image(output_col, output_row) = floor(mean_cost / num_disp);
    } // End x loop
  }// End y loop
  write_image("min_cost_image.tif", min_cost_image);
  write_image("max_cost_image.tif", max_cost_image);
  write_image("mean_cost_image.tif", mean_cost_image);
  */

/* DEBUG: Generate a scanline costs image
 int debug_row = m_num_output_rows * 0.5;
 if (m_num_output_rows > debug_row) {
   ImageView<uint8> cost_image(m_num_output_cols, m_num_disp); // TODO: Change type?
   std::fill(cost_image.data(), cost_image.data()+m_num_output_cols*m_num_disp, 255);
   for (int i = 0; i < m_num_output_cols; i++) {
     CostType* buff = get_cost_vector(i, debug_row);
     Vector4i bounds = m_disp_bound_image(i,debug_row);
     int index = 0;
     for (int dy = bounds[1]; dy < bounds[3]; dy++) {
       for (int dx = bounds[0]; dx < bounds[2]; dx++) {
         int d = xy_to_disp(dx, dy);
         cost_image(i,d) = buff[index]; // TODO: scale!
         ++index;
       }
     }
   }
   write_image("scanline_costs_block.tif",cost_image);
 }
*/

} // end compute_disparity_costs()

void SemiGlobalMatcher::two_trip_path_accumulation(ImageView<uint8> const& left_image) {

  //Timer timer_total("\tSGM Cost Propagation");

  /// Create an object to manage the temporary accumulation buffers that need to be used here.
  const int  num_paths_in_pass=4;
  const bool vertical=false;
  MultiAccumRowBuffer buf_mgr(this, num_paths_in_pass, vertical);

  // Init this buffer to bad scores representing disparities that were
  //  not in the search range for the given pixel.
  boost::shared_array<AccumCostType> full_prior_buffer;
  full_prior_buffer.reset(new AccumCostType[m_num_disp]);
  for (int i = 0; i < m_num_disp; i++)
    full_prior_buffer[i] = get_bad_accum_val();

  AccumCostType* full_prior_ptr = full_prior_buffer.get();
  AccumCostType* output_accum_ptr;
  const int last_column = m_num_output_cols - 1;
  const int last_row    = m_num_output_rows - 1;

  // Loop through all pixels in the output image for the first trip, top-left to bottom-right.
  for (int row=0; row<m_num_output_rows; ++row) {
    for (int col=0; col<m_num_output_cols; ++col) {

      //printf("Accum pass 1 col = %d, row = %d\n", col, row);

      int num_disp = get_num_disparities(col, row);
      CostType * const local_cost_ptr = get_cost_vector(col, row);
      bool debug = false;//(col==152) && (row == 12);

      // Top left
      output_accum_ptr = buf_mgr.get_output_accum_ptr(MultiAccumRowBuffer::TOP_LEFT);
      if ((row > 0) && (col > 0)) {
        // Fill in the accumulated value in the bottom buffer
        int pixel_diff = get_path_pixel_diff(left_image, col, row, 1, 1);
        AccumCostType* const prior_accum_ptr = buf_mgr.get_trailing_pixel_accum_ptr(-1, -1, MultiAccumRowBuffer::TOP_LEFT);
        evaluate_path(col, row, col-1, row-1,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr,
                       pixel_diff, debug);
      } else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

      // Top
      output_accum_ptr = buf_mgr.get_output_accum_ptr(MultiAccumRowBuffer::TOP);
      if (row > 0) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, 0, 1);
        AccumCostType* const prior_accum_ptr = buf_mgr.get_trailing_pixel_accum_ptr(0, -1, MultiAccumRowBuffer::TOP);
        evaluate_path(col, row, col, row-1,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr,
                       pixel_diff, debug);
      } else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

      // Top right
      output_accum_ptr = buf_mgr.get_output_accum_ptr(MultiAccumRowBuffer::TOP_RIGHT);
      if ((row > 0) && (col < last_column)) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, -1, 1);
        AccumCostType* const prior_accum_ptr = buf_mgr.get_trailing_pixel_accum_ptr(1, -1, MultiAccumRowBuffer::TOP_RIGHT);
        evaluate_path(col, row, col+1, row-1,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr,
                       pixel_diff, debug);
      } else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

      // Left
      output_accum_ptr = buf_mgr.get_output_accum_ptr(MultiAccumRowBuffer::LEFT);
      if (col > 0) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, 1, 0);
        AccumCostType* const prior_accum_ptr = buf_mgr.get_trailing_pixel_accum_ptr(-1, 0, MultiAccumRowBuffer::LEFT);
        evaluate_path(col, row, col-1, row,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr,
                       pixel_diff, debug);
      } else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

      buf_mgr.next_pixel();
    } // End col loop

    buf_mgr.next_row(row==m_num_output_rows-1);
  } // End row loop

  // Done with the first trip!
  buf_mgr.switch_trips();

  // Loop through all pixels in the output image for the first trip, bottom-right to top-left.
  for (int row = last_row; row >= 0; --row) {
    for (int col = last_column; col >= 0; --col) {

      int num_disp = get_num_disparities(col, row);
      CostType * const local_cost_ptr = get_cost_vector(col, row);
      bool debug = false;//((row == 244) && (col == 341));

      // Bottom right
      output_accum_ptr = buf_mgr.get_output_accum_ptr(MultiAccumRowBuffer::BOT_RIGHT);
      if ((row < last_row) && (col < last_column)) {
        // Fill in the accumulated value in the bottom buffer
        int pixel_diff = get_path_pixel_diff(left_image, col, row, -1, -1);
        AccumCostType* const prior_accum_ptr = buf_mgr.get_trailing_pixel_accum_ptr(1, 1, MultiAccumRowBuffer::BOT_RIGHT);
        evaluate_path(col, row, col+1, row+1,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr,
                       pixel_diff, debug);
      } else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

      // Bottom
      output_accum_ptr = buf_mgr.get_output_accum_ptr(MultiAccumRowBuffer::BOT);
      if (row < last_row) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, 0, -1);
        AccumCostType* const prior_accum_ptr = buf_mgr.get_trailing_pixel_accum_ptr(0, 1, MultiAccumRowBuffer::BOT);
        evaluate_path(col, row, col, row+1,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr,
                       pixel_diff, debug);
      } else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

      // Bottom left
      output_accum_ptr = buf_mgr.get_output_accum_ptr(MultiAccumRowBuffer::BOT_LEFT);
      if ((row < last_row) && (col > 0)) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, 1, -1);
        AccumCostType* const prior_accum_ptr = buf_mgr.get_trailing_pixel_accum_ptr(-1, 1, MultiAccumRowBuffer::BOT_LEFT);
        evaluate_path(col, row, col-1, row+1,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr,
                       pixel_diff, debug);
      } else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

      // Right
      output_accum_ptr = buf_mgr.get_output_accum_ptr(MultiAccumRowBuffer::RIGHT);
      if (col < last_column) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, -1, 0);
        AccumCostType* const prior_accum_ptr = buf_mgr.get_trailing_pixel_accum_ptr(1, 0, MultiAccumRowBuffer::RIGHT);
        evaluate_path(col, row, col+1, row,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr,
                       pixel_diff, debug);
      } else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

      buf_mgr.next_pixel();
    } // End col loop

    buf_mgr.next_row(row==0);
  } // End row loop

  // Done with both trips!
}



// This version of the function requires four passes and is based on the paper:
// MGM: A Significantly More Global Matching for Stereovision
void SemiGlobalMatcher::smooth_path_accumulation(ImageView<uint8> const& left_image) {

  //Timer timer_total("\tSGM Cost Propagation");

  const int PATHS_PER_PASS = 2;

  // Create objects to manage the temporary accumulation buffers that need to be used here.
  // - Two copies are needed here, one for the two horizontal passes and one for the two vertical passes
  MultiAccumRowBuffer buf_mgr_horiz(this, PATHS_PER_PASS, false);
  MultiAccumRowBuffer buf_mgr_vertical  (this, PATHS_PER_PASS, true);

  // Init this buffer to bad scores representing disparities that were
  //  not in the search range for the given pixel.
  boost::shared_array<AccumCostType> full_prior_buffer;
  full_prior_buffer.reset(new AccumCostType[m_num_disp]);
  for (int i = 0; i < m_num_disp; i++)
    full_prior_buffer[i] = get_bad_accum_val();

  AccumCostType* full_prior_ptr = full_prior_buffer.get();
  AccumCostType* output_accum_ptr;
  const int last_column = m_num_output_cols - 1;
  const int last_row    = m_num_output_rows - 1;

  // Loop through all pixels in the output image for the first trip, top-left to bottom-right.
  for (int row=0; row<m_num_output_rows; ++row) {
    for (int col=0; col<m_num_output_cols; ++col) {

      //printf("Accum pass 1 col = %d, row = %d\n", col, row);

      // TODO: Do we need to skip these pixels in the two-pass SGM function above?
      int num_disp = get_num_disparities(col, row);
      if (num_disp == 0) {
        buf_mgr_horiz.next_pixel();
        continue;
      }
      CostType * const local_cost_ptr = get_cost_vector(col, row);
      bool debug = false;//((row == 244) && (col == 341));

      boost::shared_array<AccumCostType> temp_buffer(new AccumCostType[num_disp]);

      // Left
      output_accum_ptr = buf_mgr_horiz.get_output_accum_ptr(MultiAccumRowBuffer::PASS_ONE);
      if ((row > 0) && (col > 0)) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, -1, 0);
        // Compute accumulation from the values in the left pixel
        AccumCostType* const prior_accum_ptr = buf_mgr_horiz.get_trailing_pixel_accum_ptr(-1, 0, MultiAccumRowBuffer::PASS_ONE);
        evaluate_path(col, row, col-1, row,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr,
                       pixel_diff, debug);

        // Compute accumulation from the values in the above pixel
        AccumCostType* const prior_accum_ptr2 = buf_mgr_horiz.get_trailing_pixel_accum_ptr(0, -1, MultiAccumRowBuffer::PASS_ONE);
        evaluate_path(col, row, col, row-1,
                       prior_accum_ptr2, full_prior_ptr, local_cost_ptr, temp_buffer.get(),
                       pixel_diff, debug);
        // The final accumulation values are the average of the two computations
        for (int d=0; d<num_disp; ++d)
          output_accum_ptr[d] = (output_accum_ptr[d] + temp_buffer[d])/2;
      } else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

      // Top left
      output_accum_ptr = buf_mgr_horiz.get_output_accum_ptr(MultiAccumRowBuffer::PASS_TWO);
      if ((row > 0) && (col > 0) && (col < last_column)) {

        int pixel_diff = get_path_pixel_diff(left_image, col, row, -1, -1);
        AccumCostType* const prior_accum_ptr
          = buf_mgr_horiz.get_trailing_pixel_accum_ptr(-1, -1, MultiAccumRowBuffer::PASS_TWO);
        evaluate_path(col, row, col-1, row-1,
                      prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr,
                      pixel_diff, debug);

        AccumCostType* const prior_accum_ptr2
        = buf_mgr_horiz.get_trailing_pixel_accum_ptr(1, -1, MultiAccumRowBuffer::PASS_TWO);
        evaluate_path(col, row, col+1, row-1,
                      prior_accum_ptr2, full_prior_ptr, local_cost_ptr, temp_buffer.get(),
                      pixel_diff, debug);
        for (int d=0; d<num_disp; ++d)
          output_accum_ptr[d] = (output_accum_ptr[d] + temp_buffer[d])/2;
      } else {
        // Just init to the local cost
        for (int d=0; d<num_disp; ++d)
          output_accum_ptr[d] = local_cost_ptr[d];
      }

      buf_mgr_horiz.next_pixel();
    } // End col loop

    buf_mgr_horiz.next_row(row==last_row);
  } // End row loop

  // Done with the first trip!
  buf_mgr_horiz.switch_trips();

  // Loop through all pixels in the output image for the second trip,
  // bottom-right to top-left.
  for (int row = last_row; row >= 0; --row) {
    for (int col = last_column; col >= 0; --col) {

      int num_disp = get_num_disparities(col, row);
      if (num_disp == 0) {
        buf_mgr_horiz.next_pixel();
        continue;
      }
      CostType * const local_cost_ptr = get_cost_vector(col, row);
      bool debug = false;//((row == 244) && (col == 341));

      boost::shared_array<AccumCostType> temp_buffer(new AccumCostType[num_disp]);

      // Right
      output_accum_ptr = buf_mgr_horiz.get_output_accum_ptr(MultiAccumRowBuffer::PASS_ONE);
      if ((row < last_row) && (col < last_column)) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, 1, 0);
        AccumCostType* const prior_accum_ptr = buf_mgr_horiz.get_trailing_pixel_accum_ptr(1, 0, MultiAccumRowBuffer::PASS_ONE);
        evaluate_path(col, row, col+1, row,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr,
                       pixel_diff, debug);

        AccumCostType* const prior_accum_ptr2 = buf_mgr_horiz.get_trailing_pixel_accum_ptr(0, 1, MultiAccumRowBuffer::PASS_ONE);
        evaluate_path(col, row, col, row+1,
                       prior_accum_ptr2, full_prior_ptr, local_cost_ptr, temp_buffer.get(),
                       pixel_diff, debug);
        for (int d=0; d<num_disp; ++d)
          output_accum_ptr[d] = (output_accum_ptr[d] + temp_buffer[d])/2;
      } else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];


      // Bottom right
      output_accum_ptr = buf_mgr_horiz.get_output_accum_ptr(MultiAccumRowBuffer::PASS_TWO);
      if ((row < last_row) && (col > 0) && (col < last_column)) {

        int pixel_diff = get_path_pixel_diff(left_image, col, row, 1, 1);
        AccumCostType* const prior_accum_ptr = buf_mgr_horiz.get_trailing_pixel_accum_ptr(1, 1, MultiAccumRowBuffer::PASS_TWO);
        evaluate_path(col, row, col+1, row+1,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr,
                       pixel_diff, debug);

        AccumCostType* const prior_accum_ptr2 = buf_mgr_horiz.get_trailing_pixel_accum_ptr(-1, 1, MultiAccumRowBuffer::PASS_TWO);
        evaluate_path(col, row, col-1, row+1,
                       prior_accum_ptr2, full_prior_ptr, local_cost_ptr, temp_buffer.get(),
                       pixel_diff, debug);
        for (int d=0; d<num_disp; ++d)
          output_accum_ptr[d] = (output_accum_ptr[d] + temp_buffer[d])/2;
      } else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];


      buf_mgr_horiz.next_pixel();
    } // End col loop

    buf_mgr_horiz.next_row(row==0);
  } // End row loop

  // -- Done with the horizontal passes, switch to vertical passes.

  // Loop through all pixels in the output image for the third trip, bottom-left to top-right.
  for (int col = 0; col < m_num_output_cols; ++col) {
    for (int row = last_row; row >= 0; --row) {

      int num_disp = get_num_disparities(col, row);
      //printf("Accum pass 3 col = %d, row = %d, num_disp = %d\n", col, row, num_disp);
      if (num_disp == 0) {
        buf_mgr_vertical.next_pixel();
        continue;
      }
      CostType * const local_cost_ptr = get_cost_vector(col, row);
      bool debug = false;//((row == 244) && (col == 341));

      boost::shared_array<AccumCostType> temp_buffer(new AccumCostType[num_disp]);

      // Bottom
      output_accum_ptr = buf_mgr_vertical.get_output_accum_ptr(MultiAccumRowBuffer::PASS_ONE);
      if ((row < last_row) && (col > 0)) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, 0, 1);
        AccumCostType* const prior_accum_ptr = buf_mgr_vertical.get_trailing_pixel_accum_ptr(0, 1, MultiAccumRowBuffer::PASS_ONE);
        evaluate_path(col, row, col, row+1,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr,
                       pixel_diff, debug);

        AccumCostType* const prior_accum_ptr2 = buf_mgr_vertical.get_trailing_pixel_accum_ptr(-1, 0, MultiAccumRowBuffer::PASS_ONE);
        evaluate_path(col, row, col-1, row,
                       prior_accum_ptr2, full_prior_ptr, local_cost_ptr, temp_buffer.get(),
                       pixel_diff, debug);
        for (int d=0; d<num_disp; ++d)
          output_accum_ptr[d] = (output_accum_ptr[d] + temp_buffer[d])/2;
      } else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

      // Bottom left
      output_accum_ptr = buf_mgr_vertical.get_output_accum_ptr(MultiAccumRowBuffer::PASS_TWO);
      if ((row > 0) && (row < last_row) && (col > 0)) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, -1, 1);
        AccumCostType* const prior_accum_ptr = buf_mgr_vertical.get_trailing_pixel_accum_ptr(-1, 1, MultiAccumRowBuffer::PASS_TWO);
        evaluate_path(col, row, col-1, row+1,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr,
                       pixel_diff, debug);

        AccumCostType* const prior_accum_ptr2 = buf_mgr_vertical.get_trailing_pixel_accum_ptr(-1, -1, MultiAccumRowBuffer::PASS_TWO);
        evaluate_path(col, row, col-1, row-1,
                       prior_accum_ptr2, full_prior_ptr, local_cost_ptr, temp_buffer.get(),
                       pixel_diff, debug);
        for (int d=0; d<num_disp; ++d)
          output_accum_ptr[d] = (output_accum_ptr[d] + temp_buffer[d])/2;
      } else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

      buf_mgr_vertical.next_pixel();
    } // End col loop

    buf_mgr_vertical.next_row(col==last_column);
  } // End row loop

  // Done with the first trip!
  buf_mgr_vertical.switch_trips();

  // Loop through all pixels in the output image for the fourth trip, top-right to bottom-left.
  for (int col=last_column; col>=0; --col) {
    for (int row=0; row<m_num_output_rows; ++row) {

      //printf("Accum pass 1 col = %d, row = %d\n", col, row);

      int num_disp = get_num_disparities(col, row);
      if (num_disp == 0) {
        buf_mgr_vertical.next_pixel();
        continue;
      }
      CostType * const local_cost_ptr = get_cost_vector(col, row);
      bool debug = false;//((row == 244) && (col == 341));

      boost::shared_array<AccumCostType> temp_buffer(new AccumCostType[num_disp]);

      // Top
      output_accum_ptr = buf_mgr_vertical.get_output_accum_ptr(MultiAccumRowBuffer::PASS_ONE);
      if ((row > 0) && (col < last_column)) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, 0, -1);
        AccumCostType* const prior_accum_ptr = buf_mgr_vertical.get_trailing_pixel_accum_ptr(0, -1, MultiAccumRowBuffer::PASS_ONE);
        evaluate_path(col, row, col, row-1,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr,
                       pixel_diff, debug);

        AccumCostType* const prior_accum_ptr2 = buf_mgr_vertical.get_trailing_pixel_accum_ptr(1, 0, MultiAccumRowBuffer::PASS_ONE);
        evaluate_path(col, row, col+1, row,
                       prior_accum_ptr2, full_prior_ptr, local_cost_ptr, temp_buffer.get(),
                       pixel_diff, debug);
        for (int d=0; d<num_disp; ++d)
          output_accum_ptr[d] = (output_accum_ptr[d] + temp_buffer[d])/2;
      } else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

      // Top right
      output_accum_ptr = buf_mgr_vertical.get_output_accum_ptr(MultiAccumRowBuffer::PASS_TWO);
      if ((row > 0) && (row < last_row) && (col < last_column)) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, 1, -1);
        AccumCostType* const prior_accum_ptr = buf_mgr_vertical.get_trailing_pixel_accum_ptr(1, -1, MultiAccumRowBuffer::PASS_TWO);
        evaluate_path(col, row, col+1, row-1,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr,
                       pixel_diff, debug);

        AccumCostType* const prior_accum_ptr2 = buf_mgr_vertical.get_trailing_pixel_accum_ptr(1, 1, MultiAccumRowBuffer::PASS_TWO);
        evaluate_path(col, row, col+1, row+1,
                       prior_accum_ptr2, full_prior_ptr, local_cost_ptr, temp_buffer.get(),
                       pixel_diff, debug);
        for (int d=0; d<num_disp; ++d)
          output_accum_ptr[d] = (output_accum_ptr[d] + temp_buffer[d])/2;
      } else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

      buf_mgr_vertical.next_pixel();
    } // End col loop

    buf_mgr_vertical.next_row(col==0);
  } // End row loop

  // Done with all trips!
} // End function smooth_path_accumulation

SemiGlobalMatcher::DisparityImage
SemiGlobalMatcher::semi_global_matching_func(ImageView<uint8> const& left_image,
                                             ImageView<uint8> const& right_image,
                                             ImageView<uint8> const* left_image_mask,
                                             ImageView<uint8> const* right_image_mask,
                                             DisparityImage const* prev_disparity) {

  // Compute safe bounds to search through given the disparity range and kernel size.
  // - Using inclusive bounds here.

  const int half_kernel_size = (m_kernel_size-1) / 2;

  // Figure out the maximum possible image extent that we can compute stereo for
  //  given the size of two input images and the requirement that we fully search
  //  the specified search range.
  // - The search region = the region in the left image + the search bounds, shrunk by
  //   half the cost function kernel size.
  m_min_row = half_kernel_size - m_min_disp_y; // Assumes the (0,0) pixels are aligned
  m_min_col = half_kernel_size - m_min_disp_x;
  int left_last_col  = left_image.cols () - 1;
  int left_last_row  = left_image.rows () - 1;
  int right_last_col = right_image.cols() - 1;
  int right_last_row = right_image.rows() - 1;
  m_max_row = std::min(left_last_row  -  half_kernel_size,
                       right_last_row - (half_kernel_size + m_max_disp_y));
  m_max_col = std::min(left_last_col  - half_kernel_size,
                       right_last_col - (half_kernel_size + m_max_disp_x));
  if (m_min_row < 0) m_min_row = 0;
  if (m_min_col < 0) m_min_col = 0;
  if (m_max_row > left_last_row) m_max_row = left_last_row;
  if (m_max_col > left_last_col) m_max_col = left_last_col;

  m_num_output_cols  = m_max_col - m_min_col + 1;
  m_num_output_rows  = m_max_row - m_min_row + 1;

  populate_adjacent_disp_lookup_table();

  // By default the search bounds are the same for each pixel,
  //  but set them from the prior disparity image if the user passed it in.
  populate_constant_disp_bound_image();

  if (!populate_disp_bound_image(left_image_mask, right_image_mask, prev_disparity)) {
    vw_out(WarningMessage, "stereo")
      << "Unable to compute valid search ranges for SGM input. Increase --corr-memory-mb.\n";
    // If the inputs are invalid, return a default disparity image.
    DisparityImage disparity(m_num_output_cols, m_num_output_rows);
    return invalidate_mask(disparity);
  }

  // All the hard work is done in the next few function calls
  allocate_large_buffers();
  compute_disparity_costs(left_image, right_image);
  if (m_use_mgm)
    accum_mgm_multithread(left_image);
  else
    accum_sgm_multithread(left_image);

  // Now that all the costs are calculated, fetch the best disparity for each
  // pixel. This computes integer disparities. Subpixel disparities are computed
  // in CorrelationView.tcc.
  return create_disparity_view();
}

// If some success values are 0, throw an exception that will be caught by the caller.
void SgmThrowOnFailure(std::vector<int> const& success) {
  for (int i = 0; i < (int)success.size(); i++) {
    if (success[i] == 0)
      vw_throw(vw::ArgumentErr() << "Failed to run correlation.\n");
  }
}

// Perform standard SGM path accumulation using N threads.
// Individual threads use a flag to indicate failure, as exceptions cannot be
// thrown in a thread. Then, this function throws an exception if any thread
// failed.
void SemiGlobalMatcher::accum_sgm_multithread(ImageView<uint8> const& left_image) {

  int num_threads = vw_settings().default_num_threads();
  int height = m_num_output_rows;
  int width  = m_num_output_cols;

  ImageView<uint8> const* image_ptr = &left_image;
  Vector2i image_size(width, height);

  // Start up the thread pool
  FifoWorkQueue thread_pool(num_threads);

  // Initialize a number of line buffers equal to the number of threads
  OneLineBufferManager mem_buf_mgr(num_threads, this);

  // Each direction of lines is handled separately.  This adds some additional wait time
  //  as each direction finishes up but it allows the threads to simultaneously read/write
  //  the main accumulation buffer without any collisions.

  //  Must track the success of each thread. That because exceptions cannot be
  // thrown in a thread. A value of 1 in this vector means success.
  std::vector<int> success;
  int thread_id = 0;
  int max_len = width + height; // To not use custom sizes for the vectors below
  typedef boost::shared_ptr<PixelPassTask> TaskPtrType;

  // Add lines going down
  success.resize(max_len, 1); thread_id = 0;
  for (int i = 0; i < width; i++) {
    Vector2i top_pixel(i, 0);
    PixelLineIterator line_from_top(top_pixel, PixelLineIterator::B, image_size);
    TaskPtrType task(new PixelPassTask(image_ptr, this, &mem_buf_mgr, line_from_top,
                                       &success[thread_id]));
    thread_pool.add_task(task); thread_id++;
  }
  thread_pool.join_all(); // Wait for all tasks to complete
  SgmThrowOnFailure(success);

  // Add lines going up
  success.resize(max_len, 1); thread_id = 0;
  for (int i = 0; i < width; i++) {
    Vector2i bottom_pixel(i, height-1);
    PixelLineIterator line_from_bottom(bottom_pixel, PixelLineIterator::T, image_size);
    TaskPtrType task(new PixelPassTask(image_ptr, this, &mem_buf_mgr, line_from_bottom,
                                       &success[thread_id]));
    thread_pool.add_task(task); thread_id++;
  }
  thread_pool.join_all(); // Wait for all tasks to complete
  SgmThrowOnFailure(success);

  // Add lines from the left
  success.resize(max_len, 1); thread_id = 0;
  for (int i = 0; i < height; i++) {
    Vector2i left_pixel(0, i);
    PixelLineIterator line_from_left (left_pixel, PixelLineIterator::R, image_size);
    TaskPtrType task(new PixelPassTask(image_ptr, this, &mem_buf_mgr, line_from_left,
                                       &success[thread_id]));
    thread_pool.add_task(task); thread_id++;
  }
  thread_pool.join_all(); // Wait for all tasks to complete
  SgmThrowOnFailure(success);

  // Add lines from the right
  success.resize(max_len, 1); thread_id = 0;
  for (int i = 0; i < height; i++) {
    Vector2i right_pixel(width-1, i);
    PixelLineIterator line_from_right(right_pixel, PixelLineIterator::L, image_size);
    TaskPtrType task(new PixelPassTask(image_ptr, this, &mem_buf_mgr, line_from_right,
                                       &success[thread_id]));
    thread_pool.add_task(task); thread_id++;
  }
  thread_pool.join_all(); // Wait for all tasks to complete
  SgmThrowOnFailure(success);

  // Add lines from the top left
  success.resize(max_len, 1); thread_id = 0;
  for (int i = 0; i < width; i++) {
    Vector2i top_pixel(i, 0);
    PixelLineIterator line_from_top_br(top_pixel, PixelLineIterator::BR, image_size);
    TaskPtrType task(new PixelPassTask(image_ptr, this, &mem_buf_mgr, line_from_top_br,
                                       &success[thread_id]));
    thread_pool.add_task(task); thread_id++;
  }
  for (int i = 1; i < height; i++) {
    Vector2i left_pixel(0, i);
    PixelLineIterator line_from_left_br(left_pixel, PixelLineIterator::BR, image_size);
    TaskPtrType task(new PixelPassTask(image_ptr, this, &mem_buf_mgr, line_from_left_br,
                                       &success[thread_id]));
    thread_pool.add_task(task); thread_id++;
  }
  thread_pool.join_all(); // Wait for all tasks to complete
  SgmThrowOnFailure(success);

  // Add lines from the top right
  success.resize(max_len, 1); thread_id = 0;
  for (int i = 0; i < width; i++) {
    Vector2i top_pixel(i, 0);
    PixelLineIterator line_from_top_bl(top_pixel, PixelLineIterator::BL, image_size);
    TaskPtrType task(new PixelPassTask(image_ptr, this, &mem_buf_mgr, line_from_top_bl,
                                       &success[thread_id]));
    thread_pool.add_task(task); thread_id++;
  }
  for (int i = 1; i < height; i++) {
    Vector2i right_pixel(width-1, i);
    PixelLineIterator line_from_right_bl(right_pixel, PixelLineIterator::BL, image_size);
    TaskPtrType task(new PixelPassTask(image_ptr, this, &mem_buf_mgr, line_from_right_bl,
                                       &success[thread_id]));
    thread_pool.add_task(task); thread_id++;
  }
  thread_pool.join_all(); // Wait for all tasks to complete
  SgmThrowOnFailure(success);

  // Add lines from the bottom left
  success.resize(max_len, 1); thread_id = 0;
  for (int i = 0; i < width; i++) {
    Vector2i bot_pixel(i, height-1);
    PixelLineIterator line_from_bot_tr(bot_pixel, PixelLineIterator::TR, image_size);
    TaskPtrType task(new PixelPassTask(image_ptr, this, &mem_buf_mgr, line_from_bot_tr,
                                       &success[thread_id]));
    thread_pool.add_task(task); thread_id++;
  }
  for (int i = 0; i < height-1; i++) {
    Vector2i left_pixel(0, i);
    PixelLineIterator line_from_left_tr(left_pixel, PixelLineIterator::TR, image_size);
    TaskPtrType task(new PixelPassTask(image_ptr, this, &mem_buf_mgr, line_from_left_tr,
                                       &success[thread_id]));
    thread_pool.add_task(task); thread_id++;
  }
  thread_pool.join_all(); // Wait for all tasks to complete
  SgmThrowOnFailure(success);

  // Add lines from the bottom right
  success.resize(max_len, 1); thread_id = 0;
  for (int i = 0; i < width; i++) {
    Vector2i bot_pixel(i, height-1);
    PixelLineIterator line_from_bot_tl(bot_pixel, PixelLineIterator::TL, image_size);
    TaskPtrType task(new PixelPassTask(image_ptr, this, &mem_buf_mgr, line_from_bot_tl,
                                       &success[thread_id]));
    thread_pool.add_task(task); thread_id++;
  }
  for (int i = 0; i < height-1; i++) {
    Vector2i right_pixel(width-1, i);
    PixelLineIterator line_from_right_tl(right_pixel, PixelLineIterator::TL, image_size);
    TaskPtrType task(new PixelPassTask(image_ptr, this, &mem_buf_mgr, line_from_right_tl,
                                       &success[thread_id]));
    thread_pool.add_task(task); thread_id++;
  }
  thread_pool.join_all(); // Wait for all tasks to complete
  SgmThrowOnFailure(success);

} // End function accum_sgm_multithread

// This version of the function requires four passes and is based on the paper:
// MGM: A Significantly More Global Matching for Stereovision.
// Individual threads use a flag to indicate failure, as exceptions cannot be
// thrown in a thread. Then, this function throws an exception if any thread
// failed.
void SemiGlobalMatcher::accum_mgm_multithread(ImageView<uint8> const& left_image) {

  const int PATHS_PER_PASS = 1;
  const int MAX_USABLE_THREADS = 8;

  int num_threads = vw_settings().default_num_threads();

  // In order to use more than 8 threads at a time we need a more complicated MGM
  // implementation.  Currently each of the directions gets its own thread.
  if (num_threads > MAX_USABLE_THREADS) {
    num_threads = MAX_USABLE_THREADS;
    vw_out(WarningMessage) << "MGM stereo processing is currently limited "
      << "to 8 parallel threads." << "\n";
  }

  // Create objects to manage the temporary accumulation buffers that need to be used here.
  // - A separate copy instance is used for each direction to allow multithreading
  MultiAccumRowBuffer buf_mgr_horiz_left(this, PATHS_PER_PASS, false);
  MultiAccumRowBuffer buf_mgr_horiz_right(this, PATHS_PER_PASS, false);
  MultiAccumRowBuffer buf_mgr_horiz_topleft(this, PATHS_PER_PASS, false);
  MultiAccumRowBuffer buf_mgr_horiz_botright(this, PATHS_PER_PASS, false);

  MultiAccumRowBuffer buf_mgr_vertical_top(this, PATHS_PER_PASS, true);
  MultiAccumRowBuffer buf_mgr_vertical_bot(this, PATHS_PER_PASS, true);
  MultiAccumRowBuffer buf_mgr_vertical_topright(this, PATHS_PER_PASS, true);
  MultiAccumRowBuffer buf_mgr_vertical_botleft(this, PATHS_PER_PASS, true);

  // Set some of the buffers to the reverse direction
  buf_mgr_horiz_right.switch_trips();
  buf_mgr_horiz_botright.switch_trips();
  buf_mgr_vertical_top.switch_trips();
  buf_mgr_vertical_topright.switch_trips();

  // Start up thread pool
  FifoWorkQueue thread_pool(num_threads);

  //  Must track the success of each thread. That because exceptions
  // cannot be thrown in a thread.
  std::vector<int> success(8, 1);

  // Load all eight required passes as tasks in the thread pool
  typedef boost::shared_ptr<SmoothPathAccumTask> TaskPtrType;
  TaskPtrType task_L (new SmoothPathAccumTask
                      (&buf_mgr_horiz_left, this, &left_image,
                       SmoothPathAccumTask::L, &success[0]));
  TaskPtrType task_R (new SmoothPathAccumTask
                      (&buf_mgr_horiz_right, this, &left_image,
                       SmoothPathAccumTask::R, &success[1]));
  TaskPtrType task_TL(new SmoothPathAccumTask
                      (&buf_mgr_horiz_topleft, this, &left_image,
                       SmoothPathAccumTask::TL, &success[2]));
  TaskPtrType task_BR(new SmoothPathAccumTask
                      (&buf_mgr_horiz_botright, this, &left_image,
                       SmoothPathAccumTask::BR, &success[3]));
  TaskPtrType task_T (new SmoothPathAccumTask
                      (&buf_mgr_vertical_top, this, &left_image,
                       SmoothPathAccumTask::T, &success[4]));
  TaskPtrType task_B (new SmoothPathAccumTask
                      (&buf_mgr_vertical_bot, this, &left_image,
                       SmoothPathAccumTask::B, &success[5]));
  TaskPtrType task_TR(new SmoothPathAccumTask
                      (&buf_mgr_vertical_topright, this, &left_image,
                       SmoothPathAccumTask::TR, &success[6]));
  TaskPtrType task_BL(new SmoothPathAccumTask
                      (&buf_mgr_vertical_botleft, this, &left_image,
                       SmoothPathAccumTask::BL, &success[7]));

  thread_pool.add_task(task_L);
  thread_pool.add_task(task_R);
  thread_pool.add_task(task_TL);
  thread_pool.add_task(task_BR);
  thread_pool.add_task(task_T);
  thread_pool.add_task(task_B);
  thread_pool.add_task(task_TR);
  thread_pool.add_task(task_BL);

  thread_pool.join_all(); // Wait for all tasks to complete

  // If some success values are 0, throw an exception that will be caught by the caller.
  SgmThrowOnFailure(success);

} // End function accum_mgm_multithread

/// Get a pointer to a cost vector
SemiGlobalMatcher::CostType * SemiGlobalMatcher::get_cost_vector(int col, int row) {
  size_t start_index = m_buffer_starts(col, row);
  return m_cost_buffer.get() + start_index;
}

/// Get a pointer to an accumulated cost vector
SemiGlobalMatcher::AccumCostType* SemiGlobalMatcher::get_accum_vector(int col, int row) {
  size_t start_index = m_buffer_starts(col, row);
  return m_accum_buffer.get() + start_index;
}

/// Get the pixel diff along a line at a specified output location.
int SemiGlobalMatcher::get_path_pixel_diff(ImageView<uint8> const& left_image,
                        int col, int row, int dir_x, int dir_y) const {
  // Take the offset between the output location and the input pixel coordinates.
  int a = left_image(col        +m_min_col, row        +m_min_row);
  int b = left_image((col-dir_x)+m_min_col, (row-dir_y)+m_min_row);
  return std::abs(a - b);
}

SemiGlobalMatcher::DisparityType
SemiGlobalMatcher::xy_to_disp(DisparityType dx, DisparityType dy) const {
    return (dy-m_min_disp_y)*m_num_disp_x + (dx-m_min_disp_x);
}

/// Converts from a linear disparity index to the dx, dy values it represents.
/// - This function is too slow to use inside the inner loop!
void SemiGlobalMatcher::disp_to_xy(DisparityType disp, DisparityType &dx,
                                   DisparityType &dy) const {
  dy = (disp / m_num_disp_x) + m_min_disp_y; // 2D implementation
  dx = (disp % m_num_disp_x) + m_min_disp_x;
}

/// Convert a pixel's minimum disparity index to dx, dy.
void SemiGlobalMatcher::disp_index_to_xy(int min_index, int col, int row,
                        DisparityType &dx, DisparityType &dy) const {
  // Convert the disparity index to dx and dy
  const Vector4i bounds = m_disp_bound_image(col,row);
  int d_width  = bounds[2] - bounds[0] + 1;
  dy = (min_index / d_width);
  dx = min_index - (dy*d_width) + bounds[0];
  dy += bounds[1];
}

} // end namespace stereo
} // end namespace vw

