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

/// \file SGMAccum.cc
/// Path cost accumulation routines for Semi-Global Matching

#include <queue>
#include <math.h>
#include <vw/Stereo/SGM.h>
#include <vw/Stereo/SGMAssist.h>
#include <vw/Core/Debugging.h>
#include <vw/Core/ThreadPool.h>

#if defined(VW_ENABLE_SSE) && (VW_ENABLE_SSE==1)
  #include <emmintrin.h>
  #include <smmintrin.h> // SSE4.1
#endif

namespace vw {

Mutex MultiAccumRowBuffer::m_mutex;

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

/// Get the pixel diff along a line at a specified output location.
int SemiGlobalMatcher::get_path_pixel_diff(ImageView<uint8> const& left_image,
                        int col, int row, int dir_x, int dir_y) const {
  // Take the offset between the output location and the input pixel coordinates.
  int a = left_image(col        +m_min_col, row        +m_min_row);
  int b = left_image((col-dir_x)+m_min_col, (row-dir_y)+m_min_row);
  return std::abs(a - b);
}

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


} // end namespace vw
