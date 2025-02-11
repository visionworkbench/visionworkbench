
#ifndef __SGM_ASSIST_H__
#define __SGM_ASSIST_H__

#include <vw/Stereo/SGM.h>
#include <vw/Core/Thread.h>
#include <vw/Core/ThreadPool.h>
#include <vw/Image/PixelIterator.h>
#include <vw/Core/Log.h>

/**
  This file contains supporting classes and functions for the SGM algorithm.
*/


namespace vw {

namespace stereo {

/// Class to compute parabola surface sub-pixel fits
class ParabolaFit2d {

public: // Functions

  /// Constructor
  ParabolaFit2d() {

    // We get a considerable speedup in our 2d subpixel correlation if
    // we go ahead and compute the pseudoinverse of the A matrix (where
    // each row in A is [ x^2 y^2 xy x y 1] (our 2d parabolic surface)
    // for the range of x = [-1:1] and y = [-1:1].
/* From ASP    
    static double pinvA_data[] =
      { 1.0/6,  1.0/6,  1.0/6, -1.0/3, -1.0/3, -1.0/3,  1.0/6,  1.0/6,  1.0/6,
        1.0/6, -1.0/3,  1.0/6,  1.0/6, -1.0/3,  1.0/6,  1.0/6, -1.0/3,  1.0/6,
        1.0/4,    0.0, -1.0/4,    0.0,    0.0,    0.0, -1.0/4,    0.0,  1.0/4,
        -1.0/6, -1.0/6, -1.0/6,    0.0,    0.0,   0.0,  1.0/6,  1.0/6,  1.0/6,
        -1.0/6,    0.0,  1.0/6, -1.0/6,    0.0, 1.0/6, -1.0/6,    0.0,  1.0/6,
        -1.0/9,  2.0/9, -1.0/9,  2.0/9,  5.0/9, 2.0/9, -1.0/9,  2.0/9, -1.0/9 };*/

/* From paper
    static double pinvA_data[] =
      {  1.0/6, -1.0/3,  1.0/6,  1.0/6, -1.0/3,  1.0/6,   1.0/6, -1.0/3,  1.0/6,  // = a
         1.0/6,  1.0/6,  1.0/6, -1.0/3, -1.0/3, -1.0/3,   1.0/6,  1.0/6,  1.0/6,  // = b
        -1.0/4,    0.0,  1.0/4,    0.0,    0.0,    0.0,   1.0/4,    0.0, -1.0/4,  // = c
        -1.0/6,    0.0,  1.0/6, -1.0/6,    0.0,  1.0/6,  -1.0/6,    0.0,  1.0/6,  // = d
         1.0/6,  1.0/6,  1.0/6,    0.0,    0.0,    0.0,  -1.0/6, -1.0/6, -1.0/6,  // = e
        -1.0/9,  2.0/9, -1.0/9,  2.0/9,   5.0/9, 2.0/9,  -1.0/9,  2.0/9, -1.0/9 };// = f*/

/* From MATLAB:
"
m = [1 1  1 -1 -1 1;   
     0 1  0  0 -1 1;  
     1 1 -1  1 -1 1;   
     1 0  0 -1  0 1;   
     0 0  0  0  0 1; 
     1 0  0  1  0 1;
     1 1 -1 -1  1 1;  
     0 1  0  0  1 1;   
     1 1  1  1  1 1];
 
 pinv(m)
"
    0.1667   -0.3333    0.1667    0.1667   -0.3333    0.1667    0.1667   -0.3333    0.1667
    0.1667    0.1667    0.1667   -0.3333   -0.3333   -0.3333    0.1667    0.1667    0.1667
    0.2500    0.0000   -0.2500   -0.0000   -0.0000    0.0000   -0.2500   -0.0000    0.2500
   -0.1667    0.0000    0.1667   -0.1667   -0.0000    0.1667   -0.1667   -0.0000    0.1667
   -0.1667   -0.1667   -0.1667    0.0000    0.0000   -0.0000    0.1667    0.1667    0.1667
   -0.1111    0.2222   -0.1111    0.2222    0.5556    0.2222   -0.1111    0.2222   -0.1111
*/        
    static double pinvA_data[] =
      {  1.0/6, -1.0/3,  1.0/6,  1.0/6, -1.0/3,  1.0/6,   1.0/6, -1.0/3,  1.0/6,  // = a
         1.0/6,  1.0/6,  1.0/6, -1.0/3, -1.0/3, -1.0/3,   1.0/6,  1.0/6,  1.0/6,  // = b
         1.0/4,    0.0, -1.0/4,    0.0,    0.0,    0.0,  -1.0/4,    0.0,  1.0/4,  // = c
        -1.0/6,    0.0,  1.0/6, -1.0/6,    0.0,  1.0/6,  -1.0/6,    0.0,  1.0/6,  // = d
        -1.0/6, -1.0/6, -1.0/6,    0.0,    0.0,    0.0,   1.0/6,  1.0/6,  1.0/6,  // = e
        -1.0/9,  2.0/9, -1.0/9,  2.0/9,   5.0/9, 2.0/9,  -1.0/9,  2.0/9, -1.0/9 };// = f

    m_fit_params = Matrix<double,6,9>(pinvA_data);
  }
  
  /// Find the peak of the surface fit nearby central point z5.
  /// - dx and dy are relative to z5 at coordinate (0,0)
  bool find_peak(double z1, double z2, double z3,
                 double z4, double z5, double z6,
                 double z7, double z8, double z9,
                 double &dx, double &dy) {
    Vector<double,9> z_vec;
    z_vec[0] = z1;  z_vec[1] = z2;  z_vec[2] = z3;
    z_vec[3] = z4;  z_vec[4] = z5;  z_vec[5] = z6;
    z_vec[6] = z7;  z_vec[7] = z8;  z_vec[8] = z9;
    Vector<double,6> vals(m_fit_params * z_vec);

    //  Max is at [x,y] where:
    //   dz/dx = 2ax + cy + d = 0
    //   dz/dy = 2by + cx + e = 0
    double denom = 4.0 * vals[0] * vals[1] - (vals[2] * vals[2]); // = 4ab - c^2
    if (fabs(denom) < 0.01) {
      return false;
    }
    Vector2f offset((vals[2] * vals[4] - 2.0 * vals[1] * vals[3]) / denom,   // ce - 2bd
                     (vals[2] * vals[3] - 2.0 * vals[0] * vals[4]) / denom); // cd - 2ae

    dx = offset[0];
    dy = offset[1];

    // Apply correction to dx and dy
    double sX = 0.34574;
    double sY = 0.38944;
    dx = erf(dx/(sX*sqrt(2.0))) / 2.0;
    dy = erf(dy/(sY*sqrt(2.0))) / 2.0;

    if (norm_2(Vector2(dx, dy)) >= 0.5) {
      double scale = norm_2(Vector2(dx, dy)) / 0.5;
      dx /= scale;
      dy /= scale;
    }

    return true;
  } // End find_peak()

private: // Variables

  Matrix<float,6,9> m_fit_params;

}; // End class ParabolaFit2d

/// Class to make counting up "on" mask pixel in a box more efficient.
/// - This class moves a box through an image and returns the fraction of valid
///   pixels at each position of the box.
/// - By operating iteratively this class avoids a large memory buffer allocation.
/// - It looks like we have some similar VW classes, but nothing exactly like this.
/// TODO: Currently only works with positive position offsets!
/// TODO: Should this be moved to a different location?
class IterativeMaskBoxCounter {
public:
  IterativeMaskBoxCounter(ImageView<uint8> const* right_image_mask,
                          Vector2i box_size)
    : m_right_image_mask(right_image_mask), m_box_size(box_size) {

    // Set position to a flag value
    m_curr_pos = Vector2i(-1, 0);
    m_box_area = box_size[0]*box_size[1];
    m_curr_sum = 0;
  }

  /// Compute the percentage of valid pixels from the next column over.
  double next_pixel() {
    // Handle the first pixel in each column
    if (m_curr_pos[0] < 0) {
      m_curr_pos[0] = 0;
      return recompute();
    }
    // Handle other pixels
    m_curr_pos[0] += 1;

    // Account for column we are "losing"
    m_curr_sum -= m_column_sums.front(); 
    m_column_sums.pop();

    // Sum values from the new column
    int new_col_sum = 0;
    int col = m_curr_pos[0] + m_box_size[0] - 1;

    for (int row=m_curr_pos[1]; row<m_curr_pos[1]+m_box_size[1]; ++row) {
      if (m_right_image_mask->operator()(col, row) > 0)
        ++new_col_sum;
    }
    m_column_sums.push(new_col_sum);
    m_curr_sum += new_col_sum;

    return static_cast<double>(m_curr_sum) / m_box_area;
  }

  /// Move to the next row of pixels.
  void advance_row() {
    // Update pixel position and set recompute flag
    m_curr_pos[0]  = -1;
    m_curr_pos[1] +=  1;
    m_curr_sum = 0; 
  }
private:

  ImageView<uint8> const* m_right_image_mask;
  std::queue<int> m_column_sums;
  int    m_curr_sum;
  double m_box_area;
  Vector2i m_box_size;
  Vector2i m_curr_pos;

  /// Init counts from scratch at the current position.
  /// - Returns the fraction of active pixels at this position.
  double recompute() {
    m_column_sums = std::queue<int>(); // Clear sums
    
    for (int col=m_curr_pos[0]; col<m_curr_pos[0]+m_box_size[0]; ++col) {
      int col_sum = 0;
      for (int row=m_curr_pos[1]; row<m_curr_pos[1]+m_box_size[1]; ++row) {
        if (m_right_image_mask->operator()(col,row) > 0){
          col_sum += 1;
        }
      }
      m_column_sums.push(col_sum);
      m_curr_sum += col_sum;
    }
    return static_cast<double>(m_curr_sum) / m_box_area;
  }

}; // End class IterativeMaskBoxCounter

/**
 Helper class to manage the rolling accumulation buffer temporary memory
  until the results are added to the final accumulation buffer.  Used by 
  two_trip_path_accumulation().
- This class needs to temporarily store two rows of the accumulation buffer
  with each row containing four directional pass results simultaneously.
- A dedicated class is needed to manage the complexity of maintaining and
  accessing this buffer when each pixel requires a different amount of
  space in the buffer depending on its search range.
*/
class MultiAccumRowBuffer {
public:

  /// Four or eight directions, or "passes", are processed at a time and are 
  ///  stored according to their index number here.
  enum PassIndex { TOP_LEFT  = 0,
                   TOP       = 1,
                   TOP_RIGHT = 2,
                   LEFT      = 3,
                   BOT_RIGHT = 0,
                   BOT       = 1,
                   BOT_LEFT  = 2,
                   RIGHT     = 3,
                   PASS_ONE  = 0,
                   PASS_TWO  = 1 };

  /// Return the number of elements in the buffer
  static size_t multi_buf_size(const SemiGlobalMatcher* parent_ptr,
                               const int  num_paths_in_pass,
                               const bool vertical,
                               double buf_size_factor) {

    size_t line_size = parent_ptr->m_num_output_cols;
    if (vertical)
      line_size = parent_ptr->m_num_output_rows;

    // Instantiate two single-row buffers that will be used to temporarily store
    //  accumulated cost info until it is no longer needed.
    // - Within each buffer, data is indexed in order [col][pass][disparity]
    // - The actual data size in the buffer will vary each line, so it is 
    //    initialized to be the maximum possible size.
    const size_t buffer_pixel_size = num_paths_in_pass*parent_ptr->m_num_disp;
    size_t multi_buf_size = line_size*buffer_pixel_size;
    
    // No reason for the buffer size to be larger than the entire accumulator!
    if (multi_buf_size > parent_ptr->m_buffer_lengths)
      multi_buf_size = parent_ptr->m_buffer_lengths;

    // TODO(oalexan1): Two users reported failure that can be traced back to
    // "Ran out of memory in the small buffer", and ultimately to choices made
    // here. These choices are arbitrary. Why can't the buffer be made larger or
    // grow with m_buffer_lengths? This logic is repeated in a different place
    // in this file. We should allow the buffer size to be max of 
    // SAFE_BUFFER_SIZE and a percentage of the entire accumulation buffer.

    // If the buffer is over 128 MB, reduce its size to a percentage of the
    // size of the entire accumulation buffer.
    const size_t SAFE_BUFFER_SIZE = (1024*1024*128) / sizeof(SemiGlobalMatcher::AccumCostType);
    const double MAX_PERCENTAGE   = 0.04;

    std::cout << "1parent ptr buff len: " << parent_ptr->m_buffer_lengths << std::endl;
    std::cout << "Multi Buffer size is        " << multi_buf_size << std::endl;
    std::cout << "Multi Safe buffer size is   " << SAFE_BUFFER_SIZE << std::endl;
    
    // TODO(oalexan1): Must use here instead max of safe buffer and the percentage.
    // But must test.
    if (multi_buf_size > SAFE_BUFFER_SIZE) {
      multi_buf_size = parent_ptr->m_buffer_lengths * MAX_PERCENTAGE;
      std::cout << "---1 Multi will reduce buffer size to " << multi_buf_size << std::endl;
      if (multi_buf_size < SAFE_BUFFER_SIZE)
        multi_buf_size = SAFE_BUFFER_SIZE; // Buffer can at least be this size
    }
    std::cout << "--Multi final buffer size is " << multi_buf_size << std::endl;

    // This is a bugfix. Adjust the buffer size if a previous invocation failed.
    multi_buf_size *= buf_size_factor;
    
    return multi_buf_size;
  }

  /// Construct the buffers.
  /// - num_paths_in_pass can be 1, four (8 directions), or eight (16 directions)
  MultiAccumRowBuffer(const SemiGlobalMatcher* parent_ptr,
                      const int  num_paths_in_pass,
                      const bool vertical, 
                      double buf_size_factor) {
    m_parent_ptr        = parent_ptr;  
    m_num_paths_in_pass = num_paths_in_pass;
    m_vertical          = vertical;

    const int num_rows = m_parent_ptr->m_num_output_rows;
    m_line_size = m_parent_ptr->m_num_output_cols;
    if (vertical)
      m_line_size = num_rows;

    m_multi_buf_size = multi_buf_size(parent_ptr, num_paths_in_pass, vertical,
                                   buf_size_factor);
    m_multi_buf_size_bytes = m_multi_buf_size*sizeof(SemiGlobalMatcher::AccumCostType);

    vw_out(DebugMessage, "stereo") << "MultiAccumRowBuffer - allocating buffer size (MB): " 
                                   << m_multi_buf_size_bytes/(1024*1024) << std::endl;

    // Allocate buffers that store accumulation scores
    m_bufferA.reset(new SemiGlobalMatcher::AccumCostType[m_multi_buf_size]);
    m_bufferB.reset(new SemiGlobalMatcher::AccumCostType[m_multi_buf_size]);
    m_trail_buffer = m_bufferA.get();
    m_lead_buffer  = m_bufferB.get();

    // Allocate buffers that store pixel offsets into the buffers we just allocated
    m_offsetsA.reset(new size_t[m_line_size]);
    m_offsetsB.reset(new size_t[m_line_size]);   
    m_offsets_lead  = m_offsetsA.get();
    m_offsets_trail = m_offsetsB.get();

    if (!vertical) { // horizontal
      // The first trip rasterizes top-left to bottom-right
      m_current_col = 0;
      m_current_row = 0;
      m_col_advance = 1;
      m_row_advance = 1;
    } else { // vertical
      // The first trip rasterizes bottom-left to top-right
      m_current_col = 0;
      m_current_row = num_rows-1;
      m_col_advance = 1;
      m_row_advance = -1;
    }

    // Set up lead buffers, trailing buffer is not used until the next row.
    memset(m_lead_buffer, 0, m_multi_buf_size_bytes); // Init this buffer to zero
    //std::fill(m_lead_buffer, m_lead_buffer+m_multi_buf_size, m_parent_ptr->get_bad_accum_val());
    fill_lead_offset_buffer();
  }

  /// Load buffer offsets into the lead buffer for the current row.
  void fill_lead_offset_buffer() {
    //  Convert offsets to be relative to the start of our row/col instead of
    //  from pixel (0,0).  This allows us easy access into our line buffers.
    //  Remember to multiply by the number of paths stored.
    // - Pixel info is always stored left to right, even on the second trip through the image
    const size_t* raw_offsets    = m_parent_ptr->m_buffer_starts.data();
    if (!m_vertical) { // horizontal
      size_t  new_lead_index = m_current_row*m_parent_ptr->m_num_output_cols;
      size_t  start_offset   = raw_offsets[new_lead_index]; // Offset of the first column
      for (int i=0; i<m_line_size; ++i) {
        m_offsets_lead[i] = (raw_offsets[new_lead_index+i] - start_offset) * m_num_paths_in_pass;
      }
    } else { // vertical
      // In the vertical case we need to rebuild a set of offsets to describe the column.
      size_t position = 0;
      for (int i=0; i<m_line_size; ++i) {
        m_offsets_lead[i] = position;
        size_t size = m_parent_ptr->get_num_disparities(m_current_col, i) * m_num_paths_in_pass;
        position += size;
      }
    }
  }

  /// Add the results in the leading buffer to the main class accumulation buffer.
  /// - The scores from each pass are added.
  void add_lead_buffer_to_accum() {
    Mutex::Lock locker(m_mutex); // Scoped lock so this function can't be run simultaneously

    size_t buffer_index = 0;
    SemiGlobalMatcher::AccumCostType* out_ptr = m_parent_ptr->m_accum_buffer.get();
    if (!m_vertical) { // horizontal
      for (int col=0; col<m_parent_ptr->m_num_output_cols; ++col) {
        int num_disp = m_parent_ptr->get_num_disparities(col, m_current_row);
        for (int pass=0; pass<m_num_paths_in_pass; ++pass) {
          size_t out_index = m_parent_ptr->m_buffer_starts(col, m_current_row);
          for (int d=0; d<num_disp; ++d) {
            out_ptr[out_index++] += m_lead_buffer[buffer_index++];
            //printf("row, col, pass, d = %d, %d, %d, %d ->> %d ->> %d\n", 
            //    m_current_row, col, pass, d, m_trail_buffer[buffer_index], m_parent_ptr->m_accum_buffer[out_index]);
          } // end disp loop
        } // end pass loop
      } // end col loop
    } else { // vertical
      for (int row=0; row<m_parent_ptr->m_num_output_rows; ++row) {
        int num_disp = m_parent_ptr->get_num_disparities(m_current_col, row);
        for (int pass=0; pass<m_num_paths_in_pass; ++pass) {
          size_t out_index = m_parent_ptr->m_buffer_starts(m_current_col, row);
          for (int d=0; d<num_disp; ++d) {
            out_ptr[out_index++] += m_lead_buffer[buffer_index++];
            //printf("row, col, pass, d = %d, %d, %d, %d ->> %d ->> %d\n", 
            //    m_current_row, col, pass, d, m_trail_buffer[buffer_index], m_parent_ptr->m_accum_buffer[out_index]);
          } // end disp loop
        } // end pass loop
      } // end col loop
    }
  } // end add_trail_buffer_to_accum

  /// Call when moving to the next pixel in a column
  void next_pixel() {
    if (!m_vertical) // horizontal
      m_current_col += m_col_advance;
    else // vertical
      m_current_row += m_row_advance;
  }

  /// Call when finished processing a column.
  void next_row(bool trip_finished) {

    // The first thing to do is to record the results from the last row
    add_lead_buffer_to_accum();

    if (trip_finished) // Quit early if the current trip is finished
      return;

    if (!m_vertical) { // horizontal
      // Update the position in the image
      m_current_row += m_row_advance;
      if (m_row_advance > 0) // First pass
        m_current_col = 0;
      else // Second pass
        m_current_col = m_parent_ptr->m_num_output_cols - 1;
    } else { // vertical
      // Update the position in the image
      m_current_col += m_col_advance;
      if (m_col_advance > 0) // First pass
        m_current_row = m_parent_ptr->m_num_output_rows - 1;
      else // Second pass
        m_current_row = 0;
    }

    // Swap accum buffer pointers and init the lead buffer
    std::swap(m_trail_buffer, m_lead_buffer);
    //std::fill(m_lead_buffer, m_lead_buffer+m_multi_buf_size, m_parent_ptr->get_bad_accum_val());
    memset(m_lead_buffer, 0, m_multi_buf_size_bytes); // Init this buffer to zero

    // Swap offset buffer pointers and init the lead buffer
    std::swap(m_offsets_trail, m_offsets_lead);
    fill_lead_offset_buffer();
  }

  /// Call after the first trip is finished before the second pass.
  void switch_trips() {

    if (!m_vertical) { // horizontal
      // Second trip is from bottom right to top left
      m_current_col = m_parent_ptr->m_num_output_cols - 1;
      m_current_row = m_parent_ptr->m_num_output_rows - 1;
      m_col_advance = -1; 
      m_row_advance = -1;
    } else { // vertical
      // Second trip is from top right to bottom left
      m_current_col = m_parent_ptr->m_num_output_cols - 1;
      m_current_row = 0;
      m_col_advance = -1; 
      m_row_advance = 1;
    }

    // Set up lead buffers, trailing buffer is not used until the next row.
    memset(m_lead_buffer, 0, m_multi_buf_size_bytes);    
    fill_lead_offset_buffer();
  }

  /// Get the pointer to write the output of the current pass to
  SemiGlobalMatcher::AccumCostType * get_output_accum_ptr(PassIndex pass) {
    
    int    num_disp    = m_parent_ptr->get_num_disparities(m_current_col, m_current_row);
    size_t pass_offset = num_disp*pass;
    size_t offset      = 0;
    if (m_vertical)
      offset = (m_offsets_lead[m_current_row] + pass_offset);
    else
      offset = (m_offsets_lead[m_current_col] + pass_offset);

    // Make sure there is enough memory left to support this location.
    // TODO(oalexan1): Look into this too.
    if (offset + num_disp > m_multi_buf_size) {
      std::cout << "--2 will throw: Insufficient memory in small buffers, disparity "
                << "image may be degenerate.\n";
      vw_throw(ArgumentErr() << "Insufficient memory in small buffers, "
                             << "disparity image may be degenerate.\n");
    }

    return m_lead_buffer + offset;
  }

  /// Gets the pointer to the accumulation buffer for the indicated pixel/pass
  SemiGlobalMatcher::AccumCostType*
   get_trailing_pixel_accum_ptr(int col_offset, int row_offset, PassIndex pass) {
 
    int    col         = m_current_col + col_offset;
    int    row         = m_current_row + row_offset;
    int    num_disp    = m_parent_ptr->get_num_disparities(col, row);
    size_t pass_offset = num_disp*pass;

    // Just get the index of the column in the correct buffer, then add an offset for the selected pass.
    size_t offset = 0;
    SemiGlobalMatcher::AccumCostType * output_ptr=0;
    if (!m_vertical) { // horizontal
      if (row_offset == 0) { 
        // Same row, must be the in the leading buffer
        offset     = m_offsets_lead[col] + pass_offset;
        output_ptr = m_lead_buffer;
      } else { 
        // Different row, must be in the trailing buffer
        offset     = m_offsets_trail[col] + pass_offset;
        output_ptr = m_trail_buffer;
      }
    } else { // vertical
      if (col_offset == 0) { 
        // Same row, must be the in the leading buffer
        offset     = m_offsets_lead[row] + pass_offset;
        output_ptr = m_lead_buffer;
      } else { 
        // Different row, must be in the trailing buffer
        offset     = m_offsets_trail[row] + pass_offset;
        output_ptr = m_trail_buffer;
      }
    }

    // Make sure there is enough memory left to support this location.
    if (offset + num_disp > m_multi_buf_size) {
      std::cout << "--3 will throw Insufficient memory in small buffers, disparity "
                << "image may be degenerate.\n";
      vw_throw(ArgumentErr() << "Insufficient memory in small buffers, disparity "
                             << "image may be degenerate.\n");
    }
    
    return output_ptr + offset;

  } // End function get_trailing_pixel_accum_ptr

private:

  /// Limit main accumulation buffer write access across any number of threads.
  static Mutex m_mutex; 

  const SemiGlobalMatcher* m_parent_ptr; ///< Need a handle to the parent SGM object

  bool m_vertical; ///< Raster orientation
  int  m_num_paths_in_pass; ///< Must be 4 or 8

  int    m_line_size; ///< Length of a column(horizontal) or a row(vertical) in pixels.
  size_t m_multi_buf_size, m_multi_buf_size_bytes;
  int    m_current_col, m_current_row; ///< The current position as we iterate through the pixels
  int    m_col_advance, m_row_advance; ///< These are set according to the current trip

  // These point to m_buffer_starts for the leading and trailing row respectively.
  // - Since there are N passes per pixel, the offsets are different than in the main class
  //   accumulation buffer.
  boost::shared_array<size_t> m_offsetsA, m_offsetsB;
  size_t * m_offsets_lead;  // The role of the buffers keeps swapping so these pointers are
  size_t * m_offsets_trail; //  used to keep things consistent.
  
  // Buffers which store the accumulated cost info before it is dumped to the main accum buffer
  boost::shared_array<SemiGlobalMatcher::AccumCostType> m_bufferA, m_bufferB;
  SemiGlobalMatcher::AccumCostType* m_trail_buffer; // Another set of pointers for swapping buffers.
  SemiGlobalMatcher::AccumCostType* m_lead_buffer;

}; // End class MultiAccumRowBuffer

Mutex MultiAccumRowBuffer::m_mutex;

/**
  A single line, single pass SGM accumulation buffer designed to be used by a single thread.
  - Each line-pass instance of the SGM accumulation problem is independent until the end.
  - Instances of this class must share access to the parent class accumulation buffer so that
    only one instance of the class touches it at any one time.
  - This class is much less complicated than the multi-line buffer class!
*/
class OneLineBuffer {
public:
  /// Default constructor does not initialize the object
  OneLineBuffer() {}

  /// Initializing constructor
  OneLineBuffer(SemiGlobalMatcher const* parent_ptr) {
    initialize(parent_ptr);
  }

  /// Returns the size it elements of the buffer.
  static size_t one_buf_size(const SemiGlobalMatcher* parent_ptr) {

    // Figure out the max possible line length (diagonal line down the center)    
    const int num_cols = parent_ptr->m_num_output_cols;
    const int num_rows = parent_ptr->m_num_output_rows;
    int line_size = sqrt(num_cols*num_cols + num_rows*num_rows) + 1;

    // Instantiate single-row buffer that will be used to temporarily store
    //  accumulated cost info until it can be added to the main SGM class buffer.
    // - Within each buffer, data is indexed in order [pixel][disparity]
    // - The actual data size in the buffer will vary each line, so it is 
    //    initialized to be the maximum possible size.
    size_t one_buf_size = line_size*parent_ptr->m_num_disp;

    // No reason for the buffer size to be larger than the entire accumulator!
    if (one_buf_size > parent_ptr->m_buffer_lengths)
      one_buf_size = parent_ptr->m_buffer_lengths;

    // TODO(oalexan1): Two users reported failure that can be traced back to
    // "Ran out of memory in the small buffer", and ultimately to choices made
    // here. These choices are arbitrary. Why can't the buffer be made larger or
    // grow with m_buffer_lengths? This logic is repeated in a different place
    // in this file. We should allow the buffer size to be max of
    // SAFE_BUFFER_SIZE and a percentage of the entire accumulation buffer.

    // TODO(oalexan1): Must use here max of 128 MB and the percentage. But must test.
    // If the buffer is over 64 MB, reduce its size to a percentage of the
    //  size of the entire accumulation buffer.
    size_t SAFE_BUFFER_SIZE = (1024*1024*64) / sizeof(SemiGlobalMatcher::AccumCostType);
    double MAX_PERCENTAGE   = 0.02;

    std::cout << "one 2parent ptr buff len: " << parent_ptr->m_buffer_lengths << std::endl;
    std::cout << "one Buffer size is        " << one_buf_size << std::endl;
    std::cout << "one Safe buffer size is   " << SAFE_BUFFER_SIZE << std::endl;

    if (one_buf_size > SAFE_BUFFER_SIZE) {
      one_buf_size = parent_ptr->m_buffer_lengths * MAX_PERCENTAGE;
      std::cout << "--one 2 will reduce buffer size to " << one_buf_size << std::endl;
      if (one_buf_size < SAFE_BUFFER_SIZE)
        one_buf_size = SAFE_BUFFER_SIZE; // Buffer can at least be this size
    }
    std::cout << "--one final buffer size is " << one_buf_size << std::endl;

    // This is a bugfix. Adjust the buffer size if a previous invocation failed.
    one_buf_size *= parent_ptr->m_buf_size_factor;

    return one_buf_size;
  }

  /// Initialize the buffer.
  /// - All this really does is allocate a memory buffer of the maximum possible required size.
  void initialize(const SemiGlobalMatcher* parent_ptr) {

    // Determine the buffer size
    m_num_disp          = parent_ptr->m_num_disp;
    m_one_buf_size       = one_buf_size(parent_ptr);
    m_one_buf_size_bytes = m_one_buf_size*sizeof(SemiGlobalMatcher::AccumCostType);

    vw_out(DebugMessage, "stereo") << "OneLineBuffer - allocating buffer size (MB): " 
                                   << m_one_buf_size_bytes/(1024*1024) << std::endl;

    // Allocate the accumulation buffer
    m_buffer.reset(new SemiGlobalMatcher::AccumCostType[m_one_buf_size]);

    // Set up the small buffer
    m_bad_disp_value = parent_ptr->get_bad_accum_val();
    m_full_prior_buffer.reset(new SemiGlobalMatcher::AccumCostType[m_num_disp]);
  }

  /// Clear both buffers
  void clear_buffers() {
    memset(m_buffer.get(), 0, m_one_buf_size_bytes);

    for (size_t i=0; i<m_num_disp; ++i)
      m_full_prior_buffer[i] = m_bad_disp_value;
  }

  /// Get the pointer to the start of the output accumulation buffer
  SemiGlobalMatcher::AccumCostType * get_output_accum_buf_ptr(size_t &one_buf_size) {
    one_buf_size = m_one_buf_size;
    return m_buffer.get();
  }

  // Return the smaller buffer
  SemiGlobalMatcher::AccumCostType * get_full_prior_ptr() {
    return m_full_prior_buffer.get();
  }
  
private: // Variables

  SemiGlobalMatcher::AccumCostType m_bad_disp_value;
  size_t m_one_buf_size, m_one_buf_size_bytes, m_num_disp;

  /// Buffer which store the accumulated cost info before it is dumped to the main accum buffer
  boost::shared_array<SemiGlobalMatcher::AccumCostType> m_buffer;

  /// Much smaller buffer needed by the pixel evaluation function
  boost::shared_array<SemiGlobalMatcher::AccumCostType> m_full_prior_buffer;

}; // End class OneLineBuffer





/// Assigns OneLineBuffer objects to processing threads so they do not conflict.
class OneLineBufferManager {
public:
  OneLineBufferManager(int num_threads, SemiGlobalMatcher * parent_ptr) {
  
    // Initialize a number of line buffers equal to the number of threads
    m_buffer_vec.resize(num_threads);
    m_busy_vec.resize  (num_threads);
    for (int i=0; i<num_threads; ++i) {
      // Init all buffers and mark as unused.
      m_buffer_vec[i].initialize(parent_ptr);
      m_busy_vec  [i] = false;
    }
  }

  // Request the ID of a free buffer
  size_t get_free_buffer_id() {
    Mutex::Lock locker(m_main_mutex); // Scoped lock so this function can't be run simultaneously

    // Find the next free buffer
    for (size_t i=0; i<m_busy_vec.size(); ++i) {
      if (m_busy_vec[i] == false) {
        m_busy_vec[i] = true; // Mark the buffer as in use and return the index 
        return i;
      }
    }
    vw_throw(LogicErr() << "Error: Should always be a free buffer!\n");
  }

  /// Once we have an ID, get the actual buffer
  OneLineBuffer* get_line_buffer(size_t id) {
    // Clear the buffers before we hand it off
    m_buffer_vec[id].clear_buffers();
    return &(m_buffer_vec[id]);
  }

  /// Notify the manager that we are finished using a buffer.
  void release_buffer(size_t id) {
    Mutex::Lock locker(m_main_mutex); // Scoped lock so this function can't be run simultaneously
    m_busy_vec[id] = false; // Mark the buffer as no longer in use
  }

private:
  Mutex m_main_mutex; ///< Control access to this class
  std::vector<OneLineBuffer> m_buffer_vec; ///< Available pre-allocated memory buffers.
  std::vector<bool         > m_busy_vec;   ///< A mutex to lock each of the buffers.
}; // End class OneLineBufferManager




/// Performs SGM accumulation along one line and then adds it to the parent accumulation buffer
/// - Task object which can be passed to a thread pool.
class PixelPassTask: public Task {
public:

  /// Constructor, pass in all the needed links.
  PixelPassTask(ImageView<uint8>  const* image_ptr,
                SemiGlobalMatcher      * parent_ptr,
                OneLineBufferManager   * buffer_manager_ptr,
                PixelLineIterator pixel_loc_iter, 
                int * success): 
  m_image_ptr(image_ptr), m_parent_ptr(parent_ptr), 
  m_buffer_manager_ptr(buffer_manager_ptr), m_pixel_loc_iter(pixel_loc_iter),
  m_success(success) {}

  /// Do the work in a thread
  void PixelPassDoWork() {

    // Retrieve a memory buffer to work with
    size_t buffer_id = m_buffer_manager_ptr->get_free_buffer_id();
    OneLineBuffer* buff_ptr       = m_buffer_manager_ptr->get_line_buffer(buffer_id);
    AccumCostType* full_prior_ptr = buff_ptr->get_full_prior_ptr();

    // Make a copy of the pixel iterator so we can re-use it for accum buffer addition
    PixelLineIterator pixel_loc_iter_copy(m_pixel_loc_iter);

    //const bool debug = false;
    int last_pixel_val = -1;
    int col_prev = -1, row_prev = -1; // Previous row and column

    // Get the start of the output accumulation buffer
    // - Storage here is simply num_disp for each pixel in the line, one after the other.
    size_t buffer_size   = 0;
    size_t consumed_size = 0;
    AccumCostType* computed_accum_ptr = buff_ptr->get_output_accum_buf_ptr(buffer_size);
    AccumCostType* prior_accum_ptr    = 0;

    while (pixel_loc_iter_copy.is_good()) {

      // Get current location in output image
      const int col = pixel_loc_iter_copy.col();
      const int row = pixel_loc_iter_copy.row();

      // Get corresponding location in the input image
      const int input_col = col + m_parent_ptr->m_min_col;
      const int input_row = row + m_parent_ptr->m_min_row;

      bool debug = false;//(col==152) && (row ==12);

      // Get information about the current pixel location from the parent
      int num_disp = m_parent_ptr->get_num_disparities(col, row);
      CostType * const local_cost_ptr = m_parent_ptr->get_cost_vector(col, row);

      // Make sure we don't run out of memory in the buffer
      // TODO(oalexan1): This is the problem
      consumed_size += num_disp;
      if (consumed_size > buffer_size) {
        std::cout << "---1will throw: Ran out of memory in the small buffer, disparity image may be degenerate.\n";
        vw_throw(ArgumentErr() << "Ran out of memory in the small buffer, "
                               << "disparity image may be degenerate.\n");
      }

      // Fill in the accumulated value in the bottom buffer
      int curr_pixel_val = static_cast<int>(m_image_ptr->operator()(input_col, input_row));
      int pixel_diff     = std::abs(curr_pixel_val - last_pixel_val);

      //printf("Loc %d, %d, diff = %d, num_disp = %d\n", col, row, pixel_diff, num_disp);
      //std::cout << "DEBUG " << m_parent_ptr->m_disp_bound_image(col,row) << std::endl;

      if (last_pixel_val >= 0) { // All pixels after the first
        m_parent_ptr->evaluate_path(col, row, col_prev, row_prev,
                                    prior_accum_ptr, full_prior_ptr, local_cost_ptr, computed_accum_ptr, 
                                    pixel_diff, debug);
      } else { // First pixel only, nothing to accumulate.
        for (int d = 0; d < num_disp; d++)
          computed_accum_ptr[d] = local_cost_ptr[d];
      }

      // Advance the position
      prior_accum_ptr = computed_accum_ptr; // Retain the current accumulation buffer location
      computed_accum_ptr += num_disp;       // Advance to the next accumulation buffer location
      last_pixel_val  = curr_pixel_val;     // Retain the current pixel value
      col_prev        = col;                // Retain the current pixel location
      row_prev        = row;
      pixel_loc_iter_copy++;                // Update the pixel location
    } // End loop through pixels

    // Now that we computed all the results, add them to the main buffer 
    update_accum_buffer(buff_ptr);

    // Notify the buffer manager that we are finished with the buffer
    //std::cout << "Releasing buffer " << buffer_id << std::endl;
    m_buffer_manager_ptr->release_buffer(buffer_id);
  } // end PixelPassDoWork()
  
  // Do the work while catching any exceptions in the thread, as those can abort the
  // entire program.
  virtual void operator()() {
    
    *m_success = 1; // Assume success until proven otherwise
    
    try {
      PixelPassDoWork();
    } catch (const std::exception& e) {
      *m_success = 0;
    }

  } // End operator() function

  /// Add the computed buffer results to the parent accumulation buffer
  void update_accum_buffer(OneLineBuffer * buff_ptr) {

    // Get the start of the output accumulation buffer
    // - Storage here is simply num_disp for each pixel in the line, one after the other.
    size_t buffer_size=0;
    AccumCostType* computed_accum_ptr = buff_ptr->get_output_accum_buf_ptr(buffer_size);

    // Loop through all pixels in the line
    while (m_pixel_loc_iter.is_good()) {

      // Get current location
      const int col = m_pixel_loc_iter.col();
      const int row = m_pixel_loc_iter.row();

      // Get information about the current pixel location from the parent
      int num_disp = m_parent_ptr->get_num_disparities(col, row);

      AccumCostType* output_accum_ptr = m_parent_ptr->get_accum_vector(col, row);

      // Add the computed values to the output accumulation location      
      for (int i=0; i<num_disp; ++i) {
        output_accum_ptr[i] += computed_accum_ptr[i];
      }
      // Advance through pixel position and the computed data buffer
      computed_accum_ptr += num_disp;
      m_pixel_loc_iter++;
    }
  } // End update_accum_buffer() function

private:

  typedef SemiGlobalMatcher::AccumCostType AccumCostType;
  typedef SemiGlobalMatcher::CostType      CostType;

  ImageView<uint8>  const* m_image_ptr;
  SemiGlobalMatcher      * m_parent_ptr;
  OneLineBufferManager   * m_buffer_manager_ptr;
  PixelLineIterator m_pixel_loc_iter; ///< Keeps track of the pixel position
  int * m_success; // Will be set to zero if the task fails 

}; // End class PixelPassTask

/// Task wrapper for each of the required smooth accumulation passes.
class SmoothPathAccumTask : public Task {

public:

  // Used to specify the accumulation direction
  enum Direction {TL, T, TR, L, R, BL, B, BR};

  /// Constructor
  SmoothPathAccumTask(MultiAccumRowBuffer    * buffer_ptr,
                      SemiGlobalMatcher      * parent_ptr,
                      ImageView<uint8>  const* image_ptr,
                      Direction dir, 
                      int * success):
  m_buffer_ptr(buffer_ptr), m_parent_ptr(parent_ptr),
  m_image_ptr(image_ptr), m_dir(dir), m_success(success) {

    *m_success = 1; // Assume success until proven otherwise
   
    // Init this buffer to bad scores representing disparities that were
    //  not in the search range for the given pixel. 
    m_full_prior_buffer.reset(new AccumCostType[parent_ptr->m_num_disp]);
    for (int i=0; i<parent_ptr->m_num_disp; ++i)
      m_full_prior_buffer[i] = parent_ptr->get_bad_accum_val();  

    // Allocate a buffer for the "perpendicular direction" results to be written to
    m_temp_buffer.reset(new AccumCostType[parent_ptr->m_num_disp]);

    m_last_column = parent_ptr->m_num_output_cols - 1;
    m_last_row    = parent_ptr->m_num_output_rows - 1;
  }

  /// Main task function redirects to the dedicated function
  virtual void operator()() {
    
    // Must catch exceptions and track success
    try {
      switch (m_dir) {
      case TL: task_TL(); return;
      case T:  task_T (); return;
      case TR: task_TR(); return;
      case L:  task_L (); return;
      case R:  task_R (); return;
      case BL: task_BL(); return;
      case B:  task_B (); return;
      default: task_BR(); return;
      };
    } catch (...) {
      *m_success = 0;
    }
  }

private: // Variables

  typedef SemiGlobalMatcher::AccumCostType AccumCostType;
  typedef SemiGlobalMatcher::CostType      CostType;

  MultiAccumRowBuffer    * m_buffer_ptr; ///< Two-line buffer to be used for this task
  SemiGlobalMatcher      * m_parent_ptr; ///< Pointer to the main SGM class
  ImageView<uint8>  const* m_image_ptr;
  Direction m_dir;                       ///< Direction of this task

  int m_last_column, m_last_row;

  boost::shared_array<AccumCostType> m_full_prior_buffer; ///< Working buffer needed for evaluate_task() function
  boost::shared_array<AccumCostType> m_temp_buffer;       ///< Buffer for storing the smoothed-in result.
  int * m_success; /// Will have a non-zero value if the task was successful
private: // Functions - one per direction

  /// Helper function to save some lines.
  int get_num_disp(int col, int row) {
    int num_disp = m_parent_ptr->get_num_disparities(col, row);
    if (num_disp == 0)
      m_buffer_ptr->next_pixel(); // In preparation for skipping the loop iteration
    return num_disp;
  }

  void task_L() {
    
    AccumCostType* full_prior_ptr = m_full_prior_buffer.get();
    AccumCostType* output_accum_ptr;

    // Loop through all pixels in the output image for the L trip, top-left to bottom-right.
    for (int row=0; row<m_parent_ptr->m_num_output_rows; ++row) {
      for (int col=0; col<m_parent_ptr->m_num_output_cols; ++col) {

        //printf("Accum pass 1 col = %d, row = %d\n", col, row);

        // TODO: Do we need to skip these pixels in the two-pass SGM function above?
        int num_disp = get_num_disp(col, row);
        if (num_disp == 0)
          continue;
        CostType * const local_cost_ptr = m_parent_ptr->get_cost_vector(col, row);
        bool debug = false;//((row == 244) && (col == 341));

        // Left
        output_accum_ptr = m_buffer_ptr->get_output_accum_ptr(MultiAccumRowBuffer::PASS_ONE);
        if ((row > 0) && (col > 0)) {
          int pixel_diff = m_parent_ptr->get_path_pixel_diff(*m_image_ptr, col, row, -1, 0);
          // Compute accumulation from the values in the left pixel
          AccumCostType* const prior_accum_ptr = m_buffer_ptr->get_trailing_pixel_accum_ptr(-1, 0, MultiAccumRowBuffer::PASS_ONE);
          m_parent_ptr->evaluate_path(col, row, col-1, row,
                                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                                       pixel_diff, debug);

          // Compute accumulation from the values in the above pixel
          AccumCostType* const prior_accum_ptr2 = m_buffer_ptr->get_trailing_pixel_accum_ptr(0, -1, MultiAccumRowBuffer::PASS_ONE);
          m_parent_ptr->evaluate_path(col, row, col, row-1,
                                       prior_accum_ptr2, full_prior_ptr, local_cost_ptr, m_temp_buffer.get(), 
                                       pixel_diff, debug);
          // The final accumulation values are the average of the two computations
          for (int d=0; d<num_disp; ++d)
            output_accum_ptr[d] = (output_accum_ptr[d] + m_temp_buffer[d])/2;
        }
        else // Just init to the local cost
          for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

        m_buffer_ptr->next_pixel();
      } // End col loop

      m_buffer_ptr->next_row(row==m_last_row);
    } // End row loop

  } // End task_TL

  void task_TL() {

    AccumCostType* full_prior_ptr = m_full_prior_buffer.get();
    AccumCostType* output_accum_ptr;

    // Loop through all pixels in the output image for the TL trip, top-left to bottom-right.
    for (int row=0; row<m_parent_ptr->m_num_output_rows; ++row) {
      for (int col=0; col<m_parent_ptr->m_num_output_cols; ++col) {

        //printf("Accum pass 1 col = %d, row = %d\n", col, row);

        // TODO: Do we need to skip these pixels in the two-pass SGM function above?
        int num_disp = get_num_disp(col, row);
        if (num_disp == 0)
          continue;
        CostType * const local_cost_ptr = m_parent_ptr->get_cost_vector(col, row);
        bool debug = false;//((row == 244) && (col == 341));

        // Top left
        output_accum_ptr = m_buffer_ptr->get_output_accum_ptr(MultiAccumRowBuffer::PASS_ONE);
        if ((row > 0) && (col > 0) && (col < m_last_column)) {

          int pixel_diff = m_parent_ptr->get_path_pixel_diff(*m_image_ptr, col, row, -1, -1);
          AccumCostType* const prior_accum_ptr = m_buffer_ptr->get_trailing_pixel_accum_ptr(-1, -1, MultiAccumRowBuffer::PASS_ONE);
          m_parent_ptr->evaluate_path(col, row, col-1, row-1,
                                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                                       pixel_diff, debug);

          AccumCostType* const prior_accum_ptr2 = m_buffer_ptr->get_trailing_pixel_accum_ptr(1, -1, MultiAccumRowBuffer::PASS_ONE);
          m_parent_ptr->evaluate_path(col, row, col+1, row-1,
                                       prior_accum_ptr2, full_prior_ptr, local_cost_ptr, m_temp_buffer.get(), 
                                       pixel_diff, debug);
          for (int d=0; d<num_disp; ++d)
            output_accum_ptr[d] = (output_accum_ptr[d] + m_temp_buffer[d])/2;
        }
        else // Just init to the local cost
          for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

        m_buffer_ptr->next_pixel();
      } // End col loop

      m_buffer_ptr->next_row(row==m_last_row);
    } // End row loop
  } // End task_TL

  void task_R() {
    AccumCostType* full_prior_ptr = m_full_prior_buffer.get();
    AccumCostType* output_accum_ptr;

    // Loop through all pixels in the output image for the R trip, bottom-right to top-left.
    for (int row = m_last_row; row >= 0; --row) {
      for (int col = m_last_column; col >= 0; --col) {

        int num_disp = get_num_disp(col, row);
        if (num_disp == 0)
          continue;
        CostType * const local_cost_ptr = m_parent_ptr->get_cost_vector(col, row);
        bool debug = false;//((row == 244) && (col == 341));

        // Right
        output_accum_ptr = m_buffer_ptr->get_output_accum_ptr(MultiAccumRowBuffer::PASS_ONE);
        if ((row < m_last_row) && (col < m_last_column)) {
          int pixel_diff = m_parent_ptr->get_path_pixel_diff(*m_image_ptr, col, row, 1, 0);
          AccumCostType* const prior_accum_ptr = m_buffer_ptr->get_trailing_pixel_accum_ptr(1, 0, MultiAccumRowBuffer::PASS_ONE);
          m_parent_ptr->evaluate_path(col, row, col+1, row,
                                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                                       pixel_diff, debug);

          AccumCostType* const prior_accum_ptr2 = m_buffer_ptr->get_trailing_pixel_accum_ptr(0, 1, MultiAccumRowBuffer::PASS_ONE);
          m_parent_ptr->evaluate_path(col, row, col, row+1,
                                       prior_accum_ptr2, full_prior_ptr, local_cost_ptr, m_temp_buffer.get(), 
                                       pixel_diff, debug);
          for (int d=0; d<num_disp; ++d)
            output_accum_ptr[d] = (output_accum_ptr[d] + m_temp_buffer[d])/2;                      
        }
        else // Just init to the local cost
          for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

        m_buffer_ptr->next_pixel();
      } // End col loop

      m_buffer_ptr->next_row(row==0);    
    } // End row loop
  } // End task_R

  void task_BR() {
    AccumCostType* full_prior_ptr = m_full_prior_buffer.get();
    AccumCostType* output_accum_ptr;

    // Loop through all pixels in the output image for the second trip, bottom-right to top-left.
    for (int row = m_last_row; row >= 0; --row) {
      for (int col = m_last_column; col >= 0; --col) {

        int num_disp = get_num_disp(col, row);
        if (num_disp == 0)
          continue;
        CostType * const local_cost_ptr = m_parent_ptr->get_cost_vector(col, row);
        bool debug = false;//((row == 244) && (col == 341));

        // Bottom right
        output_accum_ptr = m_buffer_ptr->get_output_accum_ptr(MultiAccumRowBuffer::PASS_ONE);
        if ((row < m_last_row) && (col > 0) && (col < m_last_column)) {

          int pixel_diff = m_parent_ptr->get_path_pixel_diff(*m_image_ptr, col, row, 1, 1);
          AccumCostType* const prior_accum_ptr = m_buffer_ptr->get_trailing_pixel_accum_ptr(1, 1, MultiAccumRowBuffer::PASS_ONE);
          m_parent_ptr->evaluate_path(col, row, col+1, row+1,
                                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                                       pixel_diff, debug);

          AccumCostType* const prior_accum_ptr2 = m_buffer_ptr->get_trailing_pixel_accum_ptr(-1, 1, MultiAccumRowBuffer::PASS_ONE);
          m_parent_ptr->evaluate_path(col, row, col-1, row+1,
                                       prior_accum_ptr2, full_prior_ptr, local_cost_ptr, m_temp_buffer.get(), 
                                       pixel_diff, debug);
          for (int d=0; d<num_disp; ++d)
            output_accum_ptr[d] = (output_accum_ptr[d] + m_temp_buffer[d])/2;
        }
        else // Just init to the local cost
          for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

        m_buffer_ptr->next_pixel();
      } // End col loop

      m_buffer_ptr->next_row(row==0);    
    } // End row loop 
  } // End task_BR

  void task_B() {
    AccumCostType* full_prior_ptr = m_full_prior_buffer.get();
    AccumCostType* output_accum_ptr;

    // Loop through all pixels in the output image for the third trip, bottom-left to top-right.
    for (int col = 0; col < m_parent_ptr->m_num_output_cols; ++col) {
      for (int row = m_last_row; row >= 0; --row) {

        int num_disp = get_num_disp(col, row);
        if (num_disp == 0)
          continue;
        CostType * const local_cost_ptr = m_parent_ptr->get_cost_vector(col, row);
        bool debug = false;//((row == 244) && (col == 341));

        // Bottom
        output_accum_ptr = m_buffer_ptr->get_output_accum_ptr(MultiAccumRowBuffer::PASS_ONE);
        if ((row < m_last_row) && (col > 0)) {
          int pixel_diff = m_parent_ptr->get_path_pixel_diff(*m_image_ptr, col, row, 0, 1);
          AccumCostType* const prior_accum_ptr = m_buffer_ptr->get_trailing_pixel_accum_ptr(0, 1, MultiAccumRowBuffer::PASS_ONE);
          m_parent_ptr->evaluate_path(col, row, col, row+1,
                                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                                       pixel_diff, debug);

          AccumCostType* const prior_accum_ptr2 = m_buffer_ptr->get_trailing_pixel_accum_ptr(-1, 0, MultiAccumRowBuffer::PASS_ONE);
          m_parent_ptr->evaluate_path(col, row, col-1, row,
                                       prior_accum_ptr2, full_prior_ptr, local_cost_ptr, m_temp_buffer.get(), 
                                       pixel_diff, debug);
          for (int d=0; d<num_disp; ++d)
            output_accum_ptr[d] = (output_accum_ptr[d] + m_temp_buffer[d])/2;
        }
        else // Just init to the local cost
          for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

        m_buffer_ptr->next_pixel();
      } // End col loop

      m_buffer_ptr->next_row(col==m_last_column);    
    } // End row loop 
  } // End task_B

  void task_BL() {
    AccumCostType* full_prior_ptr = m_full_prior_buffer.get();
    AccumCostType* output_accum_ptr;

    // Loop through all pixels in the output image for the third trip, bottom-left to top-right.
    for (int col = 0; col < m_parent_ptr->m_num_output_cols; ++col) {
      for (int row = m_last_row; row >= 0; --row) {

        int num_disp = get_num_disp(col, row);
        if (num_disp == 0)
          continue;
        CostType * const local_cost_ptr = m_parent_ptr->get_cost_vector(col, row);
        bool debug = false;//((row == 244) && (col == 341));

        // Bottom left
        output_accum_ptr = m_buffer_ptr->get_output_accum_ptr(MultiAccumRowBuffer::PASS_ONE);
        if ((row > 0) && (row < m_last_row) && (col > 0)) {
          int pixel_diff = m_parent_ptr->get_path_pixel_diff(*m_image_ptr, col, row, -1, 1);
          AccumCostType* const prior_accum_ptr = m_buffer_ptr->get_trailing_pixel_accum_ptr(-1, 1, MultiAccumRowBuffer::PASS_ONE);
          m_parent_ptr->evaluate_path(col, row, col-1, row+1,
                                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                                       pixel_diff, debug);

          AccumCostType* const prior_accum_ptr2 = m_buffer_ptr->get_trailing_pixel_accum_ptr(-1, -1, MultiAccumRowBuffer::PASS_ONE);
          m_parent_ptr->evaluate_path(col, row, col-1, row-1,
                                       prior_accum_ptr2, full_prior_ptr, local_cost_ptr, m_temp_buffer.get(), 
                                       pixel_diff, debug);
          for (int d=0; d<num_disp; ++d)
            output_accum_ptr[d] = (output_accum_ptr[d] + m_temp_buffer[d])/2;
        }
        else // Just init to the local cost
          for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

        m_buffer_ptr->next_pixel();
      } // End col loop

      m_buffer_ptr->next_row(col==m_last_column);    
    } // End row loop 
  } // End task_BL

  void task_T() {
    AccumCostType* full_prior_ptr = m_full_prior_buffer.get();
    AccumCostType* output_accum_ptr;

    // Loop through all pixels in the output image top-right to bottom-left.
    for (int col=m_last_column; col>=0; --col) {
      for (int row=0; row<m_parent_ptr->m_num_output_rows; ++row) {

        //printf("Accum pass 1 col = %d, row = %d\n", col, row);
        int num_disp = get_num_disp(col, row);
        if (num_disp == 0)
          continue;
        CostType * const local_cost_ptr = m_parent_ptr->get_cost_vector(col, row);
        bool debug = false;//((row == 244) && (col == 341));

        // Top
        output_accum_ptr = m_buffer_ptr->get_output_accum_ptr(MultiAccumRowBuffer::PASS_ONE);
        if ((row > 0) && (col < m_last_column)) {
          int pixel_diff = m_parent_ptr->get_path_pixel_diff(*m_image_ptr, col, row, 0, -1);
          AccumCostType* const prior_accum_ptr = m_buffer_ptr->get_trailing_pixel_accum_ptr(0, -1, MultiAccumRowBuffer::PASS_ONE);
          m_parent_ptr->evaluate_path(col, row, col, row-1,
                                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                                       pixel_diff, debug);

          AccumCostType* const prior_accum_ptr2 = m_buffer_ptr->get_trailing_pixel_accum_ptr(1, 0, MultiAccumRowBuffer::PASS_ONE);
          m_parent_ptr->evaluate_path(col, row, col+1, row,
                                       prior_accum_ptr2, full_prior_ptr, local_cost_ptr, m_temp_buffer.get(), 
                                       pixel_diff, debug);
          for (int d=0; d<num_disp; ++d)
            output_accum_ptr[d] = (output_accum_ptr[d] + m_temp_buffer[d])/2;
        }
        else // Just init to the local cost
          for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

        m_buffer_ptr->next_pixel();
      } // End col loop

      m_buffer_ptr->next_row(col==0);
    } // End row loop
  } // End task_T

  void task_TR() {
    AccumCostType* full_prior_ptr = m_full_prior_buffer.get();
    AccumCostType* output_accum_ptr;

    // Loop through all pixels in the output image top-right to bottom-left.
    for (int col=m_last_column; col>=0; --col) {
      for (int row=0; row<m_parent_ptr->m_num_output_rows; ++row) {

        //printf("Accum pass 1 col = %d, row = %d\n", col, row);
        int num_disp = get_num_disp(col, row);
        if (num_disp == 0)
          continue;
        CostType * const local_cost_ptr = m_parent_ptr->get_cost_vector(col, row);
        bool debug = false;//((row == 244) && (col == 341));

        // Top right
        output_accum_ptr = m_buffer_ptr->get_output_accum_ptr(MultiAccumRowBuffer::PASS_ONE);
        if ((row > 0) && (row < m_last_row) && (col < m_last_column)) {
          int pixel_diff = m_parent_ptr->get_path_pixel_diff(*m_image_ptr, col, row, 1, -1);
          AccumCostType* const prior_accum_ptr = m_buffer_ptr->get_trailing_pixel_accum_ptr(1, -1, MultiAccumRowBuffer::PASS_ONE);
          m_parent_ptr->evaluate_path(col, row, col+1, row-1,
                                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                                       pixel_diff, debug);

          AccumCostType* const prior_accum_ptr2 = m_buffer_ptr->get_trailing_pixel_accum_ptr(1, 1, MultiAccumRowBuffer::PASS_ONE);
          m_parent_ptr->evaluate_path(col, row, col+1, row+1,
                                       prior_accum_ptr2, full_prior_ptr, local_cost_ptr, m_temp_buffer.get(), 
                                       pixel_diff, debug);
          for (int d=0; d<num_disp; ++d)
            output_accum_ptr[d] = (output_accum_ptr[d] + m_temp_buffer[d])/2;
        }
        else // Just init to the local cost
          for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

        m_buffer_ptr->next_pixel();
      } // End col loop

      m_buffer_ptr->next_row(col==0);
    } // End row loop
  } // End task_TR

}; // End class SmoothPathAccumTask



} // End namespace stereo
} // End namespace vw


#endif // __SGM_ASSIST_H_
