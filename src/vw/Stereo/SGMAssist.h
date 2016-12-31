
#ifndef __SGM_ASSIST_H__
#define __SGM_ASSIST_H__

#include <vw/Stereo/SGM.h>
#include <vw/Core/Thread.h>
#include <vw/Core/ThreadPool.h>
#include <vw/Image/PixelIterator.h>

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

    m_fit_params = Matrix<double,6,9>( pinvA_data );
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
    double denom = 4.0 * vals[0] * vals[1] - ( vals[2] * vals[2] ); // = 4ab - c^2
    if (fabs(denom) < 0.01) {
      //std::cout << "DENOM DEBUG!!!\n";
      //std::cout << "vals  = " << vals << std::endl;
      //std::cout << "denom = " << denom << std::endl;
      return false;
    }
    Vector2f offset( ( vals[2] * vals[4] - 2.0 * vals[1] * vals[3] ) / denom,   // ce - 2bd
                     ( vals[2] * vals[3] - 2.0 * vals[0] * vals[4] ) / denom ); // cd - 2ae
    /*
    if ( norm_2(offset) > 0.99 ) {
      //std::cout << "DEBUG!!!\n";
      //std::cout << "offset  = " << offset << std::endl;
      //std::cout << "vals  = " << vals << std::endl;
      //std::cout << "denom = " << denom << std::endl;
      //vw_throw( NoImplErr() << "DEBUG!\n" );
      return false;
    }
    dx = offset[0];
    dy = offset[1];*/
    
    
    dx = offset[0];
    dy = offset[1];
        
    // APPLY CORRECTION
    double sX = 0.34574;
    double sY = 0.38944;
    dx = erf(dx/(sX*sqrt(2.0))) / 2.0;
    dy = erf(dy/(sY*sqrt(2.0))) / 2.0;
    
    //std::cout << offset << " -F-> " << dx << ", " << dy << std::endl;
    
    if ( norm_2(Vector2(dx, dy)) >= 0.5 ) {
      double scale = norm_2(Vector2(dx, dy)) / 0.5;
      dx /= scale;
      dy /= scale;
    }
    
    
    
// TODO: Check afterwards.
    //if ( norm_2(offset) > 0.99 ) {
    if ( norm_2(Vector2(dx, dy)) >= 1.0 ) {
      std::cout << "DEBUG!!!\n";
      std::cout << "offset  = " << offset << std::endl;
      std::cout << "vals  = " << vals << std::endl;
      std::cout << "denom = " << denom << std::endl;
      //vw_throw( NoImplErr() << "DEBUG!\n" );
      return false;
    }
    
    
    return true;
  } // End find_peak()

private: // Variables

  Matrix<float,6,9> m_fit_params;

}; // End class ParabolaFit2d



//==================================================================================



/// Class to make counting up "on" mask pixel in a box more efficient.
/// - This class performs a convolution of two images where it
///   computes the a percentage of nonzero pixels in the window.
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
    //std::cout << "Pop " << m_column_sums.front() << ", sum = " << m_curr_sum << std::endl;
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
    //std::cout << "Push " << new_col_sum << ", sum = " << m_curr_sum << std::endl;
    
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
  
  double recompute() {
    m_column_sums = std::queue<int>();
    //std::cout << "Cleared sum\n";
    
    for (int col=m_curr_pos[0]; col<m_curr_pos[0]+m_box_size[0]; ++col) {
      int col_sum = 0;
      for (int row=m_curr_pos[1]; row<m_curr_pos[1]+m_box_size[1]; ++row) {
        if (m_right_image_mask->operator()(col,row) > 0){
          col_sum += 1;
        }
      }
      m_column_sums.push(col_sum);
      m_curr_sum += col_sum;
      //std::cout << "Push " << col_sum << ", sum = " << m_curr_sum << std::endl;
    }
    return static_cast<double>(m_curr_sum) / m_box_area;
  }
  
}; // End class IterativeMaskBoxCounter


//==================================================================================


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

  /// Construct the buffers.
  /// - num_paths_in_pass must be four (8 directions) or eight (16 directions)
  MultiAccumRowBuffer(const SemiGlobalMatcher* parent_ptr,
                      const int  num_paths_in_pass=4,
                      const bool vertical=false) {
    m_parent_ptr        = parent_ptr;  
    m_num_paths_in_pass = num_paths_in_pass;
    m_vertical          = vertical;
    
    const int num_cols = m_parent_ptr->m_num_output_cols;
    const int num_rows = m_parent_ptr->m_num_output_rows;
    
    m_line_size = num_cols;
    if (vertical)
      m_line_size = num_rows;

    //std::cout << "m_line_size = " << m_line_size << std::endl;

    // Instantiate two single-row buffers that will be used to temporarily store
    //  accumulated cost info until it is no longer needed.
    // - Within each buffer, data is indexed in order [col][pass][disparity]
    // - The actual data size in the buffer will vary each line, so it is 
    //    initialized to be the maximum possible size.
    const size_t buffer_pixel_size = m_num_paths_in_pass*m_parent_ptr->m_num_disp;
    m_buffer_size       = m_line_size*buffer_pixel_size;
    m_buffer_size_bytes = m_buffer_size*sizeof(SemiGlobalMatcher::AccumCostType);
   
    //std::cout << "Allocating buffer: " << m_buffer_size_bytes << std::endl;
   
    // Allocate buffers that store accumulation scores
    m_bufferA.reset(new SemiGlobalMatcher::AccumCostType[m_buffer_size]);
    m_bufferB.reset(new SemiGlobalMatcher::AccumCostType[m_buffer_size]);
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
    memset(m_lead_buffer, 0, m_buffer_size_bytes); // Init this buffer to zero
    //std::fill(m_lead_buffer, m_lead_buffer+m_buffer_size, m_parent_ptr->get_bad_accum_val());
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
      for (int i=0; i<m_line_size; ++i)
        m_offsets_lead[i] = (raw_offsets[new_lead_index+i] - start_offset) * m_num_paths_in_pass;
    } else { // vertical
      // In the vertical case we need to rebuild a set of offsets to describe the column.
      size_t position = 0;
      for (int i=0; i<m_line_size; ++i) {
        m_offsets_lead[i] = position;
        size_t size = m_parent_ptr->get_num_disparities(m_current_col, i) * m_num_paths_in_pass;
        position += size;
        //std::cout << "offset lead " << i << " = " << position << std::endl;
      }
    }
  }

  /// Add the results in the leading buffer to the main class accumulation buffer.
  /// - The scores from each pass are added.
  void add_lead_buffer_to_accum() {
    size_t buffer_index = 0;
    SemiGlobalMatcher::AccumCostType* out_ptr = m_parent_ptr->m_accum_buffer.get();
    if (!m_vertical) { // horizontal
      for (int col=0; col<m_parent_ptr->m_num_output_cols; ++col) {
        int num_disps = m_parent_ptr->get_num_disparities(col, m_current_row);
        for (int pass=0; pass<m_num_paths_in_pass; ++pass) {
          size_t out_index = m_parent_ptr->m_buffer_starts(col, m_current_row);
          for (int d=0; d<num_disps; ++d) {
            out_ptr[out_index++] += m_lead_buffer[buffer_index++];
            //printf("row, col, pass, d = %d, %d, %d, %d ->> %d ->> %d\n", 
            //    m_current_row, col, pass, d, m_trail_buffer[buffer_index], m_parent_ptr->m_accum_buffer[out_index]);
          } // end disp loop
        } // end pass loop
      } // end col loop
    } else { // vertical
      for (int row=0; row<m_parent_ptr->m_num_output_rows; ++row) {
        int num_disps = m_parent_ptr->get_num_disparities(m_current_col, row);
        for (int pass=0; pass<m_num_paths_in_pass; ++pass) {
          size_t out_index = m_parent_ptr->m_buffer_starts(m_current_col, row);
          for (int d=0; d<num_disps; ++d) {
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
    //std::fill(m_lead_buffer, m_lead_buffer+m_buffer_size, m_parent_ptr->get_bad_accum_val());
    memset(m_lead_buffer, 0, m_buffer_size_bytes); // Init this buffer to zero
    
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
    memset(m_lead_buffer, 0, m_buffer_size_bytes);    
    fill_lead_offset_buffer();
  }
  
  /// Get the pointer to write the output of the current pass to
  SemiGlobalMatcher::AccumCostType * get_output_accum_ptr(PassIndex pass) {
    int    num_disps   = m_parent_ptr->get_num_disparities(m_current_col, m_current_row);
    size_t pass_offset = num_disps*pass;
    size_t offset      = (m_offsets_lead[m_current_col] + pass_offset);
    if (m_vertical)
      offset = (m_offsets_lead[m_current_row] + pass_offset);
    //std::cout << "output accum offset = " <<offset-pass_offset
    //          << " + pass " << pass_offset <<", num_disps = " << num_disps << std::endl;
    return m_lead_buffer + offset;
  }
  
  /// Gets the pointer to the accumulation buffer for the indicated pixel/pass
  SemiGlobalMatcher::AccumCostType * get_trailing_pixel_accum_ptr(int col_offset, int row_offset, PassIndex pass) {
 
    int    col         = m_current_col + col_offset;
    int    row         = m_current_row + row_offset;
    int    num_disps   = m_parent_ptr->get_num_disparities(col, row);
    size_t pass_offset = num_disps*pass;

    // Just get the index of the column in the correct buffer, then add an offset for the selected pass.
    
    if (!m_vertical) { // horizontal
      if (row_offset == 0) { 
        // Same row, must be the in the leading buffer
        return m_lead_buffer + (m_offsets_lead[col] + pass_offset);
      } else { 
        // Different row, must be in the trailing buffer
        return m_trail_buffer + (m_offsets_trail[col] + pass_offset);
      }
    } else { // vertical
      if (col_offset == 0) { 
        // Same row, must be the in the leading buffer
        //std::cout << "offset = " << pass_offset << ", lead loc = " << m_offsets_lead[row] << std::endl;
        return m_lead_buffer + (m_offsets_lead[row] + pass_offset);
      } else { 
        // Different row, must be in the trailing buffer
        //std::cout << "offset = " << pass_offset << ", trail loc = " << m_offsets_trail[row] << std::endl;
        return m_trail_buffer + (m_offsets_trail[row] + pass_offset);
      }
    }
  } // End function get_trailing_pixel_accum_ptr

private:

  const SemiGlobalMatcher* m_parent_ptr; ///< Need a handle to the parent SGM object

  bool m_vertical; ///< Raster orientation
  int  m_num_paths_in_pass; ///< Must be 4 or 8

  int    m_line_size; ///< Length of a column(horizontal) or a row(vertical) in pixels.
  size_t m_buffer_size, m_buffer_size_bytes;
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



/**
  A single line, single pass SGM accumulation buffer designed to be used by a single thread.
  - Each line-pass instance of the SGM accumulation problem is independent until the end.
  - Instances of this class must share access to the parent class accumulation buffer so that
    only one instance of the class touches it at any one time.
  - This class is much less complicated than the multi-line buffer class!
*/
class OneLineBuffer{
public:
  /// Default constructor does not initialize the object
  OneLineBuffer() {}

  /// Initializing constructor
  OneLineBuffer(SemiGlobalMatcher const* parent_ptr) {
    initialize(parent_ptr);
  }

  /// Initialize the buffer.
  /// - All this really does is allocate a memory buffer of the maximum possible required size.
  void initialize(const SemiGlobalMatcher* parent_ptr) {

    // Figure out the max possible line length (diagonal line down the center)    
    const int num_cols = parent_ptr->m_num_output_cols;
    const int num_rows = parent_ptr->m_num_output_rows;
    int line_size = sqrt(num_cols*num_cols + num_rows*num_rows) + 1;
    //std::cout << "m_line_size = " << m_line_size << std::endl;

    // Instantiate single-row buffer that will be used to temporarily store
    //  accumulated cost info until it can be addded to the main SGM class buffer.
    // - Within each buffer, data is indexed in order [pixel][disparity]
    // - The actual data size in the buffer will vary each line, so it is 
    //    initialized to be the maximum possible size.
    m_num_disp          = parent_ptr->m_num_disp;
    m_buffer_size       = line_size*m_num_disp;
    m_buffer_size_bytes = m_buffer_size*sizeof(SemiGlobalMatcher::AccumCostType);
   
    std::cout << "OneLineBuffer - allocating buffer size: " << m_buffer_size_bytes << std::endl;
   
    // Allocate the accumulation buffer
    m_buffer.reset(new SemiGlobalMatcher::AccumCostType[m_buffer_size]);
    
    // Set up the small buffer
    m_bad_disp_value = parent_ptr->get_bad_accum_val();
    m_full_prior_buffer.reset(new SemiGlobalMatcher::AccumCostType[m_num_disp]);
  }

  /// Clear both buffers
  void clear_buffers() {
    memset(m_buffer.get(), 0, m_buffer_size_bytes);
    
    for (size_t i=0; i<m_num_disp; ++i)
      m_full_prior_buffer[i] = m_bad_disp_value;
  }

  /// Get the pointer to the start of the output accumulation buffer
  SemiGlobalMatcher::AccumCostType * get_output_accum_ptr() {
    return m_buffer.get();
  }

  // Return the smaller buffer
  SemiGlobalMatcher::AccumCostType * get_full_prior_ptr() {
    return m_full_prior_buffer.get();
  }
  
private: // Variables

  SemiGlobalMatcher::AccumCostType m_bad_disp_value;
  size_t m_buffer_size, m_buffer_size_bytes, m_num_disp;
  
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
    vw_throw( LogicErr() << "Error: Should always be a free buffer!\n" );
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
class PixelPassTask : public Task {
public:

  /// Constructor, pass in all the needed links.
  PixelPassTask(ImageView<uint8>  const* image_ptr,
                SemiGlobalMatcher      * parent_ptr,
                OneLineBufferManager   * buffer_manager_ptr,
                PixelLineIterator pixel_loc_iter)
    : m_image_ptr(image_ptr), m_parent_ptr(parent_ptr), m_buffer_manager_ptr(buffer_manager_ptr),
      m_pixel_loc_iter(pixel_loc_iter) {
  }

  /// Do the work!
  virtual void operator()() {

    // Retrive a memory buffer to work with
    size_t buffer_id = m_buffer_manager_ptr->get_free_buffer_id();
    OneLineBuffer                   * buff_ptr       = m_buffer_manager_ptr->get_line_buffer(buffer_id);
    SemiGlobalMatcher::AccumCostType* full_prior_ptr = buff_ptr->get_full_prior_ptr();

    // Make a copy of the pixel iterator so we can re-use it for accum buffer addition
    PixelLineIterator pixel_loc_iter_copy(m_pixel_loc_iter);

    //std::cout << "Starting task with buffer id " << buffer_id 
    //          << " and location " << pixel_loc_iter_copy.to_string() << std::endl;

    //const bool debug = false;
    int last_pixel_val = -1;
    int col_prev = -1, row_prev = -1; // Previous row and column
    
    // Get the start of the output accumulation buffer
    // - Storage here is simply num_disps for each pixel in the line, one after the other.
    SemiGlobalMatcher::AccumCostType* computed_accum_ptr = buff_ptr->get_output_accum_ptr();
    SemiGlobalMatcher::AccumCostType* prior_accum_ptr    = 0;

    while (pixel_loc_iter_copy.is_good()) {

      // Get current location
      const int col = pixel_loc_iter_copy.col();
      const int row = pixel_loc_iter_copy.row();

      bool debug = false;//col==0;

      // Get information about the current pixel location from the parent
      int num_disp = m_parent_ptr->get_num_disparities(col, row);
      SemiGlobalMatcher::CostType * const local_cost_ptr = m_parent_ptr->get_cost_vector(col, row);

      // Fill in the accumulated value in the bottom buffer
      int curr_pixel_val = static_cast<int>(m_image_ptr->operator()(col, row));
      int pixel_diff     = abs(curr_pixel_val - last_pixel_val);
      
      //printf("Loc %d, %d, diff = %d, num_disp = %d\n", col, row, pixel_diff, num_disp);
      //std::cout << "DEBUG " << m_parent_ptr->m_disp_bound_image(col,row) << std::endl;
      
      if (last_pixel_val >= 0) { // All pixels after the first
        m_parent_ptr->evaluate_path( col, row, col_prev, row_prev,
                                    prior_accum_ptr, full_prior_ptr, local_cost_ptr, computed_accum_ptr, 
                                    pixel_diff, debug );
      } else { // First pixel only, nothing to accumulate.
        for (int d=0; d<num_disp; ++d) 
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

  } // End operator() function

  /// Add the computed buffer results to the parent accumulation buffer
  void update_accum_buffer(OneLineBuffer * buff_ptr) {

    // Get the start of the output accumulation buffer
    // - Storage here is simply num_disps for each pixel in the line, one after the other.
    SemiGlobalMatcher::AccumCostType* computed_accum_ptr = buff_ptr->get_output_accum_ptr();

    // Loop through all pixels in the line
    while (m_pixel_loc_iter.is_good()) {

      // Get current location
      const int col = m_pixel_loc_iter.col();
      const int row = m_pixel_loc_iter.row();

      // Get information about the current pixel location from the parent
      int num_disp = m_parent_ptr->get_num_disparities(col, row);
      
      //printf("Loc %d, %d, num_disp = %d\n", col, row, num_disp);
      
      SemiGlobalMatcher::AccumCostType* output_accum_ptr = m_parent_ptr->get_accum_vector(col, row);

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

  ImageView<uint8>  const* m_image_ptr;
  SemiGlobalMatcher      * m_parent_ptr;
  OneLineBufferManager   * m_buffer_manager_ptr;
  PixelLineIterator m_pixel_loc_iter; ///< Keeps track of the pixel position

}; // End class OneLineBuffer




} // End namespace stereo
} // End namespace vw


#endif // __SGM_ASSIST_H_
