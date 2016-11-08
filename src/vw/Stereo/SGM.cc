
#include <queue>
#include <vw/Stereo/SGM.h>
#include <vw/Core/Debugging.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/PixelMask.h>
#include <vw/Cartography/GeoReferenceUtils.h>

#if defined(VW_ENABLE_SSE) && (VW_ENABLE_SSE==1)
  #include <emmintrin.h>
  #include <smmintrin.h> // SSE4.1
#endif

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
  
};


//=========================================================================

void SemiGlobalMatcher::set_parameters(CostFunctionType cost_type,
                                       bool use_mgm,
                                       int min_disp_x, int min_disp_y,
                                       int max_disp_x, int max_disp_y,
                                       int kernel_size, uint16 p1, uint16 p2,
                                       int ternary_census_threshold) {
  m_cost_type   = cost_type;
  m_use_mgm     = use_mgm;
  m_min_disp_x  = min_disp_x;
  m_min_disp_y  = min_disp_y;
  m_max_disp_x  = max_disp_x;
  m_max_disp_y  = max_disp_y;
  m_kernel_size = kernel_size;
  m_ternary_census_threshold = ternary_census_threshold;
  
  m_num_disp_x = m_max_disp_x - m_min_disp_x + 1;
  m_num_disp_y = m_max_disp_y - m_min_disp_y + 1;
  size_t size_check = m_num_disp_x * m_num_disp_y;;
  if (size_check > (size_t)std::numeric_limits<DisparityType>::max())
    vw_throw( NoImplErr() << "Number of disparities is too large for data type!\n" );
  m_num_disp   = m_num_disp_x * m_num_disp_y;   

  if (p1 > 0) // User provided
    m_p1 = p1; 
  else { // Choose based on the kernel size and cost type
    switch(cost_type) {
    case CENSUS_TRANSFORM:
      switch(kernel_size){
      case 3:  m_p1 = 3;  break; // 0 - 8
      case 5:  m_p1 = 15; break; // 0 - 24
      case 7:  m_p1 = 30; break; // 0 - 48
      default: m_p1 = 3; break; // Unsupported!
      };
      break;
    case TERNARY_CENSUS_TRANSFORM:
      switch(kernel_size){
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
      switch(kernel_size){
      case 3:  m_p2 = 70;   break; // 0 - 8
      case 5:  m_p2 = 750;  break; // 0 - 24
      case 7:  m_p2 = 1500; break; // 0 - 48
      default: m_p2 = 22;   break; // Unsupported!
      };
      break;
    case TERNARY_CENSUS_TRANSFORM:
      switch(kernel_size){
      case 5:  m_p2 = 1500; break; // 0 - 48
      case 7:  m_p2 = 2000; break; // 0 - 64
      case 9:  m_p2 = 2000; break; // 0 - 64
      default: m_p2 = 30;   break; // Unsupported!
      };
      break;
    default: // MAD
      m_p1 = 250; // 0 - 255 range
      break;
    };
  } // End p2 cases
  
} // End set_parameters


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
                                                  DisparityImage const* prev_disparity,
                                                  int search_buffer) {

  //Timer timer_total("Populate disparity bounds");

  vw_out(VerboseDebugMessage, "stereo") << "disparity bound image size = " << bounding_box(m_disp_bound_image) << std::endl;

  // The masks are assumed to be the same size as the output image.
  // TODO: Check or automatically compute the left valid mask size!
  bool left_mask_valid = false, right_mask_valid = false;
  if (left_image_mask) {
    vw_out(VerboseDebugMessage, "stereo") << "Left  mask image size:" << bounding_box(*left_image_mask ) << std::endl;
    
    if ( (left_image_mask->cols() == m_disp_bound_image.cols()) && 
         (left_image_mask->rows() == m_disp_bound_image.rows())   )
      left_mask_valid = true;
    else
      vw_throw( LogicErr() << "Left mask size does not match the output size!\n" );
  }
  if (right_image_mask) {
    vw_out(VerboseDebugMessage, "stereo") << "Right mask image size:" << bounding_box(*right_image_mask) << std::endl;
    
    if ( (right_image_mask->cols() >= m_disp_bound_image.cols()+m_num_disp_x-1) && 
         (right_image_mask->rows() >= m_disp_bound_image.rows()+m_num_disp_y-1)   )
      right_mask_valid = true;
    else
      vw_throw( LogicErr() << "Right mask size does not match the output size!\n" );
  }

  // The low-res disparity image must be half-resolution.
  const int SCALE_UP = 2; 
  
  // Require that the right image mask is valid for this percentage of the
  // search range for each pixel in the left image.
  // - If this value is too low, many border pixels will be assigned the full
  //   search range and significantly slow down SGM!
  const double MIN_MASK_OVERLAP = 0.95;

  // There needs to be some "room" in the disparity search space for
  // us to discard prior results on the edge as not trustworthy predictors.
  // In other words, don't mark any pixels as edge if the one dimension has search space 1!
  const bool check_x_edge = ((m_max_disp_x - m_min_disp_x) > 1);
  const bool check_y_edge = ((m_max_disp_y - m_min_disp_y) > 1);

  double area = 0, percent_trusted = 0, percent_masked = 0;

  // This class will check the right image mask in an efficient manner.
  IterativeMaskBoxCounter right_mask_checker(right_image_mask, Vector2i(m_num_disp_x, m_num_disp_y));

  // Loop through the output disparity image and compute a search range for each pixel
  int r_in, c_in;
  int dx_scaled, dy_scaled;
  PixelMask<Vector2i> input_disp;
  Vector4i bounds;
  for (int r=0; r<m_disp_bound_image.rows(); ++r) {
    r_in = r / SCALE_UP;
    for (int c=0; c<m_disp_bound_image.cols(); ++c) {

      // TODO: This will fail if not used with positive search ranges!!!!!!!!!!!!!!!!!!!!!!!
      // Verify that there is sufficient overlap with the right image mask
      if (right_mask_valid) {
        double right_percent  = right_mask_checker.next_pixel();
        bool   right_check_ok = (right_percent >= MIN_MASK_OVERLAP);
        // If none of the right mask pixels were valid, flag this pixel as invalid.
        if (!right_check_ok) {
          m_disp_bound_image(c,r) = Vector4i(0,0,-1,-1); // Zero search area.
          ++percent_masked;
          continue;
        }
      } // End right mask handling

      // If the left mask is invalid here, flag the pixel as invalid.
      // - Do this second so that our right image pixel tracker stays up to date.
      if (left_mask_valid && (left_image_mask->operator()(c,r) == 0)) {
        m_disp_bound_image(c,r) = Vector4i(0,0,-1,-1); // Zero search area.
        ++percent_masked;
        continue;
      }

      // If a previous disparity was provided, see if we have a valid disparity 
      // estimate for this pixel.
      bool good_disparity = false;
      if (prev_disparity) {
        c_in = c / SCALE_UP;
        // Verify that the pixel we want exists
        if ( (c_in >= prev_disparity->cols()) || (r_in >= prev_disparity->rows()) ) {
          vw_throw( LogicErr() << "Size error!\n" );
        }

        input_disp = prev_disparity->operator()(c_in,r_in);
        
        // Disparity values on the edge of our 2D search range are not considered trustworthy!
        dx_scaled = input_disp[0] * SCALE_UP; 
        dy_scaled = input_disp[1] * SCALE_UP;
        
        bool on_edge = (  ( check_x_edge && ((dx_scaled <= m_min_disp_x) || (dx_scaled >= m_max_disp_x)) )
                       || ( check_y_edge && ((dy_scaled <= m_min_disp_y) || (dy_scaled >= m_max_disp_y)) ) );

        good_disparity = (is_valid(input_disp) && !on_edge);
        
      } // End prev disparity check
      
      if (good_disparity) {
      
        // We are more confident in the prior disparity, search nearby.
        bounds[0]  = dx_scaled - search_buffer; // Min x
        bounds[2]  = dx_scaled + search_buffer; // Max X
        bounds[1]  = dy_scaled - search_buffer; // Min y
        bounds[3]  = dy_scaled + search_buffer; // Max y
      
        // Constrain to global limits
        if (bounds[0] < m_min_disp_x) bounds[0] = m_min_disp_x;
        if (bounds[1] < m_min_disp_y) bounds[1] = m_min_disp_y;
        if (bounds[2] > m_max_disp_x) bounds[2] = m_max_disp_x;
        if (bounds[3] > m_max_disp_y) bounds[3] = m_max_disp_y;
      
        percent_trusted += 1.0;
      } else {
        // Not a trusted prior disparity, search the entire range!
        bounds = Vector4i(m_min_disp_x, m_min_disp_y, m_max_disp_x, m_max_disp_y); // DEBUG
      }
      
      m_disp_bound_image(c,r) = bounds;
      area += (bounds[3]-bounds[1]+1)*(bounds[2]-bounds[0]+1);
    } // End col loop
    right_mask_checker.advance_row();
  } // End row loop
  
  // Compute some statistics for help improving the speed
  double num_pixels = m_disp_bound_image.rows()*m_disp_bound_image.cols();
  area            /= num_pixels;
  percent_trusted /= num_pixels;
  percent_masked  /= num_pixels;
  
  double max_search_area    = (m_max_disp_x-m_min_disp_x+1)*(m_max_disp_y-m_min_disp_y+1);
  double percent_full_range = 1.0 - (percent_trusted+percent_masked);
  vw_out(InfoMessage, "stereo") << "Max pixel search area  = "            << max_search_area    << std::endl;
  vw_out(InfoMessage, "stereo") << "Mean pixel search area = "            << area               << std::endl;
  vw_out(InfoMessage, "stereo") << "Percent trusted prior disparities = " << percent_trusted    << std::endl;
  vw_out(InfoMessage, "stereo") << "Percent masked pixels  = "            << percent_masked     << std::endl;
  vw_out(InfoMessage, "stereo") << "Percent full search range pixels  = " << percent_full_range << std::endl;
  
  // Return false if the image cannot be processed
  if ((area <= 0) || (percent_masked >= 100))
    return false;
  return true;

} // End populate_disp_bound_image


void SemiGlobalMatcher::allocate_large_buffers() {

  //Timer timer_total("Memory allocation");

  // Init the starts data storage
  m_buffer_starts.set_size(m_num_output_cols, m_num_output_rows);

  vw_out(DebugMessage, "stereo") << "SGM: Num pixels = "      << m_num_output_rows * m_num_output_cols << std::endl;
  vw_out(DebugMessage, "stereo") << "SGM: Num disparities = " << m_num_disp << std::endl;
  
  // For each pixel, record the starting offset and add the disparity search area 
  //  of this pixel to the running offset total.
  size_t total_offset = 0;
  for (int r=0; r<m_num_output_rows; ++r) {
    for (int c=0; c<m_num_output_cols; ++c) {   
      m_buffer_starts(c,r) = total_offset;
      total_offset += get_num_disparities(c, r);
    }
  }
  // Finished computing the pixel offsets.

  vw_out(DebugMessage, "stereo") << "SGM: Total disparity search area = " << total_offset << std::endl;

  // TODO: Check available memory first!

  const size_t cost_buffer_num_bytes = total_offset * sizeof(CostType);  
  vw_out(DebugMessage, "stereo") << "SGM: Allocating buffer of size: " << cost_buffer_num_bytes/(1024*1024) << " MB\n";

  m_cost_buffer.reset(new CostType[total_offset]);

  const size_t accum_buffer_num_bytes = total_offset * sizeof(AccumCostType);  
  vw_out(DebugMessage, "stereo") << "SGM: Allocating buffer of size: " << accum_buffer_num_bytes/(1024*1024) << " MB\n";

  // Allocate the requested memory and init all to zero
  m_accum_buffer.reset(new AccumCostType[total_offset]);
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
void SemiGlobalMatcher::evaluate_path( int col, int row, int col_p, int row_p,
                       AccumCostType* const prior,
                       AccumCostType*       full_prior_buffer,
                       CostType     * const local,
                       AccumCostType*       output,
                       int path_intensity_gradient, bool debug ) {

  // Decrease p2 (jump cost) with increasing disparity along the path
  AccumCostType p2_mod = m_p2;
  if (path_intensity_gradient > 0)
    p2_mod /= path_intensity_gradient;
  if (p2_mod < m_p1)
    p2_mod = m_p1;

  //int num_disparities   = get_num_disparities(col,   row  ); // Can be input arg
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
  if (debug) {
    std::cout << "Prior pixel = ("<<col_p<<","<<row_p<<")\n";
    std::cout << "Bounds  : " << pixel_disp_bounds << std::endl;
    std::cout << "Bounds_P: " << pixel_disp_bounds_p << std::endl;
  
    std::cout << "min_prior = " <<  min_prior << std::endl;
    
    std::cout << "Priors: \n";
    int i=0;    
    for (int dy=pixel_disp_bounds_p[1]; dy<=pixel_disp_bounds_p[3]; ++dy) {
      for (int dx=pixel_disp_bounds_p[0]; dx<=pixel_disp_bounds_p[2]; ++dx) {
        std::cout << prior[i] << " ";
        ++i;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Costs: \n";
    i=0;    
    for (int dy=pixel_disp_bounds[1]; dy<=pixel_disp_bounds[3]; ++dy) {
      for (int dx=pixel_disp_bounds[0]; dx<=pixel_disp_bounds[2]; ++dx) {
        std::cout << int(local[i]) << " ";
        ++i;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    
    std::cout << "Full prior buffer: \n";
    i = 0;
    for (int r=0; r<m_num_disp_y; ++r) {
      for (int c=0; c<m_num_disp_x; ++c) {
        std::cout << full_prior_buffer[i] << " ";
        ++i;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
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


      ++packed_d;
      ++full_d;
    }
  } // End loop through this disparity  
  
  if(debug) {
    int min_val   = 99999;
    int min_index = 0;
    
    std::cout << "Output: \n";
    int i=0;    
    for (int dy=pixel_disp_bounds[1]; dy<=pixel_disp_bounds[3]; ++dy) {
      for (int dx=pixel_disp_bounds[0]; dx<=pixel_disp_bounds[2]; ++dx) {
        std::cout << output[i] << " ";
        if (output[i] < min_val) {
          min_val   = output[i];
          min_index = i;
        }
        ++i;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    
    DisparityType dx, dy;
    disp_to_xy(min_index, dx, dy);
    std::cout << "Min value = " << min_val << std::endl;
    std::cout << "Disp = " << dx-4 <<", " << dy-4 << std::endl;
    std::cout << "Min index = " << min_index << std::endl<< std::endl;
    std::cout << "========================================\n\n";
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
void SemiGlobalMatcher::evaluate_path( int col, int row, int col_p, int row_p,
                       AccumCostType* const prior,
                       AccumCostType*       full_prior_buffer,
                       CostType     * const local,
                       AccumCostType*       output,
                       int path_intensity_gradient, bool debug ) {

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
  //std::cout << "pixel_disp_bounds = " << pixel_disp_bounds << std::endl;
  
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
      if (sse_index == SSE_BUFF_LEN){

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


SemiGlobalMatcher::AccumCostType 
SemiGlobalMatcher::get_accum_vector_min(int col, int row,
                                        DisparityType &dx, DisparityType &dy) {
  // Get the array
  AccumCostType const* vec = get_accum_vector(col, row);
  const int num_disp = get_num_disparities(col, row);
  
  // Get the minimum index of the array
  int min_index = 0;
  AccumCostType value = std::numeric_limits<AccumCostType>::max();
  for (int i=0; i<num_disp; ++i) {
    if (vec[i] < value) {
      value = vec[i];
      min_index = i;
    }
  }
  
  // Convert the disparity index to dx and dy
  const Vector4i bounds = m_disp_bound_image(col,row);
  int d_width  = bounds[2] - bounds[0] + 1;
  dy = (min_index / d_width);
  dx = min_index - (dy*d_width) + bounds[0];
  dy += bounds[1];
  
  //printf("%d, %d, %d -> %d, %d\n", value, min_index, d_width, dx, dy);
  
  return value;
}

SemiGlobalMatcher::DisparityImage
SemiGlobalMatcher::create_disparity_view() {
  // Init output vector
  DisparityImage disparity( m_num_output_cols, m_num_output_rows );
  // For each element in the accumulated costs matrix, 
  //  select the disparity with the lowest accumulated cost.
  //Timer timer("Calculate Disparity Minimum");
  DisparityType dx, dy;
  for ( int j = 0; j < m_num_output_rows; j++ ) {
    for ( int i = 0; i < m_num_output_cols; i++ ) {
      
      int num_disp = get_num_disparities(i, j);
      if (num_disp > 0) {
        // Valid pixel, choose the best disparity with a winner-take-all (WTA) method.
        // - Would a fancier selection method improve our results?
        get_accum_vector_min(i, j, dx, dy);
        disparity(i,j) = DisparityImage::pixel_type(dx, dy);
        
      } else { // Pixels with no search area were never valid.
        disparity(i,j) = DisparityImage::pixel_type();
        invalidate(disparity(i,j));
      }

      /*
      if ((i >= 15) && (i <= 15) && (j==40)) { // DEBUG!
        printf("ACC costs (%d,%d): %d, %d\n", i, j, dx, dy);

        AccumCostType* vec = get_accum_vector(i, j);
        const int num_disp = get_num_disparities(i, j);
        for (int k=0; k<num_disp; ++k)
          std::cout << vec[k] << " " ;
      }
      */
    }
  }
  return disparity;
}



// TODO: Clean up and move!
class HistClass{

public: // Functions
  HistClass(int num_bins, double min, double max): m_num_bins(num_bins), m_min(min), m_max(max) {
    m_data.assign(num_bins, 0);
  };
  
  void add(double val) {
    int bin = static_cast<int>(round( (m_num_bins - 1) * ( (val - m_min)/(m_max - m_min) ) ));
    if ((bin < 0) || (bin >= m_num_bins)) {
      std::cout << "Bad bin = " << bin << ", val = " << val << std::endl;
    }
    else
      m_data[bin]++;
  }
  
  double bin(int i) const {return m_data[i];}
  
  void write(std::string const& path) const {
    std::ofstream f(path.c_str());
    for (int i=0; i<m_num_bins; ++i) {
      f << m_data[i] << std::endl;
    }
    f.close();
  }
  
private: // Functions

  int    m_num_bins;
  double m_min, m_max;
  std::vector<uint64> m_data;
};

ImageView<PixelMask<Vector2f> > SemiGlobalMatcher::
create_disparity_view_subpixel(DisparityImage const& integer_disparity) {

  //Timer timer("Calculate Subpixel Disparity");

  typedef  PixelMask<Vector2f> p_type;
  ImageView<p_type> disparity(m_num_output_cols, m_num_output_rows);
  ParabolaFit2d fitter;
  
  //// DEBUG
  //HistClass hist_dx(201, -1.0, 1.0), hist_dy(201, -1.0, 1.0);
  
  // For each element in the accumulated costs matrix, 
  //  select the disparity with the lowest accumulated cost.
  double percent_bad = 0;
  double delta_x, delta_y;
  for ( int j = 0; j < m_num_output_rows; j++ ) {
    for ( int i = 0; i < m_num_output_cols; i++ ) {
      
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
      
      // Get the best disparity from the integer image.
      int dx = integer_pixel[0];
      int dy = integer_pixel[1];
      int min_index = (dy-bounds[1])*width + (dx-bounds[0]);

      // Fetch the 8 adjacent accumulation vector values

      // Don't interpolate out of bounds
      int x_left  = -1;
      int x_right =  1;
      int y_up    = -width;
      int y_down  =  width;
      if (dx == bounds[0]) x_left  = 0;
      if (dx == bounds[2]) x_right = 0;
      if (dy == bounds[1]) y_up    = 0;
      if (dy == bounds[3]) y_down  = 0;  

      // Apply subpixel correction and apply
      AccumCostType const* accum_vec = get_accum_vector(i, j);
      bool valid = 
        fitter.find_peak( accum_vec[min_index+x_left+y_up  ],  accum_vec[min_index+y_up  ], accum_vec[min_index+x_right+y_up  ],
                          accum_vec[min_index+x_left       ],  accum_vec[min_index       ], accum_vec[min_index+x_right       ],
                          accum_vec[min_index+x_left+y_down],  accum_vec[min_index+y_down], accum_vec[min_index+x_right+y_down],
                          delta_x, delta_y);
      if (valid) {
        disparity(i,j) = p_type(dx+delta_x, dy+delta_y);
        //hist_dx.add(delta_x);
        //hist_dy.add(delta_y);
      }
      else {
        disparity(i,j) = p_type(dx, dy);
        percent_bad += 1.0;
      }

      /*
      if ((i >= 15) && (i <= 15) && (j==40)) { // DEBUG!
        printf("ACC costs (%d,%d): %d, %d\n", i, j, dx, dy);

        AccumCostType* vec = get_accum_vector(i, j);
        const int num_disp = get_num_disparities(i, j);
        for (int k=0; k<num_disp; ++k)
          std::cout << vec[k] << " " ;
      }
      */
    }
  }
  percent_bad /= (double)(m_num_output_rows*m_num_output_cols);
  //std::cout << "Percent bad = " << percent_bad << std::endl;
  vw_out(DebugMessage, "stereo") << "Subpixel interpolation failure percentage: " << percent_bad << std::endl;
  //write_image( "subpixel_disp.tif", disparity );
  //hist_dx.write("delta_x.csv");
  //hist_dy.write("delta_y.csv");
  
  return disparity;
  
}


// TODO: Replace with ASP implementation?
SemiGlobalMatcher::CostType SemiGlobalMatcher::get_cost_block(ImageView<uint8> const& left_image,
               ImageView<uint8> const& right_image,
               int left_x, int left_y, int right_x, int right_y, bool debug) {
  if (m_kernel_size == 1) { // Special handling for single pixel case
    int diff = static_cast<int>(left_image (left_x,  left_y )) - 
               static_cast<int>(right_image(right_x, right_y));
    return static_cast<CostType>(abs(diff));
  }

  // Block mean of abs dists
  const int half_kernel_size = (m_kernel_size-1) / 2;
  int sum=0, diff=0;
  for (int j=-half_kernel_size; j<=half_kernel_size; ++j) {
    for (int i=-half_kernel_size; i<=half_kernel_size; ++i) {
      diff = static_cast<int>(left_image (left_x +i, left_y +j)) - 
             static_cast<int>(right_image(right_x+i, right_y+j));
      sum += abs(diff);
    }
  }
  CostType result = sum / static_cast<int>(m_kernel_size*m_kernel_size);
  //printf("sum = %d, result = %d\n", sum, result);
  return static_cast<CostType>(result);
}

void SemiGlobalMatcher::fill_costs_block(ImageView<uint8> const& left_image,
                                         ImageView<uint8> const& right_image){
  // Make sure we don't go out of bounds here due to the disparity shift and kernel.
  size_t cost_index = 0;
  for ( int r = m_min_row; r <= m_max_row; r++ ) { // For each row in left
    int output_row = r - m_min_row;
    //int input_row  = r;
    for ( int c = m_min_col; c <= m_max_col; c++ ) { // For each column in left
      int output_col = c - m_min_col;
      //int input_col  = c;
      
      Vector4i pixel_disp_bounds = m_disp_bound_image(output_col, output_row);

      // Only compute costs in the search radius for this pixel    
      for ( int dy = pixel_disp_bounds[1]; dy <= pixel_disp_bounds[3]; dy++ ) { // For each disparity
        for ( int dx = pixel_disp_bounds[0]; dx <= pixel_disp_bounds[2]; dx++ ) {          
          
          CostType cost = get_cost_block(left_image, right_image, c, r, c+dx,r+dy, false);
          m_cost_buffer[cost_index] = cost;
          ++cost_index;
        }    
      } // End disparity loops   
    } // End x loop
  }// End y loop
  
  
}


void SemiGlobalMatcher::fill_costs_census3x3(ImageView<uint8> const& left_image,
                                             ImageView<uint8> const& right_image){
  const int half_kernel = (m_kernel_size - 1) / 2;
  const int padding     = 2*half_kernel;
                                                
  // Compute the census value for each pixel.
  // - ROI handling could be fancier but this is simple and works.
  // - The 0,0 pixels in the left and right images are assumed to be aligned.
  ImageView<uint8> left_census (left_image.cols()-padding,  left_image.rows()-padding ), 
                   right_census(right_image.cols()-padding, right_image.rows()-padding);

  if (m_cost_type == CENSUS_TRANSFORM) {
    for ( int r = 0; r < left_census.rows(); r++ )
      for ( int c = 0; c < left_census.cols(); c++ )
        left_census(c,r) = get_census_value_3x3(left_image, c+half_kernel, r+half_kernel);
    for ( int r = 0; r < right_census.rows(); r++ )
      for ( int c = 0; c < right_census.cols(); c++ )
        right_census(c,r) = get_census_value_3x3(right_image, c+half_kernel, r+half_kernel);
  } else {
    vw_throw(NoImplErr() << "The ternary sensus transform not available in size 3!\n");
  } 
  get_hamming_distance_costs(left_census, right_census);
}

void SemiGlobalMatcher::fill_costs_census5x5(ImageView<uint8> const& left_image,
                                             ImageView<uint8> const& right_image){
  const int half_kernel = (m_kernel_size - 1) / 2;
  const int padding     = 2*half_kernel;
  
  // Compute the census value for each pixel.
  // - ROI handling could be fancier but this is simple and works.
  // - The 0,0 pixels in the left and right images are assumed to be aligned.
  ImageView<uint32> left_census (left_image.cols()-padding,  left_image.rows()-padding ), 
                    right_census(right_image.cols()-padding, right_image.rows()-padding);

  if (m_cost_type == CENSUS_TRANSFORM) {
    for ( int r = 0; r < left_census.rows(); r++ )
      for ( int c = 0; c < left_census.cols(); c++ )
        left_census(c,r) = get_census_value_5x5(left_image, c+half_kernel, r+half_kernel);
    for ( int r = 0; r < right_census.rows(); r++ )
      for ( int c = 0; c < right_census.cols(); c++ )
        right_census(c,r) = get_census_value_5x5(right_image, c+half_kernel, r+half_kernel);
  } else { // TERNARY_CENSUS_TRANSFORM
    for ( int r = 0; r < left_census.rows(); r++ )
      for ( int c = 0; c < left_census.cols(); c++ )
        left_census(c,r) = get_census_value_ternary_5x5(left_image, c+half_kernel, r+half_kernel, m_ternary_census_threshold);
    for ( int r = 0; r < right_census.rows(); r++ )
      for ( int c = 0; c < right_census.cols(); c++ )
        right_census(c,r) = get_census_value_ternary_5x5(right_image, c+half_kernel, r+half_kernel, m_ternary_census_threshold);
  }
  
  get_hamming_distance_costs(left_census, right_census);
}

void SemiGlobalMatcher::fill_costs_census7x7(ImageView<uint8> const& left_image,
                                             ImageView<uint8> const& right_image){
  const int half_kernel = (m_kernel_size - 1) / 2;
  const int padding     = 2*half_kernel;
  
  // Compute the census value for each pixel.
  // - ROI handling could be fancier but this is simple and works.
  // - The 0,0 pixels in the left and right images are assumed to be aligned.
  ImageView<uint64> left_census (left_image.cols()-padding,  left_image.rows()-padding ), 
                    right_census(right_image.cols()-padding, right_image.rows()-padding);
                   
  if (m_cost_type == CENSUS_TRANSFORM) {
    for ( int r = 0; r < left_census.rows(); r++ )
      for ( int c = 0; c < left_census.cols(); c++ )
        left_census(c,r) = get_census_value_7x7(left_image, c+half_kernel, r+half_kernel);
    for ( int r = 0; r < right_census.rows(); r++ )
      for ( int c = 0; c < right_census.cols(); c++ )
        right_census(c,r) = get_census_value_7x7(right_image, c+half_kernel, r+half_kernel);
  } else { // TERNARY_CENSUS_TRANSFORM
    for ( int r = 0; r < left_census.rows(); r++ )
      for ( int c = 0; c < left_census.cols(); c++ )
        left_census(c,r) = get_census_value_ternary_7x7(left_image, c+half_kernel, r+half_kernel, m_ternary_census_threshold);
    for ( int r = 0; r < right_census.rows(); r++ )
      for ( int c = 0; c < right_census.cols(); c++ )
        right_census(c,r) = get_census_value_ternary_7x7(right_image, c+half_kernel, r+half_kernel, m_ternary_census_threshold);
  }
  
  get_hamming_distance_costs(left_census, right_census);
}

void SemiGlobalMatcher::fill_costs_census9x9(ImageView<uint8> const& left_image,
                                             ImageView<uint8> const& right_image){
  const int half_kernel = (m_kernel_size - 1) / 2;
  const int padding     = 2*half_kernel;
  
  // Compute the census value for each pixel.
  // - ROI handling could be fancier but this is simple and works.
  // - The 0,0 pixels in the left and right images are assumed to be aligned.
  ImageView<uint64> left_census (left_image.cols()-padding,  left_image.rows()-padding ), 
                    right_census(right_image.cols()-padding, right_image.rows()-padding);
                   
  if (m_cost_type == CENSUS_TRANSFORM) {
    vw_throw(NoImplErr() << "The Census transform not available in size 9!\n");
  } else { // TERNARY_CENSUS_TRANSFORM
    for ( int r = 0; r < left_census.rows(); r++ )
      for ( int c = 0; c < left_census.cols(); c++ )
        left_census(c,r) = get_census_value_ternary_9x9(left_image, c+half_kernel, r+half_kernel, m_ternary_census_threshold);
    for ( int r = 0; r < right_census.rows(); r++ )
      for ( int c = 0; c < right_census.cols(); c++ )
        right_census(c,r) = get_census_value_ternary_9x9(right_image, c+half_kernel, r+half_kernel, m_ternary_census_threshold);
  }
  
  get_hamming_distance_costs(left_census, right_census);
}

void SemiGlobalMatcher::compute_disparity_costs(ImageView<uint8> const& left_image,
                                                ImageView<uint8> const& right_image) {  
  //Timer timer("\tSGM Cost Calculation");
  if ((m_cost_type == CENSUS_TRANSFORM) || (m_cost_type == TERNARY_CENSUS_TRANSFORM)) {
    switch(m_kernel_size) {
    case 3:  fill_costs_census3x3(left_image, right_image); break;
    case 5:  fill_costs_census5x5(left_image, right_image); break;
    case 7:  fill_costs_census7x7(left_image, right_image); break;
    case 9:  fill_costs_census9x9(left_image, right_image); break;
    default: vw_throw( NoImplErr() << "Census transform is only available in size 3, 5, and 7!\n" );
    };
  }
  else { // Use the default mean of diff cost function
    // Replace this with ASP's efficient existing cost functions?
    fill_costs_block(left_image, right_image);
  }
  
/*
 int debug_row = m_num_output_rows * 0.5;
 if (m_num_output_rows > debug_row) {
   std::cout << "debug row = " << debug_row << std::endl;
   ImageView<uint8> cost_image( m_num_output_cols, m_num_disp ); // TODO: Change type?
   std::fill(cost_image.data(), cost_image.data()+m_num_output_cols*m_num_disp, 255);
   for ( int i = 0; i < m_num_output_cols; i++ ) {
     CostType* buff = get_cost_vector(i, debug_row);
     Vector4i bounds = m_disp_bound_image(i,debug_row);
     int index = 0;
     for ( int dy = bounds[1]; dy < bounds[3]; dy++ ) {
       for ( int dx = bounds[0]; dx < bounds[2]; dx++ ) {
         int d = xy_to_disp(dx, dy);
         //std::cout << "cost = " << int(buff[d]) << std::endl;
         cost_image(i,d) = buff[index]; // TODO: scale!
         ++index;
       }
     }
   }
   write_image("scanline_costs_block.tif",cost_image);
   std::cout << "Done writing line dump.\n";
 }
*/

} // end compute_disparity_costs() 



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




void SemiGlobalMatcher::two_trip_path_accumulation(ImageView<uint8> const& left_image) {

  //Timer timer_total("\tSGM Cost Propagation");

  /// Create an object to manage the temporary accumulation buffers that need to be used here.
  MultiAccumRowBuffer buff_manager(this);
  
  // Init this buffer to bad scores representing disparities that were
  //  not in the search range for the given pixel.
  boost::shared_array<AccumCostType> full_prior_buffer;
  full_prior_buffer.reset(new AccumCostType[m_num_disp]);
  for (int i=0; i<m_num_disp; ++i)
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
      bool debug = false;//((row == 244) && (col == 341));
      
      // Top left
      output_accum_ptr = buff_manager.get_output_accum_ptr(MultiAccumRowBuffer::TOP_LEFT);
      if ((row > 0) && (col > 0)) {
        // Fill in the accumulated value in the bottom buffer
        int pixel_diff = get_path_pixel_diff(left_image, col, row, 1, 1);
        AccumCostType* const prior_accum_ptr = buff_manager.get_trailing_pixel_accum_ptr(-1, -1, MultiAccumRowBuffer::TOP_LEFT);
        evaluate_path( col, row, col-1, row-1,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );
      }
      else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

      // Top
      output_accum_ptr = buff_manager.get_output_accum_ptr(MultiAccumRowBuffer::TOP);
      if (row > 0) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, 0, 1);
        AccumCostType* const prior_accum_ptr = buff_manager.get_trailing_pixel_accum_ptr(0, -1, MultiAccumRowBuffer::TOP);
        evaluate_path( col, row, col, row-1,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );
      }
      else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];
      
      // Top right
      output_accum_ptr = buff_manager.get_output_accum_ptr(MultiAccumRowBuffer::TOP_RIGHT);
      if ((row > 0) && (col < last_column)) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, -1, 1);
        AccumCostType* const prior_accum_ptr = buff_manager.get_trailing_pixel_accum_ptr(1, -1, MultiAccumRowBuffer::TOP_RIGHT);
        evaluate_path( col, row, col+1, row-1,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );
      }
      else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];
      
      // Left
      output_accum_ptr = buff_manager.get_output_accum_ptr(MultiAccumRowBuffer::LEFT);
      if (col > 0) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, 1, 0);
        AccumCostType* const prior_accum_ptr = buff_manager.get_trailing_pixel_accum_ptr(-1, 0, MultiAccumRowBuffer::LEFT);
        evaluate_path( col, row, col-1, row,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );
      }
      else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

      buff_manager.next_pixel();
    } // End col loop
    
    buff_manager.next_row(row==m_num_output_rows-1);
  } // End row loop
  
  // Done with the first trip!
  buff_manager.switch_trips();

  // Loop through all pixels in the output image for the first trip, bottom-right to top-left.
  for (int row = last_row; row >= 0; --row) {
    for (int col = last_column; col >= 0; --col) {
    
      int num_disp = get_num_disparities(col, row);
      CostType * const local_cost_ptr = get_cost_vector(col, row);
      bool debug = false;//((row == 244) && (col == 341));
              
      // Bottom right
      output_accum_ptr = buff_manager.get_output_accum_ptr(MultiAccumRowBuffer::BOT_RIGHT);
      if ((row < last_row) && (col < last_column)) {
        // Fill in the accumulated value in the bottom buffer
        int pixel_diff = get_path_pixel_diff(left_image, col, row, -1, -1);
        AccumCostType* const prior_accum_ptr = buff_manager.get_trailing_pixel_accum_ptr(1, 1, MultiAccumRowBuffer::BOT_RIGHT);
        evaluate_path( col, row, col+1, row+1,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );
      }
      else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];
      
      // Bottom
      output_accum_ptr = buff_manager.get_output_accum_ptr(MultiAccumRowBuffer::BOT);
      if (row < last_row) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, 0, -1);
        AccumCostType* const prior_accum_ptr = buff_manager.get_trailing_pixel_accum_ptr(0, 1, MultiAccumRowBuffer::BOT);
        evaluate_path( col, row, col, row+1,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );
      }
      else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];
      
      // Bottom left
      output_accum_ptr = buff_manager.get_output_accum_ptr(MultiAccumRowBuffer::BOT_LEFT);
      if ((row < last_row) && (col > 0)) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, 1, -1);
        AccumCostType* const prior_accum_ptr = buff_manager.get_trailing_pixel_accum_ptr(-1, 1, MultiAccumRowBuffer::BOT_LEFT);
        evaluate_path( col, row, col-1, row+1,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );
      }
      else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];
      
      // Right
      output_accum_ptr = buff_manager.get_output_accum_ptr(MultiAccumRowBuffer::RIGHT);
      if (col < last_column) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, -1, 0);
        AccumCostType* const prior_accum_ptr = buff_manager.get_trailing_pixel_accum_ptr(1, 0, MultiAccumRowBuffer::RIGHT);
        evaluate_path( col, row, col+1, row,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );
      }
      else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

      buff_manager.next_pixel();
    } // End col loop
    
    buff_manager.next_row(row==0);    
  } // End row loop

  // Done with both trips!
}





// This version of the function requires four passes and is based on the paper:
// MGM: A Significantly More Global Matching for Stereovision
void SemiGlobalMatcher::smooth_path_accumulation(ImageView<uint8> const& left_image) {

  //Timer timer_total("\tSGM Cost Propagation");

  int paths_per_pass = 2;

  // Create objects to manage the temporary accumulation buffers that need to be used here.
  // - Two copies are needed here, one for the two horizontal passes and one for the two vertical passes
  MultiAccumRowBuffer buff_manager_horizontal(this, paths_per_pass, false);
  MultiAccumRowBuffer buff_manager_vertical  (this, paths_per_pass, true);
  
  // Init this buffer to bad scores representing disparities that were
  //  not in the search range for the given pixel.
  boost::shared_array<AccumCostType> full_prior_buffer;
  full_prior_buffer.reset(new AccumCostType[m_num_disp]);
  for (int i=0; i<m_num_disp; ++i)
    full_prior_buffer[i] = get_bad_accum_val();  

  AccumCostType* full_prior_ptr = full_prior_buffer.get();
  AccumCostType* output_accum_ptr;
  const int last_column = m_num_output_cols - 1;
  const int last_row    = m_num_output_rows - 1;

  std::cout << "Starting first trip...\n";

  // Loop through all pixels in the output image for the first trip, top-left to bottom-right.
  for (int row=0; row<m_num_output_rows; ++row) {
    for (int col=0; col<m_num_output_cols; ++col) {
    
      //printf("Accum pass 1 col = %d, row = %d\n", col, row);
    
      int num_disp = get_num_disparities(col, row);
      if (num_disp == 0) {
        buff_manager_horizontal.next_pixel();
        continue;
      }
      CostType * const local_cost_ptr = get_cost_vector(col, row);
      bool debug = false;//((row == 244) && (col == 341));
      
      boost::shared_array<AccumCostType> temp_buffer(new AccumCostType[num_disp]);

      // Left
      output_accum_ptr = buff_manager_horizontal.get_output_accum_ptr(MultiAccumRowBuffer::PASS_ONE);
      if ((row > 0) && (col > 0)) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, -1, 0);
        AccumCostType* const prior_accum_ptr = buff_manager_horizontal.get_trailing_pixel_accum_ptr(-1, 0, MultiAccumRowBuffer::PASS_ONE);
        evaluate_path( col, row, col-1, row,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );
                       
        AccumCostType* const prior_accum_ptr2 = buff_manager_horizontal.get_trailing_pixel_accum_ptr(0, -1, MultiAccumRowBuffer::PASS_ONE);
        evaluate_path( col, row, col, row-1,
                       prior_accum_ptr2, full_prior_ptr, local_cost_ptr, temp_buffer.get(), 
                       pixel_diff, debug );
       for (int d=0; d<num_disp; ++d)
          output_accum_ptr[d] = (output_accum_ptr[d] + temp_buffer[d])/2;
      
      // Top left
      output_accum_ptr = buff_manager_horizontal.get_output_accum_ptr(MultiAccumRowBuffer::PASS_TWO);
      if ((row > 0) && (col > 0) && (col < last_column)) {

        int pixel_diff = get_path_pixel_diff(left_image, col, row, -1, -1);
        AccumCostType* const prior_accum_ptr = buff_manager_horizontal.get_trailing_pixel_accum_ptr(-1, -1, MultiAccumRowBuffer::PASS_TWO);
        evaluate_path( col, row, col-1, row-1,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );
                       
        AccumCostType* const prior_accum_ptr2 = buff_manager_horizontal.get_trailing_pixel_accum_ptr(1, -1, MultiAccumRowBuffer::PASS_TWO);
        evaluate_path( col, row, col+1, row-1,
                       prior_accum_ptr2, full_prior_ptr, local_cost_ptr, temp_buffer.get(), 
                       pixel_diff, debug );
       for (int d=0; d<num_disp; ++d)
          output_accum_ptr[d] = (output_accum_ptr[d] + temp_buffer[d])/2;
      }
      else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];
                      
      }
      else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

      buff_manager_horizontal.next_pixel();
    } // End col loop
    
    buff_manager_horizontal.next_row(row==last_row);
  } // End row loop

  
  // Done with the first trip!
  buff_manager_horizontal.switch_trips();  


  std::cout << "Starting second trip...\n";

  // Loop through all pixels in the output image for the first trip, bottom-right to top-left.
  for (int row = last_row; row >= 0; --row) {
    for (int col = last_column; col >= 0; --col) {
    
      int num_disp = get_num_disparities(col, row);
      if (num_disp == 0) {
        buff_manager_horizontal.next_pixel();
        continue;
      }
      CostType * const local_cost_ptr = get_cost_vector(col, row);
      bool debug = false;//((row == 244) && (col == 341));
      
      boost::shared_array<AccumCostType> temp_buffer(new AccumCostType[num_disp]);

      // Right
      output_accum_ptr = buff_manager_horizontal.get_output_accum_ptr(MultiAccumRowBuffer::PASS_ONE);
      if ((row < last_row) && (col < last_column)) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, 1, 0);
        AccumCostType* const prior_accum_ptr = buff_manager_horizontal.get_trailing_pixel_accum_ptr(1, 0, MultiAccumRowBuffer::PASS_ONE);
        evaluate_path( col, row, col+1, row,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );
                       
        AccumCostType* const prior_accum_ptr2 = buff_manager_horizontal.get_trailing_pixel_accum_ptr(0, 1, MultiAccumRowBuffer::PASS_ONE);
        evaluate_path( col, row, col, row+1,
                       prior_accum_ptr2, full_prior_ptr, local_cost_ptr, temp_buffer.get(), 
                       pixel_diff, debug );
       for (int d=0; d<num_disp; ++d)
          output_accum_ptr[d] = (output_accum_ptr[d] + temp_buffer[d])/2;                      
      }
      else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];

      
      // Bottom right
      output_accum_ptr = buff_manager_horizontal.get_output_accum_ptr(MultiAccumRowBuffer::PASS_TWO);
      if ((row < last_row) && (col > 0) && (col < last_column)) {

        int pixel_diff = get_path_pixel_diff(left_image, col, row, 1, 1);
        AccumCostType* const prior_accum_ptr = buff_manager_horizontal.get_trailing_pixel_accum_ptr(1, 1, MultiAccumRowBuffer::PASS_TWO);
        evaluate_path( col, row, col+1, row+1,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );

        AccumCostType* const prior_accum_ptr2 = buff_manager_horizontal.get_trailing_pixel_accum_ptr(-1, 1, MultiAccumRowBuffer::PASS_TWO);
        evaluate_path( col, row, col-1, row+1,
                       prior_accum_ptr2, full_prior_ptr, local_cost_ptr, temp_buffer.get(), 
                       pixel_diff, debug );
        for (int d=0; d<num_disp; ++d)
          output_accum_ptr[d] = (output_accum_ptr[d] + temp_buffer[d])/2;
      }
      else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];


      buff_manager_horizontal.next_pixel();
    } // End col loop
    
    buff_manager_horizontal.next_row(row==0);    
  } // End row loop

  // -- Done with the horizontal passes, switch to vertical passes.

  std::cout << "Starting third trip...\n";

  // Loop through all pixels in the output image for the first trip, bottom-left to top-right.
  for (int col = 0; col < m_num_output_cols; ++col) {
    for (int row = last_row; row >= 0; --row) {
    
      int num_disp = get_num_disparities(col, row);
      //printf("Accum pass 3 col = %d, row = %d, num_disp = %d\n", col, row, num_disp);
      if (num_disp == 0) {
        buff_manager_vertical.next_pixel();
        continue;
      }
      CostType * const local_cost_ptr = get_cost_vector(col, row);
      bool debug = false;//((row == 244) && (col == 341));
     
      boost::shared_array<AccumCostType> temp_buffer(new AccumCostType[num_disp]);
      
      // Bottom
      output_accum_ptr = buff_manager_vertical.get_output_accum_ptr(MultiAccumRowBuffer::PASS_ONE);
      if ((row < last_row) && (col > 0)) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, 0, 1);
        AccumCostType* const prior_accum_ptr = buff_manager_vertical.get_trailing_pixel_accum_ptr(0, 1, MultiAccumRowBuffer::PASS_ONE);
        evaluate_path( col, row, col, row+1,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );
                       
        AccumCostType* const prior_accum_ptr2 = buff_manager_vertical.get_trailing_pixel_accum_ptr(-1, 0, MultiAccumRowBuffer::PASS_ONE);
        evaluate_path( col, row, col-1, row,
                       prior_accum_ptr2, full_prior_ptr, local_cost_ptr, temp_buffer.get(), 
                       pixel_diff, debug );
       for (int d=0; d<num_disp; ++d)
          output_accum_ptr[d] = (output_accum_ptr[d] + temp_buffer[d])/2;
                       
      }
      else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];
      
      // Bottom left
      output_accum_ptr = buff_manager_vertical.get_output_accum_ptr(MultiAccumRowBuffer::PASS_TWO);
      if ((row > 0) && (row < last_row) && (col > 0)) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, -1, 1);
        AccumCostType* const prior_accum_ptr = buff_manager_vertical.get_trailing_pixel_accum_ptr(-1, 1, MultiAccumRowBuffer::PASS_TWO);
        evaluate_path( col, row, col-1, row+1,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );
                       
        AccumCostType* const prior_accum_ptr2 = buff_manager_vertical.get_trailing_pixel_accum_ptr(-1, -1, MultiAccumRowBuffer::PASS_TWO);
        evaluate_path( col, row, col-1, row-1,
                       prior_accum_ptr2, full_prior_ptr, local_cost_ptr, temp_buffer.get(), 
                       pixel_diff, debug );
        for (int d=0; d<num_disp; ++d)
          output_accum_ptr[d] = (output_accum_ptr[d] + temp_buffer[d])/2;
      }
      else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];
      

      buff_manager_vertical.next_pixel();
    } // End col loop
    
    buff_manager_vertical.next_row(col==last_column);    
  } // End row loop

  // Done with the first trip!
  buff_manager_vertical.switch_trips();

  std::cout << "Starting fourth trip...\n";

  // Loop through all pixels in the output image top-right to bottom-left.
  for (int col=last_column; col>=0; --col) {
    for (int row=0; row<m_num_output_rows; ++row) {
    
      //printf("Accum pass 1 col = %d, row = %d\n", col, row);
    
      int num_disp = get_num_disparities(col, row);
      if (num_disp == 0) {
        buff_manager_vertical.next_pixel();
        continue;
      }
      CostType * const local_cost_ptr = get_cost_vector(col, row);
      bool debug = false;//((row == 244) && (col == 341));
      
      boost::shared_array<AccumCostType> temp_buffer(new AccumCostType[num_disp]);
      
      // Top
      output_accum_ptr = buff_manager_vertical.get_output_accum_ptr(MultiAccumRowBuffer::PASS_ONE);
      if ((row > 0) && (col < last_column)) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, 0, -1);
        AccumCostType* const prior_accum_ptr = buff_manager_vertical.get_trailing_pixel_accum_ptr(0, -1, MultiAccumRowBuffer::PASS_ONE);
        evaluate_path( col, row, col, row-1,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );
                       
        AccumCostType* const prior_accum_ptr2 = buff_manager_vertical.get_trailing_pixel_accum_ptr(1, 0, MultiAccumRowBuffer::PASS_ONE);
        evaluate_path( col, row, col+1, row,
                       prior_accum_ptr2, full_prior_ptr, local_cost_ptr, temp_buffer.get(), 
                       pixel_diff, debug );
       for (int d=0; d<num_disp; ++d)
          output_accum_ptr[d] = (output_accum_ptr[d] + temp_buffer[d])/2;
      }
      else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];
      
      // Top right
      output_accum_ptr = buff_manager_vertical.get_output_accum_ptr(MultiAccumRowBuffer::PASS_TWO);
      if ((row > 0) && (row < last_row) && (col < last_column)) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, 1, -1);
        AccumCostType* const prior_accum_ptr = buff_manager_vertical.get_trailing_pixel_accum_ptr(1, -1, MultiAccumRowBuffer::PASS_TWO);
        evaluate_path( col, row, col+1, row-1,
                       prior_accum_ptr, full_prior_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );
                       
        AccumCostType* const prior_accum_ptr2 = buff_manager_vertical.get_trailing_pixel_accum_ptr(1, 1, MultiAccumRowBuffer::PASS_TWO);
        evaluate_path( col, row, col+1, row+1,
                       prior_accum_ptr2, full_prior_ptr, local_cost_ptr, temp_buffer.get(), 
                       pixel_diff, debug );
       for (int d=0; d<num_disp; ++d)
          output_accum_ptr[d] = (output_accum_ptr[d] + temp_buffer[d])/2;
      }
      else // Just init to the local cost
        for (int d=0; d<num_disp; ++d) output_accum_ptr[d] = local_cost_ptr[d];
      

      buff_manager_vertical.next_pixel();
    } // End col loop
    
    buff_manager_vertical.next_row(col==0);
  } // End row loop
  

  // Done with all trips!
} // End function smooth_path_accumulation







SemiGlobalMatcher::DisparityImage
SemiGlobalMatcher::semi_global_matching_func( ImageView<uint8> const& left_image,
                                              ImageView<uint8> const& right_image,
                                              ImageView<uint8> const* left_image_mask,
                                              ImageView<uint8> const* right_image_mask,
                                              DisparityImage const* prev_disparity,
                                              int search_buffer) {
                                              
  // Compute safe bounds to search through given the disparity range and kernel size.
  
  const int half_kernel_size = (m_kernel_size-1) / 2;

  // Using inclusive bounds here
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

  vw_out(DebugMessage, "stereo") << "Computed SGM cost bounding box: " << std::endl;
  vw_out(DebugMessage, "stereo") << "Left image size = ("<<left_image.cols()<<","<<left_image.rows()
                               <<"), right image size = ("<<right_image.cols()<<","<<right_image.rows()<<")\n";
  vw_out(DebugMessage, "stereo") << "min_row = "<< m_min_row <<", min_col = "<< m_min_col <<
                                  ", max_row = "<< m_max_row <<", max_col = "<< m_max_col <<
                                  ", output_height = "<< m_num_output_rows <<
                                  ", output_width = "<< m_num_output_cols <<"\n";

  populate_adjacent_disp_lookup_table();

  // By default the search bounds are the same for each pixel,
  //  but set them from the prior disparity image if the user passed it in.
  populate_constant_disp_bound_image();
  
  if (!populate_disp_bound_image(left_image_mask, right_image_mask, prev_disparity, search_buffer)) {
    vw_out(WarningMessage, "stereo") << "No valid pixels found in SGM input!.\n";
    // If the inputs are invalid, return a default disparity image.
    DisparityImage disparity( m_num_output_cols, m_num_output_rows );
    return invalidate_mask(disparity);
  }

  // All the hard work is done in the next few function calls!

  allocate_large_buffers();

  compute_disparity_costs(left_image, right_image);

  if (m_use_mgm)
    smooth_path_accumulation(left_image);
  else
    two_trip_path_accumulation(left_image);

  // Now that all the costs are calculated, fetch the best disparity for each pixel.
  //create_disparity_view_subpixel(); // DEBUG
  return create_disparity_view();
}



} // end namespace stereo
} // end namespace vw

