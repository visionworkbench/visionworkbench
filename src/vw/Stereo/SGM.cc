
#include <vw/Stereo/SGM.h>
#include <vw/Core/Debugging.h>
#include <vw/Cartography/GeoReferenceUtils.h>
/*
#if 1 //defined(VW_ENABLE_SSE) && (VW_ENABLE_SSE==1)
  #include <emmintrin.h>
  #include <smmintrin.h> // SSE4.1
#endif
*/
namespace vw {

namespace stereo {

// It is up to the user to perform bounds checking before using these functions!
// TODO: Move them somewhere, optimize the speed.
uint8 get_census_value_3x3(ImageView<uint8> const& image, int col, int row) {
  // This will be an 8 bit sequence.
  uint8 output = 0;
  uint8 center = image(col, row);
  if (image(col-1, row-1) > center) output += 128;
  if (image(col  , row-1) > center) output +=  64;
  if (image(col+1, row-1) > center) output +=  32;
  if (image(col-1, row  ) > center) output +=  16;
  if (image(col+1, row  ) > center) output +=   8;
  if (image(col-1, row+1) > center) output +=   4;
  if (image(col  , row+1) > center) output +=   2;
  if (image(col+1, row+1) > center) output +=   1;
  return output;
}
uint32 get_census_value_5x5(ImageView<uint8> const& image, int col, int row) {
  // This will be a 24 bit sequence.
  uint32 output = 0;
  uint32 addend = 1;
  uint32 center = image(col, row);
  for (int r=row+2; r>=row-2; --r) {
    for (int c=col+2; c>=col-2; --c) {
      if (r == c) continue;
      if (image(c,r) > center)
        output += addend;
      addend *=2;
    }
  }
  return output;
}

// TODO: Consolidate with code in Matcher.h!

/// Simple, unoptimized code for computing the hamming distance of two bytes.
size_t hamming_distance(unsigned char a, unsigned char b) {
    unsigned char dist = 0;
    unsigned char val = a ^ b; // XOR

    // Count the number of bits set
    while (val != 0) {
        // A bit is set, so increment the count and clear the bit
        ++dist;
        val &= val - 1;
    }
    return dist; // Return the number of differing bits
}

size_t hamming_distance(unsigned int a, unsigned int b) {
    unsigned int dist = 0;
    unsigned int val = a ^ b; // XOR

    // Count the number of bits set
    while (val != 0) {
        // A bit is set, so increment the count and clear the bit
        ++dist;
        val &= val - 1;
    }
    return dist; // Return the number of differing bits
}

//=========================================================================

void SemiGlobalMatcher::set_parameters(int min_disp_x, int min_disp_y,
                                      int max_disp_x, int max_disp_y,
                                      int kernel_size, uint16 p1, uint16 p2) {
  m_min_disp_x  = min_disp_x;
  m_min_disp_y  = min_disp_y;
  m_max_disp_x  = max_disp_x;
  m_max_disp_y  = max_disp_y;
  m_kernel_size = kernel_size;
  
  m_num_disp_x = m_max_disp_x - m_min_disp_x + 1;
  m_num_disp_y = m_max_disp_y - m_min_disp_y + 1;
  size_t size_check = m_num_disp_x * m_num_disp_y;;
  if (size_check > (size_t)std::numeric_limits<DisparityType>::max())
    vw_throw( NoImplErr() << "Number of disparities is too large for data type!\n" );
  m_num_disp   = m_num_disp_x * m_num_disp_y;   
  
  if (p1 > 0) // User provided
    m_p1 = p1; 
  else { // Choose based on the kernel size
    switch(kernel_size){
    case 3:  m_p1 = 2;  break; // Census transform 0-8
    case 5:  m_p1 = 7;  break; // Census transform 0-24
    default: m_p1 = 20; break; // 0-255 scale
    };
  }
  if (p2 > 0) // User provided
    m_p2 = p2;
  else {
    switch(kernel_size){
    case 3:  m_p2 = 50;  break; // Census transform 0-8
    case 5:  m_p2 = 100; break; // Census transform 0-24
    default: m_p2 = 250; break; // 0-255 scale
    };
  }
}


void SemiGlobalMatcher::populate_constant_disp_bound_image() {
  // Allocate the image
  m_disp_bound_image.set_size(m_num_output_cols, m_num_output_rows);
  // Fill it up with an identical vector
  Vector4i bounds_vector(m_min_disp_x, m_min_disp_y, m_max_disp_x, m_max_disp_y);
  size_t buffer_size = m_num_output_cols*m_num_output_rows;
  std::fill(m_disp_bound_image.data(), m_disp_bound_image.data()+buffer_size, bounds_vector); 
}

void SemiGlobalMatcher::populate_disp_bound_image(DisparityImage const* prev_disparity,
                                                  int search_buffer) {

  std::cout << "m_disp_bound_image" << bounding_box(m_disp_bound_image) << std::endl;

  const int SCALE_UP = 2;

  // There needs to be some "room" in the disparity search space for
  // us to discard prior results on the edge as not trustworthy predictors.
  const bool check_x_edge = ((m_max_disp_x - m_min_disp_x) > 1);
  const bool check_y_edge = ((m_max_disp_y - m_min_disp_y) > 1);

  double area = 0, percent_trusted = 0;

  int r_in, c_in;
  PixelMask<Vector2i> input_disp;
  Vector4i bounds;
  for (int r=0; r<m_disp_bound_image.rows(); ++r) {
    r_in = r / SCALE_UP;
    for (int c=0; c<m_disp_bound_image.cols(); ++c) {
      c_in = c / SCALE_UP;
      if ( (c_in >= prev_disparity->cols()) || 
           (r_in >= prev_disparity->rows())   ) {
        vw_throw( LogicErr() << "Size error!\n" );
      }
      input_disp = prev_disparity->operator()(c_in,r_in);
      
      // Disparity values on the edge of our 2D search range are not trustworthy!
      int dx_scaled = input_disp[0] * SCALE_UP; 
      int dy_scaled = input_disp[1] * SCALE_UP;
      
      bool on_edge = (  ( check_x_edge && ((dx_scaled <= m_min_disp_x) || (dx_scaled >= m_max_disp_x)) )
                     || ( check_y_edge && ((dy_scaled <= m_min_disp_y) || (dy_scaled >= m_max_disp_y)) ) );

      bool good_disparity = (is_valid(input_disp) && !on_edge);
      
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
    }
  }
  // Compute some statistics for help improving the speed
  double num_pixels = m_disp_bound_image.rows()*m_disp_bound_image.cols();
  area            /= num_pixels;
  percent_trusted /= num_pixels;
  double max_search_area = (m_max_disp_x-m_min_disp_x+1)*(m_max_disp_y-m_min_disp_y+1);
  std::cout << "Max pixel search area = " << max_search_area << std::endl;
  std::cout << "Mean pixel search area = " << area << std::endl;
  std::cout << "Percent trusted prior disparities = " << percent_trusted << std::endl;  
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
      m_adjacent_disp_lookup[table_pos+0] = y_less_o*m_num_disp_x + x_less_o;
      m_adjacent_disp_lookup[table_pos+1] = y_less_o*m_num_disp_x + x_o;
      m_adjacent_disp_lookup[table_pos+2] = y_less_o*m_num_disp_x + x_more_o;
      m_adjacent_disp_lookup[table_pos+3] = y_o     *m_num_disp_x + x_less_o;
      m_adjacent_disp_lookup[table_pos+4] = y_o     *m_num_disp_x + x_more_o;
      m_adjacent_disp_lookup[table_pos+5] = y_more_o*m_num_disp_x + x_less_o;
      m_adjacent_disp_lookup[table_pos+6] = y_more_o*m_num_disp_x + x_o;
      m_adjacent_disp_lookup[table_pos+7] = y_more_o*m_num_disp_x + x_more_o;
      
      ++d;
    }
  }
} // End function populate_adjacent_disp_lookup_table

void SemiGlobalMatcher::evaluate_path( int col, int row, int col_p, int row_p,
                       AccumCostType* const prior, // Accumulated costs leading up to this pixel
                       CostType     * const local, // The disparity costs of the current pixel
                       AccumCostType*       output,
                       int path_intensity_gradient, bool debug ) {

  // Decrease p2 (jump cost) with increasing disparity along the path
  AccumCostType p2_mod = m_p2;
  if (path_intensity_gradient > 0)
    p2_mod /= path_intensity_gradient;
  if (p2_mod < m_p1)
    p2_mod = m_p1;

  //// Output vector must be pre-initialized to default value!
  //AccumCostType min_score = output[0];

  Vector4i pixel_disp_bounds   = m_disp_bound_image(col, row);
  //Vector4i pixel_disp_bounds_p = m_disp_bound_image(col_p, row_p);
  
  // TODO: Get this from the previous iteration? SLOW
  //DisparityType min_prev_disparity_index = 0;
  AccumCostType min_prior = prior[0];
  for (DisparityType i=1; i<m_num_disp; ++i) {
    if (prior[i] < min_prior) {
      //min_prev_disparity_index = i;
      min_prior  = prior[i];
    }
  }

  AccumCostType min_prev_disparity_cost = min_prior + p2_mod;
  
//  if (debug) {
//    printf("P2_mod = %d, min_prev cost = %d, min_prev_index = %d\n", 
//           p2_mod, min_prev_disparity_cost, min_prev_disparity_index);
//  }

  // Loop through disparities for this pixel 
  for (int dy=pixel_disp_bounds[1]; dy<=pixel_disp_bounds[3]; ++dy) {
  
    // Get initial linear storage index for this dy row
    int d = (dy-m_min_disp_y)*m_num_disp_x +
            (pixel_disp_bounds[0] - m_min_disp_x);

    for (int dx=pixel_disp_bounds[0]; dx<=pixel_disp_bounds[2]; ++dx) {

      // Start with the cost for the same disparity
      AccumCostType lowest_combined_cost = prior[d];

      // Compare to the eight adjacent disparities using the lookup table
      const int LOOKUP_TABLE_WIDTH = 8;
      const int lookup_index = d*LOOKUP_TABLE_WIDTH;
/*      
//#if defined(VW_ENABLE_SSE) && (VW_ENABLE_SSE==1)
#if false

      // Load the prior values to compare into a buffer
      uint16 pre_sse_buff[8] __attribute__ ((aligned (16)));
      
      pre_sse_buff[0] = prior[m_adjacent_disp_lookup[lookup_index  ]];
      pre_sse_buff[1] = prior[m_adjacent_disp_lookup[lookup_index+1]];
      pre_sse_buff[2] = prior[m_adjacent_disp_lookup[lookup_index+2]];
      pre_sse_buff[3] = prior[m_adjacent_disp_lookup[lookup_index+3]];
      pre_sse_buff[4] = prior[m_adjacent_disp_lookup[lookup_index+4]];
      pre_sse_buff[5] = prior[m_adjacent_disp_lookup[lookup_index+5]];
      pre_sse_buff[6] = prior[m_adjacent_disp_lookup[lookup_index+6]];
      pre_sse_buff[7] = prior[m_adjacent_disp_lookup[lookup_index+7]];

      __m128i _sse2_reg = _mm_load_si128( (__m128i*) pre_sse_buff );

      _mm_minpos_epu16(_sse2_reg);
      
      _mm_store_si128( (__m128i*) pre_sse_buff, _sse2_reg );

      AccumCostType lowest_adjacent_cost = pre_sse_buff[0];

#else
*/
      AccumCostType lowest_adjacent_cost = prior[m_adjacent_disp_lookup[lookup_index]];
      lowest_adjacent_cost = std::min(lowest_adjacent_cost, prior[m_adjacent_disp_lookup[lookup_index+1]]);
      lowest_adjacent_cost = std::min(lowest_adjacent_cost, prior[m_adjacent_disp_lookup[lookup_index+2]]);
      lowest_adjacent_cost = std::min(lowest_adjacent_cost, prior[m_adjacent_disp_lookup[lookup_index+3]]);
      lowest_adjacent_cost = std::min(lowest_adjacent_cost, prior[m_adjacent_disp_lookup[lookup_index+4]]);
      lowest_adjacent_cost = std::min(lowest_adjacent_cost, prior[m_adjacent_disp_lookup[lookup_index+5]]);
      lowest_adjacent_cost = std::min(lowest_adjacent_cost, prior[m_adjacent_disp_lookup[lookup_index+6]]);
      lowest_adjacent_cost = std::min(lowest_adjacent_cost, prior[m_adjacent_disp_lookup[lookup_index+7]]);
      
//#endif
      
      // Now add the adjacent penalty cost and compare to the local cost
      lowest_adjacent_cost += m_p1;
      lowest_combined_cost = std::min(lowest_combined_cost, lowest_adjacent_cost);
      
      // Compare to the lowest prev disparity cost regardless of location
      lowest_combined_cost = std::min(lowest_combined_cost, min_prev_disparity_cost);
      /*
      if (debug && (d == 101)) {
        std::cout << "101\n";
        std::cout << "local = " <<  local[d] << std::endl;
        std::cout << "same = " <<  prior[d] << std::endl;
        for (int i=0; i<8; ++i)
          std::cout << i << " = " << m_adjacent_disp_lookup[lookup_index+i]<< " = "  
                    << prior[m_adjacent_disp_lookup[lookup_index+i]] +P_A << std::endl;
        std::cout << "min_prev_disparity_cost = " << min_prev_disparity_cost << std::endl;
      }
      */
      
      // The output cost = local cost + lowest combined cost - min_prior
      // - Subtracting out min_prior avoids overflow.
      output[d] = local[d] + lowest_combined_cost - min_prior;
      //if (output[d] < min_score)
      //  min_score = output[d];
    
      ++d;
    }
  } // End loop through this disparity
  
/*  
  if (debug) {
    printf("PATH costs (%d,%d): ", col, row);
    for (int i=0; i<m_num_disp; ++i)
      std::cout << output[i] << " ";
    std::cout << std::endl;
  }
*/  
  //return min_score;
}



SemiGlobalMatcher::DisparityImage
SemiGlobalMatcher::create_disparity_view( boost::shared_array<AccumCostType> const accumulated_costs ) {
  // Init output vector
  DisparityImage disparity( m_num_output_cols, m_num_output_rows );
  // For each element in the accumulated costs matrix, 
  //  select the disparity with the lowest accumulated cost.
  //Timer timer("\tCalculate Disparity Minimum");
  DisparityType dx, dy;
  for ( int j = 0; j < m_num_output_rows; j++ ) {
    for ( int i = 0; i < m_num_output_cols; i++ ) {
      
      const AccumCostType * vec_ptr = get_accum_vector(accumulated_costs, i,j);
      DisparityType disp = find_min_index(vec_ptr);
      
      disp_to_xy(disp, dx, dy); // TODO: Use a lookup table?
      disparity(i,j) = DisparityImage::pixel_type(dx, dy);      
      
      //if ((i >= 195) && (i <= 196) && (j==203)) {
      //  printf("ACC costs (%d,%d): %d, %d, %d\n", i, j, disp, dx, dy);
      //  print_disparity_vector(vec_ptr);
      //}
    }
  }
  printf("Done creating SGM result of size: %d, %d\n", disparity.cols(), disparity.rows());
  return disparity;
}

// TODO: Use our existing block matching code to do this!
// - Do this later, cost propagation is the real bottleneck.
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
  for ( int r = m_min_row; r < m_max_row; r++ ) { // For each row in left
    for ( int c = m_min_col; c < m_max_col; c++ ) { // For each column in left
        
      size_t d=0;
      for ( int dy = m_min_disp_y; dy <= m_max_disp_y; dy++ ) { // For each disparity
        for ( int dx = m_min_disp_x; dx <= m_max_disp_x; dx++ ) {
          bool debug = false;
          size_t cost_index = get_cost_index(c-m_min_col, r-m_min_row, d);
          CostType cost = get_cost_block(left_image, right_image, c, r, c+dx,r+dy, debug);
          m_cost_buffer[cost_index] = cost;
          ++d; // Disparity values are stored in a vector for convenience.
        }    
      } // End disparity loops   
    } // End x loop
  }// End y loop
}

// TODO: Speed up the census computations!

void SemiGlobalMatcher::fill_costs_census3x3   (ImageView<uint8> const& left_image,
                                                ImageView<uint8> const& right_image){
  // Compute the census value for each pixel.
  // - ROI handling could be fancier but this is simple and works.
  // - The 0,0 pixels in the left and right images are assumed to be aligned.
  ImageView<uint8> left_census (left_image.cols()-2,  left_image.rows()-2 ), 
                   right_census(right_image.cols()-2, right_image.rows()-2);
                   
  for ( int r = 0; r < left_census.rows(); r++ )
    for ( int c = 0; c < left_census.cols(); c++ )
      left_census(c,r) = get_census_value_3x3(left_image, c+1, r+1);
  for ( int r = 0; r < right_census.rows(); r++ )
    for ( int c = 0; c < right_census.cols(); c++ )
      right_census(c,r) = get_census_value_3x3(right_image, c+1, r+1);

  
  // Now compute the disparity costs for each pixel.
  // Make sure we don't go out of bounds here due to the disparity shift and kernel.
  for ( int r = m_min_row; r < m_max_row; r++ ) { // For each row in left
    for ( int c = m_min_col; c < m_max_col; c++ ) { // For each column in left
      size_t d=0;
      for ( int dy = m_min_disp_y; dy <= m_max_disp_y; dy++ ) { // For each disparity
        for ( int dx = m_min_disp_x; dx <= m_max_disp_x; dx++ ) {
          
          size_t   cost_index = get_cost_index(c-m_min_col, r-m_min_row, d);
          CostType cost = hamming_distance(left_census(c-1,r-1), right_census(c+dx-1, r+dy-1));
          m_cost_buffer[cost_index] = cost;
          ++d; // Disparity values are stored in a vector for convenience.
        }    
      } // End disparity loops   
    } // End x loop
  }// End y loop
}

void SemiGlobalMatcher::fill_costs_census5x5   (ImageView<uint8> const& left_image,
                                                ImageView<uint8> const& right_image){
  // Compute the census value for each pixel.
  // - ROI handling could be fancier but this is simple and works.
  // - The 0,0 pixels in the left and right images are assumed to be aligned.
  ImageView<uint8> left_census (left_image.cols()-4,  left_image.rows()-4 ), 
                   right_census(right_image.cols()-4, right_image.rows()-4);
                   
  for ( int r = 0; r < left_census.rows(); r++ )
    for ( int c = 0; c < left_census.cols(); c++ )
      left_census(c,r) = get_census_value_5x5(left_image, c+2, r+2);
  for ( int r = 0; r < right_census.rows(); r++ )
    for ( int c = 0; c < right_census.cols(); c++ )
      right_census(c,r) = get_census_value_5x5(right_image, c+2, r+2);

  
  // Now compute the disparity costs for each pixel.
  // Make sure we don't go out of bounds here due to the disparity shift and kernel.
  for ( int r = m_min_row; r < m_max_row; r++ ) { // For each row in left
    for ( int c = m_min_col; c < m_max_col; c++ ) { // For each column in left
      size_t d=0;
      for ( int dy = m_min_disp_y; dy <= m_max_disp_y; dy++ ) { // For each disparity
        for ( int dx = m_min_disp_x; dx <= m_max_disp_x; dx++ ) {
          
          size_t   cost_index = get_cost_index(c-m_min_col, r-m_min_row, d);
          CostType cost = hamming_distance(left_census(c-2,r-2), right_census(c+dx-2, r+dy-2));
          m_cost_buffer[cost_index] = cost;
          ++d; // Disparity values are stored in a vector for convenience.
        }    
      } // End disparity loops   
    } // End x loop
  }// End y loop
}



void SemiGlobalMatcher::compute_disparity_costs(ImageView<uint8> const& left_image,
                                                ImageView<uint8> const& right_image) {

  const CostType MAX_COST = std::numeric_limits<CostType>::max();

  // Init storage for processing all costs. W*H*D. D= NUM_DISPS
  // - Init all costs to the max possible value.
  size_t num_output_pixels = m_num_output_cols*m_num_output_rows;
  size_t num_cost_elements = num_output_pixels*m_num_disp;
  
  m_cost_buffer.reset(new CostType[num_cost_elements]);
  for (size_t i=0; i<num_cost_elements; ++i)
    m_cost_buffer[i] = MAX_COST;
  
  // Now that the buffers are initialized, choose how to fill them based on the kernel size.
  {
    Timer timer("\tCost Calculation");

    switch(m_kernel_size) {
    case 3:  fill_costs_census3x3(left_image, right_image); break;
    case 5:  fill_costs_census5x5(left_image, right_image); break;
    default: fill_costs_block    (left_image, right_image); break;
    };    
  } // End non-census case

  std::cout << "Done computing local costs.\n";

  // DEBUG!!!!

  //if ((m_num_output_cols > 196) && (m_num_output_rows > 203)) {
  //  std::cout << "Local cost vector for 196, 203:\n";
  //  print_disparity_vector(get_cost_vector(196,203));
  //}
/*
  int debug_row = 203;//_num_output_rows * 0.065;
  if (m_num_output_cols > debug_row) {
    std::cout << "debug row = " << debug_row << std::endl;
    ImageView<uint8> cost_image( m_num_output_cols, m_num_disp ); // TODO: Change type?
    for ( int i = 0; i < m_num_output_cols; i++ ) {
      for ( int d = 0; d < m_num_disp; d++ ) {
        size_t cost_index = get_cost_index(i, debug_row, d);
        cost_image(i,d) = m_cost_buffer[cost_index]; // TODO: scale!
      }
    }
    write_image("scanline_costs_block.tif",cost_image);
    std::cout << "Done writing line dump.\n";
  }
  */
  
} // end compute_disparity_costs() 



void SemiGlobalMatcher::two_pass_path_accumulation(ImageView<uint8> const& left_image) {

  // Instantiate two single-row buffers that will be used to temporarily store
  //  accumulated cost info until it is no longer needed.
  // - Within each buffer, data is indexed in order [col][pass][disparity]
  const size_t NUM_PATHS_IN_PASS = 4;
  const size_t buffer_pixel_size = NUM_PATHS_IN_PASS*m_num_disp;
  const size_t buffer_size       = m_num_output_cols*buffer_pixel_size;
  const size_t buffer_size_bytes = buffer_size*sizeof(AccumCostType);
  
  const int last_column = m_num_output_cols - 1;
  const int last_row    = m_num_output_rows - 1;

  // Allocate both buffers
  boost::shared_array<AccumCostType> bufferA, bufferB;
  bufferA.reset(new AccumCostType[buffer_size]);
  bufferB.reset(new AccumCostType[buffer_size]);
    
  AccumCostType* top_buffer = bufferA.get();
  AccumCostType* bot_buffer = bufferB.get();

  AccumCostType BAD_VAL = 2*std::numeric_limits<CostType>::max();

  // First pass, raster top left to bottom right.
  memset(top_buffer, 0, buffer_size_bytes);
  
  for (int row=0; row<m_num_output_rows; ++row) {
  
    // Init the bottom buffer to the default value
    std::fill(bot_buffer, bot_buffer+buffer_size, BAD_VAL);
  
    for (int col=0; col<m_num_output_cols; ++col) {
    
      // Set some pointers for this pixel
      CostType     * const local_cost_ptr   = get_cost_vector(col, row);
      AccumCostType*       output_accum_ptr = bot_buffer + col*buffer_pixel_size;
      bool debug = false;
      
      // Top left
      if ((row > 0) && (col > 0)) {
        // Fill in the accumulated value in the bottom buffer
        int pixel_diff = get_path_pixel_diff(left_image, col, row, 1, 1);
        AccumCostType* const prior_accum_ptr = top_buffer + (col-1)*buffer_pixel_size;
        evaluate_path( col, row, col-1, row-1,
                       prior_accum_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );
      }
      else // Just init to the local cost
        copy_cost_vector(local_cost_ptr, output_accum_ptr);

      output_accum_ptr += m_num_disp; // Move to the next path accumulation location
      
      // Top
      if (row > 0) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, 0, 1);
        AccumCostType* const prior_accum_ptr = top_buffer + col*buffer_pixel_size+1*m_num_disp;
        evaluate_path( col, row, col, row-1,
                       prior_accum_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );
      }
      else // Just init to the local cost
        copy_cost_vector(local_cost_ptr, output_accum_ptr);
      output_accum_ptr += m_num_disp;
      
      // Top right
      if ((row > 0) && (col < last_column)) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, -1, 1);
        AccumCostType* const prior_accum_ptr = top_buffer + (col+1)*buffer_pixel_size+2*m_num_disp;
        evaluate_path( col, row, col+1, row-1,
                       prior_accum_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );
      }
      else // Just init to the local cost
        copy_cost_vector(local_cost_ptr, output_accum_ptr);
      output_accum_ptr += m_num_disp;
      
      // Left
      if (col > 0) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, 1, 0);
        AccumCostType* const prior_accum_ptr = output_accum_ptr - buffer_pixel_size;
        evaluate_path( col, row, col-1, row,
                       prior_accum_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );
      }
      else // Just init to the local cost
        copy_cost_vector(local_cost_ptr, output_accum_ptr);
    
    } // End col loop
    
    // Sum up the contents of the bottom row of the buffer into m_accum_buffer
    update_accum_buffer_row(bot_buffer, row);
    
    std::swap(top_buffer, bot_buffer); // Swap the buffers
    
  } // End row loop

  // Second pass, raster bottom left to top right.
  // - Note that the roles of the top and bottom buffers are reversed here
  memset(bot_buffer, 0, buffer_size_bytes);

  for (int row = last_row; row >= 0; --row) {
  
    // Init the top buffer to the default value
    std::fill(top_buffer, top_buffer+buffer_size, BAD_VAL);
  
    for (int col = last_column; col >= 0; --col) {
    
      // Set some pointers for this pixel
      CostType     * const local_cost_ptr   = get_cost_vector(col, row);
      AccumCostType*       output_accum_ptr = top_buffer + col*buffer_pixel_size;
      bool debug = false;
              
      // Bottom right
      if ((row < last_row) && (col < last_column)) {
        // Fill in the accumulated value in the bottom buffer
        int pixel_diff = get_path_pixel_diff(left_image, col, row, -1, -1);
        AccumCostType* const prior_accum_ptr = bot_buffer + (col+1)*buffer_pixel_size;
        evaluate_path( col, row, col+1, row+1,
                       prior_accum_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );
      }
      else // Just init to the local cost
        copy_cost_vector(local_cost_ptr, output_accum_ptr);
      output_accum_ptr += m_num_disp; // Move to the next path accumulation location
      
      // Bottom
      if (row < last_row) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, 0, -1);
        AccumCostType* const prior_accum_ptr = bot_buffer + col*buffer_pixel_size+1*m_num_disp;
        evaluate_path( col, row, col, row+1,
                       prior_accum_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );
      }
      else // Just init to the local cost
        copy_cost_vector(local_cost_ptr, output_accum_ptr);
      output_accum_ptr += m_num_disp;
      
      // Bottom left
      if ((row < last_row) && (col > 0)) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, 1, -1);
        AccumCostType* const prior_accum_ptr = bot_buffer + (col-1)*buffer_pixel_size+2*m_num_disp;
        evaluate_path( col, row, col-1, row+1,
                       prior_accum_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );
      }
      else // Just init to the local cost
        copy_cost_vector(local_cost_ptr, output_accum_ptr);
      output_accum_ptr += m_num_disp;
      
      // Right
      if (col < last_column) {
        int pixel_diff = get_path_pixel_diff(left_image, col, row, -1, 0);
        AccumCostType* const prior_accum_ptr = output_accum_ptr + buffer_pixel_size;
        evaluate_path( col, row, col+1, row,
                       prior_accum_ptr, local_cost_ptr, output_accum_ptr, 
                       pixel_diff, debug );
      }
      else // Just init to the local cost
        copy_cost_vector(local_cost_ptr, output_accum_ptr);
    
    } // End col loop
    
    // Sum up the contents of the top row of the buffer into m_accum_buffer
    update_accum_buffer_row(top_buffer, row);
    
    std::swap(top_buffer, bot_buffer); // Swap the buffers
    
  } // End row loop

  // Done with both passes!
}





SemiGlobalMatcher::DisparityImage
SemiGlobalMatcher::semi_global_matching_func( ImageView<uint8> const& left_image,
                                              ImageView<uint8> const& right_image,
                                              DisparityImage const* prev_disparity,
                                              int search_buffer) {
                                              
  // Compute safe bounds to search through given the disparity range and kernel size.
  
  const int half_kernel_size = (m_kernel_size-1) / 2;

  m_min_row = half_kernel_size - m_min_disp_y; // Assumes the (0,0) pixels are aligned
  m_min_col = half_kernel_size - m_min_disp_x;
  m_max_row = std::min(left_image.rows()  -  half_kernel_size,
                       right_image.rows() - (half_kernel_size + m_max_disp_y) );
  m_max_col = std::min(left_image.cols()  -  half_kernel_size,
                       right_image.cols() - (half_kernel_size + m_max_disp_x) );
  if (m_min_row < 0) m_min_row = 0;
  if (m_min_col < 0) m_min_col = 0;
  if (m_max_row > left_image.rows()) m_max_row = left_image.rows();
  if (m_max_col > left_image.cols()) m_max_col = left_image.cols();

  m_num_output_cols  = m_max_col - m_min_col + 1;
  m_num_output_rows  = m_max_row - m_min_row + 1;

  std::cout << "Computed cost bounding box: " << std::endl;
  printf("Left image size = (%d,%d), right image size = (%d, %d)\n", 
         left_image.cols(), left_image.rows(), right_image.cols(), right_image.rows());
  printf("min_row = %d, min_col = %d, max_row = %d, max_col = %d, output_height = %d, output_width = %d\n", 
          m_min_row, m_min_col, m_max_row, m_max_col, m_num_output_rows, m_num_output_cols);
  
  m_buffer_step_size = m_num_output_cols * m_num_disp;

  populate_adjacent_disp_lookup_table();

  // By default the search bounds are the same for each image,
  //  but set them from the prior disparity image if the user passed it in.
  populate_constant_disp_bound_image();
  if (prev_disparity) {
    std::cout << "Updating bound image from previous disparity.\n";
    populate_disp_bound_image(prev_disparity, search_buffer);
  }

  // ==== Compute the cost values ====
  compute_disparity_costs(left_image, right_image);

  // ==== Accumulate the global costs across paths ====

  size_t num_output_pixels = m_num_output_cols*m_num_output_rows;
  size_t num_cost_elements = num_output_pixels*m_num_disp;  

  // TODO: Check available memory first!
  // Set up large buffer required by SGM algorithm
  std::cout << "Num pixels      = " << m_num_output_rows * m_num_output_cols << std::endl;
  std::cout << "Num disparities = " << m_num_disp << std::endl;
  std::cout << "Allocating buffer of size: " << (num_cost_elements*sizeof(AccumCostType))/(1024*1024) << " MB\n";
  //if (m_num_output_rows > 1000)
  //  vw_throw( NoImplErr() << "DEBUG!\n" );
  m_accum_buffer.reset(new AccumCostType[num_cost_elements]);
  for (size_t i=0; i<num_cost_elements; ++i)
    m_accum_buffer[i] = 0;

  std::cout << "Done allocating memory.\n";

  {
    Timer timer_total("\tCost Propagation");
    
    // This function will go through all eight accumulation paths in two passes.
    two_pass_path_accumulation(left_image);
    
    /*
    // Accumulate costs in each of the 8 cardinal directions.
    // - For each one: Fill the storage structure with zeros.
    //               : Accumulate costs in this direction.
    //               : Write debug image.
    //               : Accumulate the final direction costs into the total accumulated costs matrix.
    std::cout << "Processing lines: \n";
    {
      Timer timer("\tCost Propagation [1,0]");
      iterate_direction<1,0>( left_image, dir_accumulated_costs );
      //write_image("effect_1_0.tif", create_disparity_view( dir_accumulated_costs ) );
      inplace_sum_views( m_accum_buffer, dir_accumulated_costs );
    }
    std::cout << "1, ";
    
    iterate_direction<-1,0>( left_image, dir_accumulated_costs );
    //write_image("effect_n1_0.tif", create_disparity_view( dir_accumulated_costs ) );
    inplace_sum_views( m_accum_buffer, dir_accumulated_costs );
    std::cout << "2, ";
    
    iterate_direction<0,1>( left_image, dir_accumulated_costs );
    //write_image("effect_0_1.tif", create_disparity_view( dir_accumulated_costs ) );
    inplace_sum_views( m_accum_buffer, dir_accumulated_costs );
    std::cout << "3, ";
    
    iterate_direction<0,-1>( left_image, dir_accumulated_costs );
    //write_image("effect_0_n1.tif", create_disparity_view( dir_accumulated_costs ) );
    inplace_sum_views( m_accum_buffer, dir_accumulated_costs );
    std::cout << "4, ";
    
    iterate_direction<1,1>( left_image, dir_accumulated_costs );
    //write_image("effect_1_1.tif", create_disparity_view( dir_accumulated_costs ) );
    inplace_sum_views( m_accum_buffer, dir_accumulated_costs );
    std::cout << "5, ";
    
    iterate_direction<-1,-1>( left_image, dir_accumulated_costs );
    //write_image("effect_n1_n1.tif", create_disparity_view( dir_accumulated_costs ) );
    inplace_sum_views( m_accum_buffer, dir_accumulated_costs );
    std::cout << "6, ";
    
    iterate_direction<1,-1>( left_image, dir_accumulated_costs );
    //write_image("effect_1_n1.tif", create_disparity_view( dir_accumulated_costs ) );
    inplace_sum_views( m_accum_buffer, dir_accumulated_costs );
    std::cout << "7, ";
    
    iterate_direction<-1,1>( left_image, dir_accumulated_costs );
    //write_image("effect_n1_1.tif", create_disparity_view( dir_accumulated_costs ) );
    inplace_sum_views( m_accum_buffer, dir_accumulated_costs );
    std::cout << "8\n";
    */
  }

  // Now that all the costs are calculated, fetch the best disparity for each pixel.
  return create_disparity_view( m_accum_buffer );
}



} // end namespace stereo
} // end namespace vw

