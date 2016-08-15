
#include <vw/Stereo/SGM.h>
#include <vw/Core/Debugging.h>
#include <vw/Cartography/GeoReferenceUtils.h>

namespace vw {

namespace stereo {


void SemiGlobalMatcher::set_parameters(int min_disp_x, int min_disp_y,
                                      int max_disp_x, int max_disp_y,
                                      int kernel_size) {
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
    
    
    
}


void SemiGlobalMatcher::populate_constant_disp_bound_image() {
  // Allocate the image
  m_disp_bound_image.set_size(m_num_output_cols, m_num_output_rows);
  // Fill it up with an identical vector
  Vector4i bounds_vector(m_min_disp_x, m_min_disp_y, m_max_disp_x, m_max_disp_y);
  size_t buffer_size = m_num_output_cols*m_num_output_rows;
  std::fill(m_disp_bound_image.data(), m_disp_bound_image.data()+buffer_size, bounds_vector); 
}

void SemiGlobalMatcher::populate_disp_bound_image(std::vector<stereo::SearchParam> const *zones) {
  // populate_disp_bound_image should already be allocated with constant values.
  
  // Loop through all vzones
  size_t num_zones = zones->size();
  for (size_t i=0; i<num_zones; ++i) {
  
    // Get disparity vector for this zone
    stereo::SearchParam zone = zones->operator[](i);
    Vector4i zone_disp(zone.disparity_range().min().x(),
                       zone.disparity_range().min().y(),
                       zone.disparity_range().max().x(),
                       zone.disparity_range().max().y() );
                       
    // Copy this vector to the region specified by the zone.
    // - Note that the zone coordinates are in original left image coordinates, and
    //   we need to fill in the (smaller) output image region.
    int min_col = zone.image_region().min().x() - m_min_col; 
    int min_row = zone.image_region().min().y() - m_min_row;
    int max_col = zone.image_region().max().x() - m_min_col;
    int max_row = zone.image_region().max().y() - m_min_row;
    if (min_col < 0) min_col = 0;
    if (min_row < 0) min_row = 0;
    if (max_col >= m_num_output_cols) max_col = m_num_output_cols-1;
    if (max_row >= m_num_output_rows) max_row = m_num_output_rows-1;
    for (int row=min_row; row<max_row; ++row) {
      for (int col=min_col; col<max_col; ++col) {
        m_disp_bound_image(col, row) = zone_disp;
      }
    } // End loop through zone region
  } // End loop through zones
}

void SemiGlobalMatcher::populate_disp_bound_image(DisparityImage const* prev_disparity,
                                                  int search_buffer) {

  std::cout << "m_disp_bound_image" << bounding_box(m_disp_bound_image) << std::endl;

  const int SCALE_UP = 2;

  double area = 0, percent_trusted = 0;

  // TODO: Revisit bounds with kernel sizes!
  // TODO: Expand the bounds when the pixel is invalid!
  // TOOD: Also handle case where we are at max search range.
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
      // TODO: May need to relax this for the Y direction!
      bool on_edge = ( (input_disp[0] == m_min_disp_x) ||
                       (input_disp[1] == m_min_disp_y) ||
                       (input_disp[0] == m_max_disp_x) ||
                       (input_disp[1] == m_max_disp_y)   );
      bool good_disparity = (is_valid(input_disp) && !on_edge);
      
      if (good_disparity) {
        // We are more confident in the prior disparity, search nearby.
        bounds[0]  = input_disp[0]*SCALE_UP - search_buffer; // Min x
        bounds[2]  = input_disp[0]*SCALE_UP + search_buffer; // Max X
        bounds[1]  = input_disp[1]*SCALE_UP - search_buffer; // Min y
        bounds[3]  = input_disp[1]*SCALE_UP + search_buffer; // Max y
        
        // Constrain to global limits
        if (bounds[0] < m_min_disp_x) bounds[0] = m_min_disp_x;
        if (bounds[1] < m_min_disp_y) bounds[1] = m_min_disp_y;
        if (bounds[2] > m_max_disp_x) bounds[2] = m_max_disp_x;
        if (bounds[3] > m_max_disp_y) bounds[3] = m_max_disp_y;
      
        percent_trusted += 1.0;
        //std::cout << input_disp << " --> " << bounds << std::endl;
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
/*
SemiGlobalMatcher::AccumCostType SemiGlobalMatcher::get_disparity_dist(DisparityType d1, DisparityType d2) {

  //// This is the traditional way for 1D disparities
  //return(abs(d_p - d));
  
  // This method is for 2D disparities
  // - Currently using truncated euclidean distance!
  DisparityType dx1, dy1, dx2, dy2;
  disp_to_xy(d1, dx1, dy1);
  disp_to_xy(d2, dx2, dy2);
  AccumCostType delX = dx2 - dx1;
  AccumCostType delY = dy2 - dy1;
  //return sqrt(delX*delX + delY*delY);
  //return delX*delX + delY*delY;
  
  return abs(delX) + abs(delY); // How does this method compare?
}*/


void SemiGlobalMatcher::populate_adjacent_disp_lookup_table() {

  std::cout << "Generating lookup table...\n";
  
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
  std::cout << "Finished generating lookup table.\n";
}

void SemiGlobalMatcher::evaluate_path( int col, int row, int col_p, int row_p,
                       AccumCostType* const prior, // Accumulated costs leading up to this pixel
                       CostType     * const local, // The disparity costs of the current pixel
                       AccumCostType*       output,
                       int path_intensity_gradient, bool debug ) {

  // TODO: Consider intensity diff
  AccumCostType P_A = 20; // Dist = 1
  AccumCostType P_B = 20; // Dist = 2
  AccumCostType P_C = 200; // Dist > 2
  
  if (path_intensity_gradient > 0)
    P_C /= path_intensity_gradient;
  if (P_C < P_B)
    P_C = P_B;

  // Init the output costs to a large value so that we don't
  //  use disparity combinations that we don't compute.
  std::vector<AccumCostType> curr_cost(m_num_disp);
  for (int i=0; i<m_num_disp; ++i)
    curr_cost[i] = 2*std::numeric_limits<CostType>::max();
  
  Vector4i pixel_disp_bounds   = m_disp_bound_image(col,   row);
  Vector4i pixel_disp_bounds_p = m_disp_bound_image(col_p, row_p);
  

  // TODO: Get this from the previous iteration!
  DisparityType min_prev_disparity_index = 0;
  CostType      min_prev_disparity_cost  = prior[0];
  for (DisparityType i=1; i<m_num_disp; ++i) {
    if (prior[i] < min_prev_disparity_cost) {
      min_prev_disparity_index = i;
      min_prev_disparity_cost  = prior[i];
    }
  }
  
  
  AccumCostType penalty;

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
      int lookup_index = d*LOOKUP_TABLE_WIDTH;
      AccumCostType lowest_adjacent_cost = prior[m_adjacent_disp_lookup[lookup_index+0]];
      lowest_adjacent_cost = std::min(lowest_adjacent_cost, prior[m_adjacent_disp_lookup[lookup_index+1]]);
      lowest_adjacent_cost = std::min(lowest_adjacent_cost, prior[m_adjacent_disp_lookup[lookup_index+2]]);
      lowest_adjacent_cost = std::min(lowest_adjacent_cost, prior[m_adjacent_disp_lookup[lookup_index+3]]);
      lowest_adjacent_cost = std::min(lowest_adjacent_cost, prior[m_adjacent_disp_lookup[lookup_index+4]]);
      lowest_adjacent_cost = std::min(lowest_adjacent_cost, prior[m_adjacent_disp_lookup[lookup_index+5]]);
      lowest_adjacent_cost = std::min(lowest_adjacent_cost, prior[m_adjacent_disp_lookup[lookup_index+6]]);
      lowest_adjacent_cost = std::min(lowest_adjacent_cost, prior[m_adjacent_disp_lookup[lookup_index+7]]);
      // Now add the adjacent penalty cost
      lowest_combined_cost = std::min(lowest_combined_cost, lowest_adjacent_cost+P_A);
      
      // Compare to the lowest prev disparity cost regardless of location
      lowest_combined_cost = std::min(lowest_combined_cost, min_prev_disparity_cost+P_C);
      
      
      // The output cost = local cost + lowest combined cost
      curr_cost[d] = local[d] + lowest_combined_cost;
    
      ++d;
    }
  } // End loop through this disparity
  

  // Normalize by subtracting min of prior cost, this avoids overflow.
  
  AccumCostType min_prior = get_accum_vector_min(prior);
  //print_disparity_vector(prior);
  //std::cout << "min_prior = " << min_prior << " output = " ;
  for (int i=0; i<m_num_disp; ++i) {
    output[i] = (curr_cost[i] - min_prior);
    //std::cout << output[i] << " ";
  }
  //std::cout << "\n" ;
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
      //print_disparity_vector(vec_ptr);
      DisparityType disp = find_min_index(vec_ptr);
      
      disp_to_xy(disp, dx, dy);
      disparity(i,j) = DisparityImage::pixel_type(dx, dy);      
      
      //if ((i == 219) && (j==173))
      //printf("ACC costs (%d,%d): %d, %d, %d\n", i, j, disp, dx, dy);
    }
  }
  printf("Done creating SGM result of size: %d, %d\n", disparity.cols(), disparity.rows());
  return disparity;
}

// TODO: Use our existing block matching code to do this!
//  - Do this later, cost propagation is the real bottleneck.
SemiGlobalMatcher::CostType SemiGlobalMatcher::get_cost(ImageView<uint8> const& left_image,
               ImageView<uint8> const& right_image,
               int left_x, int left_y, int right_x, int right_y, bool debug) {
  //// Single pixel diff
  //CostType cost = abs( int(left_image(left_x,left_y)) - int(right_image(right_x,right_y)) );
  //return cost;
  
  
  // Block mean of abs dists
  const int half_kernel_size = (m_kernel_size-1) / 2;
  int sum = 0, diff=0;
  for (int j=-half_kernel_size; j<=half_kernel_size; ++j) {
    for (int i=-half_kernel_size; i<=half_kernel_size; ++i) {
      sum += abs(int(left_image(left_x+i, left_y+j)) - int(right_image(right_x+i, right_y+j)));
    }
  }
  CostType result = CostType(sum / (m_kernel_size*m_kernel_size));;
  //printf("sum = %d, result = %d\n", sum, result);
  return CostType(sum);
}

SemiGlobalMatcher::DisparityImage
SemiGlobalMatcher::semi_global_matching_func( ImageView<uint8> const& left_image,
                                              ImageView<uint8> const& right_image,
                                              DisparityImage const* prev_disparity,
                                              int search_buffer) {

  const CostType MAX_COST       = std::numeric_limits<     CostType>::max();
  const CostType MAX_ACCUM_COST = std::numeric_limits<AccumCostType>::max();



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
  //  but set them from the zones if the user passed them in.
  populate_constant_disp_bound_image();
  //if (zones)
  //  populate_disp_bound_image(zones);
  if (prev_disparity) {
    std::cout << "Updating bound image from previous disparity.\n";
    populate_disp_bound_image(prev_disparity, search_buffer);
  }


  // ==== Compute the cost values ====


  // Init storage for processing all costs. W*H*D. D= NUM_DISPS
  // - Init all costs to the max possible value.
  size_t num_output_pixels = m_num_output_cols*m_num_output_rows;
  size_t num_cost_elements = num_output_pixels*m_num_disp;
  
  m_cost_buffer.reset(new CostType[num_cost_elements]);
  for (size_t i=0; i<num_cost_elements; ++i)
    m_cost_buffer[i] = MAX_COST;
  
  std::cout << "Done filling cost matrices\n";
  
  // For each pixel, for each disparity, compute the cost of that choice.
  // - TODO: Ony compute the specified disparity costs for each pixel
  // - TODO: Optimize this if it becames a speed concern.
  {
    Timer timer("\tCost Calculation");

    // Make sure we don't go out of bounds here due to the disparity shift and kernel.
    for ( int j = m_min_row; j < m_max_row; j++ ) { // For each row in left
      for ( int i = m_min_col; i < m_max_col; i++ ) { // For each column in left
          
        size_t d=0;
        //CostType min_cost = std::numeric_limits<CostType>::max();
        for ( int dy = m_min_disp_y; dy <= m_max_disp_y; dy++ ) { // For each disparity
          for ( int dx = m_min_disp_x; dx <= m_max_disp_x; dx++ ) {
            
            bool debug = ((i == 73) && (j == 159));
            size_t cost_index = get_cost_index(i-m_min_col, j-m_min_col, d);
            CostType cost = get_cost(left_image, right_image, i, j, i+dx,j+dy, debug);
            m_cost_buffer[cost_index] = cost;
            
            ++d; // Disparity values are stored in a vector for convenience.
          }    
        } // End disparity loops   

      } // End x loop
    }// End y loop
    
    //std::cout << "Loc 219, 173 costs --> " << costs(219,173) << std::endl;
  }

  std::cout << "Done computing local costs.\n";

  // Note: For test images, correct disp: dx=2, dy=1 --> d_index=18
/*
  // DEBUG!!!!
  int debug_row = m_num_output_rows / 2;
  ImageView<uint8> cost_image( m_num_output_cols, m_num_disp ); // TODO: Change type?
  for ( int i = 0; i < m_num_output_cols; i++ ) {
    for ( int d = 0; d < m_num_disp; d++ ) {
      size_t cost_index = get_cost_index(i, debug_row, d);
      cost_image(i,d) = m_cost_buffer[cost_index]; // TODO: scale!
    }
  }
  write_image("scanline_costs_block.tif",cost_image);
  std::cout << "Done writing line dump.\n";
*/

  // ==== Accumulate the global costs across paths ====
  
  // TODO: Possible to remove the temporary buffer?
  boost::shared_array<AccumCostType> dir_accumulated_costs; 
  dir_accumulated_costs.reset(new AccumCostType[num_cost_elements]);
  m_accum_buffer.reset(       new AccumCostType[num_cost_elements]);
  for (size_t i=0; i<num_cost_elements; ++i) {
    m_accum_buffer       [i] = 0;
    dir_accumulated_costs[i] = MAX_ACCUM_COST;
  }

  {
    Timer timer_total("\tCost Propagation");
    
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
    inplace_sum_views( m_accum_buffer, dir_accumulated_costs );
    std::cout << "2, ";
    
    iterate_direction<0,1>( left_image, dir_accumulated_costs );
    inplace_sum_views( m_accum_buffer, dir_accumulated_costs );
    std::cout << "3, ";
    
    iterate_direction<0,-1>( left_image, dir_accumulated_costs );
    inplace_sum_views( m_accum_buffer, dir_accumulated_costs );
    std::cout << "4, ";
    
    iterate_direction<1,1>( left_image, dir_accumulated_costs );
    inplace_sum_views( m_accum_buffer, dir_accumulated_costs );
    std::cout << "5, ";
    
    iterate_direction<-1,-1>( left_image, dir_accumulated_costs );
    inplace_sum_views( m_accum_buffer, dir_accumulated_costs );
    std::cout << "6, ";
    
    iterate_direction<1,-1>( left_image, dir_accumulated_costs );
    inplace_sum_views( m_accum_buffer, dir_accumulated_costs );
    std::cout << "7, ";
    
    iterate_direction<-1,1>( left_image, dir_accumulated_costs );
    inplace_sum_views( m_accum_buffer, dir_accumulated_costs );
    std::cout << "8\n";
    
  }

  // Now that all the costs are calculated, fetch the best disparity for each pixel.
  return create_disparity_view( m_accum_buffer );
}


} // end namespace stereo
} // end namespace vw

