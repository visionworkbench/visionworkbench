
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
      bounds[0]  = input_disp[0]*SCALE_UP - search_buffer; // Min x
      bounds[2]  = input_disp[0]*SCALE_UP + search_buffer; // Max X
      bounds[1]  = input_disp[1]*SCALE_UP - search_buffer; // Min y
      bounds[3]  = input_disp[1]*SCALE_UP + search_buffer; // Max y
      
      // Constrain to global limits
      if (bounds[0] < m_min_disp_x) bounds[0] = m_min_disp_x;
      if (bounds[1] < m_min_disp_y) bounds[1] = m_min_disp_y;
      if (bounds[2] > m_max_disp_x) bounds[2] = m_max_disp_x;
      if (bounds[3] > m_max_disp_y) bounds[3] = m_max_disp_y;
      
      std::cout << input_disp << " --> " << bounds << std::endl;
      
      //bounds = Vector4i(m_min_disp_x, m_min_disp_y, m_max_disp_x, m_max_disp_y); // DEBUG
      
      m_disp_bound_image(c,r) = bounds;
    }
  }
  
}

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

}

void SemiGlobalMatcher::evaluate_path( int col, int row, int col_p, int row_p,
                       AccumCostType* const prior, // Accumulated costs leading up to this pixel
                       CostType     * const local, // The disparity costs of the current pixel
                       AccumCostType*       output,
                       int path_intensity_gradient, bool debug ) {
                       
  // Two parameters controlling disparity smoothness, see the paper.
  // The penalty type is the same as the cost type
  const AccumCostType PENALTY1 = 15;
  const AccumCostType PENALTY2 = 100;
                       
  // Init the output costs to the local costs
  std::vector<AccumCostType> curr_cost(m_num_disp);
  for (int i=0; i<m_num_disp; ++i)
    curr_cost[i] = local[i];
  
  Vector4i pixel_disp_bounds   = m_disp_bound_image(col, row);
  Vector4i pixel_disp_bounds_p = m_disp_bound_image(col_p, row_p);

  // Loop through disparities for this pixel 
  for (int dy=pixel_disp_bounds[1]; dy<=pixel_disp_bounds[3]; ++dy) {
  
    // Get initial linear storage index for this dy row
    int d = (dy-m_min_disp_y)*m_num_disp_x +
            (pixel_disp_bounds[0] - m_min_disp_x);
  
    for (int dx=pixel_disp_bounds[0]; dx<=pixel_disp_bounds[2]; ++dx) {
      
      // For this disparity, find the prior disparity with the lowest cost increase.
      AccumCostType lowest_combined_cost = std::numeric_limits<int32>::max();

      // Loop through disparities in prior pixel
      for (int dy_p=pixel_disp_bounds_p[1]; dy_p<=pixel_disp_bounds_p[3]; ++dy_p) {
      
        AccumCostType y_disp_diff = abs(dy_p - dy); // Precompute part of distance
        
        // Get initial linear storage index for this dy row
        int d_p = (dy_p-m_min_disp_y)*m_num_disp_x +
                  (pixel_disp_bounds_p[0] - m_min_disp_x);
        
        for (int dx_p=pixel_disp_bounds_p[0]; dx_p<=pixel_disp_bounds_p[2]; ++dx_p) {
        
          // This is the costliest portion of the entire algorithm since it is executed so many times!
        
          // The cost to reach the disparity at the previous pixel       
          AccumCostType prior_cost     = prior[d_p];
          // Physical disparity distance from that disparity to this one   
          AccumCostType disparity_diff = abs(dx_p - dx) + y_disp_diff; // TODO: Experiment with this!
        
          // TODO: Consider intensity diff
          const AccumCostType P_A = 20; // TODO: This needs to be related to the local cost sizes!
          AccumCostType penalty  = disparity_diff * P_A;
         
          lowest_combined_cost = std::min(lowest_combined_cost, prior_cost+penalty);

          //if ((debug)){
          //  printf("d: %d, d_p %d, curr_cost = %d, lcc = %d, prior_cost = %d, penalty = %d\n", 
          //         d, d_p, curr_cost[d], lowest_combined_cost, prior_cost, penalty);
          //}
        
          ++d_p;
        }
      } // End loop through other pixel disparities
      
      // The output cost = local cost + lowest combined cost
      curr_cost[d] += lowest_combined_cost;
    
      ++d;
    }
  } // End loop through this disparity
  

  // Normalize by subtracting min of prior cost, this avoids overflow.
  
  AccumCostType min_prior = get_accum_vector_min(prior);
  //print_disparity_vector(prior);
  //std::cout << "min_prior = " << min_prior << " output = " ;
  for (int i=0; i<m_num_disp; ++i) {
    output[i] = curr_cost[i] - min_prior;
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

  // By default the search bounds are the same for each image,
  //  but set them from the zones if the user passed them in.
  populate_constant_disp_bound_image();
  //if (zones)
  //  populate_disp_bound_image(zones);
  //if (prev_disparity)
  //  populate_disp_bound_image(prev_disparity, search_buffer);


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
  // - OPT: Could use different disparity ranges for each pixel to reduce the amount of work done.
  {
    Timer timer("\tCost Calculation");


    size_t d=0;    
    for ( int dy = m_min_disp_y; dy <= m_max_disp_y; dy++ ) { // For each disparity
      for ( int dx = m_min_disp_x; dx <= m_max_disp_x; dx++ ) {
  
        // Make sure we don't go out of bounds here due to the disparity shift and kernel.
        for ( int j = m_min_row; j < m_max_row; j++ ) { // For each row in left
          for ( int i = m_min_col; i < m_max_col; i++ ) { // For each column in left
          
            // Currently this computes only the SINGLE PIXEL difference.
            // Could improve results by converting to block matching around the center pixel.
          
            // TODO: Add or subtract the disparity?
            bool debug = ((i == 73) && (j == 159));
            size_t cost_index = get_cost_index(i-m_min_col, j-m_min_col, d);
            m_cost_buffer[cost_index] = get_cost(left_image, right_image, i, j, i+dx,j+dy, debug);
            
            //if (debug) {
            //  printf("Done with loc %d, %d -->%d\n", i, j, costs(i,j)[d]);
            //}
            
          } // End x loop
        }// End y loop
        //printf("Done with disparity %d, %d\n", dx, dy);
        ++d; // Disparity values are stored in a vector for convenience.
      }    
    } // End disparity loops
    
    //std::cout << "Loc 219, 173 costs --> " << costs(219,173) << std::endl;
  }

  std::cout << "Done computing local costs.\n";

  // Note: For test images, correct disp: dx=2, dy=1 --> d_index=18
/*
  // DEBUG!!!!
  ImageView<uint8> cost_image( m_num_output_cols, m_num_disp );
  for ( int i = 0; i < m_num_output_cols; i++ ) {
    for ( int d = 0; d < m_num_disp; d++ ) {
      size_t cost_index = get_cost_index(i, 185, d);
      cost_image(i,d) = m_cost_buffer[cost_index]; // TODO: scale!
    }
  }
  write_image("scanline_costs_block.tif",cost_image);
  std::cout << "Done writing line dump.\n";
*/

  // ==== Accumulate the global costs across paths ====
  
  boost::shared_array<AccumCostType> dir_accumulated_costs; // TODO: Eliminate this!
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

