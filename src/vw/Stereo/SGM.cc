
#include <vw/Stereo/SGM.h>
#include <vw/Core/Debugging.h>
#include <vw/Cartography/GeoReferenceUtils.h>

namespace vw {

namespace stereo {


void SemiGlobalMatcher::setParameters(int min_disp_x, int min_disp_y,
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

void SemiGlobalMatcher::evaluate_path( 
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
  
  for ( int32 d = 0; d < m_num_disp; d++ ) { // Loop through disparities in current pixel
    
    // For this disparity, find the prior disparity with the lowest cost increase.
    AccumCostType lowest_combined_cost = std::numeric_limits<int32>::max();
    
    for ( int32 d_p = 0; d_p < m_num_disp; d_p++ ) { // Loop through disparities in prior pixel
    
      AccumCostType prior_cost     = prior[d_p];                 // The cost to reach the disparity at the previous pixel
      AccumCostType disparity_diff = get_disparity_dist(d, d_p); // Physical disparity distance from that disparity to this one
      
      // This penalty scheme is too simplistic for large possible disparity jumps which we can get with 2D disparity!
      // Get something reasonable here, don't want to spend too much time tweaking up front.
    
      // TODO!
      /*
      AccumCostType penalty = 0;
      if ( disparity_diff == 0 ) {
        // No penalty, zero integer disparity
      } else if ( disparity_diff <= 2 ) {
        // Small penalty, eight-connected pixel disparity
        penalty = PENALTY1;
      } else {
         // Large penalty, >1 pixel disparity
         // - The larger the intensity difference between the previous pixel and this one,
         //   the smaller the assessed penalty is for a disparity jump.
         penalty = std::max(PENALTY1, path_intensity_gradient ? PENALTY2/path_intensity_gradient : PENALTY2);
      }
      */
     
     // TODO: Consider intensity diff?
     const AccumCostType P_A = 20; // TODO: This needs to be related to the local cost sizes!
     AccumCostType penalty  = disparity_diff * P_A;
      
     
     
      lowest_combined_cost = std::min(lowest_combined_cost, prior_cost+penalty);

      if ((debug)){
        printf("d: %d, d_p %d, curr_cost = %d, lcc = %d, prior_cost = %d, penalty = %d\n", 
               d, d_p, curr_cost[d], lowest_combined_cost, prior_cost, penalty);
      }
      

           
    } // End loop through prior disparity
    
    // The output cost = local cost + lowest combined cost
    curr_cost[d] += lowest_combined_cost;
    
    
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
  ImageView<Vector<DisparityType,2> > disparity( m_num_output_cols, m_num_output_rows );
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
      int SCALE = 1; // For debug display!
      disparity(i,j) = SCALE*Vector<DisparityType,2>(dx, dy);      
      
      //if ((i == 219) && (j==173))
      //printf("ACC costs (%d,%d): %d, %d, %d\n", i, j, disp, dx, dy);
    }
  }
  std::cout << "Done creating view\n";
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
                                              ImageView<uint8> const& right_image ) {

  std::cout << "Init!\n";
  
  const CostType MAX_COST       = std::numeric_limits<     CostType>::max();
  const CostType MAX_ACCUM_COST = std::numeric_limits<AccumCostType>::max();

  m_num_output_cols  = left_image.cols();
  m_num_output_rows  = left_image.rows();
  m_buffer_step_size = m_num_output_cols * m_num_disp;

  // Calc one result for each left image pixel.
  Vector2i size(m_num_output_cols, m_num_output_rows);

  // Init storage for processing all costs. W*H*D. D= NUM_DISPS
  // - Init all costs to the max possible value.
  size_t num_output_pixels = size.x() * size.y();
  size_t num_cost_elements = num_output_pixels*m_num_disp;
  
  m_cost_buffer.reset(new CostType[num_cost_elements]);
  for (size_t i=0; i<num_cost_elements; ++i)
    m_cost_buffer[i] = MAX_COST;

  
  std::cout << "Done filling cost matrices\n";
  
  // For each pixel, for each disparity, compute the cost of that choice.
  // - OPT: Could use different disparity ranges for each pixel to reduce the amount of work done.
  {
    Timer timer("\tCost Calculation");

    // Compute safe bounds to search through.
    // - Any location we don't compute a cost for defaults to the max cost,
    //   essentially ignoring that location.
    // - 
    
    const int half_kernel_size = (m_kernel_size-1) / 2;

    int min_row = half_kernel_size - m_min_disp_y; // Assumes the (0,0) pixels are aligned
    int min_col = half_kernel_size - m_min_disp_x;
    int max_row = std::min(left_image.rows()  -  half_kernel_size,
                           right_image.rows() - (half_kernel_size + m_max_disp_y) );
    int max_col = std::min(left_image.cols()  -  half_kernel_size,
                           right_image.cols() - (half_kernel_size + m_max_disp_x) );
    if (min_row < 0) min_row = 0;
    if (min_col < 0) min_col = 0;
    if (max_row > left_image.rows()) min_row = left_image.rows();
    if (max_col > left_image.cols()) min_col = left_image.cols();

    std::cout << "Computed cost bounding box: " << std::endl;
    printf("min_row = %d, min_col = %d, max_row = %d, max_col = %d\n", min_row, min_col, max_row, max_col);

    size_t d=0;    
    for ( int dy = m_min_disp_y; dy <= m_max_disp_y; dy++ ) { // For each disparity
      for ( int dx = m_min_disp_x; dx <= m_max_disp_x; dx++ ) {
  
        // Make sure we don't go out of bounds here due to the disparity shift and kernel.
        for ( int j = min_row; j < max_row; j++ ) { // For each row in left
          for ( int i = min_col; i < max_col; i++ ) { // For each column in left
          
            // Currently this computes only the SINGLE PIXEL difference.
            // Could improve results by converting to block matching around the center pixel.
          
            // TODO: Add or subtract the disparity?
            bool debug = ((i == 73) && (j == 159));
            size_t cost_index = get_cost_index(i, j, d);
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

  // Now that we have the local costs, accumulate the global costs.
  
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
      write_image("effect_1_0.tif", create_disparity_view( dir_accumulated_costs ) );
      inplace_sum_views( m_accum_buffer, dir_accumulated_costs );
    }
    
    std::cout << "1, ";
    {
      //Timer timer("\tCost Propagation [-1,0]");
      iterate_direction<-1,0>( left_image, dir_accumulated_costs );
      //write_image("effect_-1_0.tif", create_disparity_view( dir_accumulated_costs ) );
      inplace_sum_views( m_accum_buffer, dir_accumulated_costs );
    }
    std::cout << "2, ";
    {
      //Timer timer("\tCost Propagation [0,1]");
      iterate_direction<0,1>( left_image, dir_accumulated_costs );
      //write_image("effect_0_1.tif", create_disparity_view( dir_accumulated_costs ) );
      inplace_sum_views( m_accum_buffer, dir_accumulated_costs );
    }
    std::cout << "3, ";
    {
      //Timer timer("\tCost Propagation [0,-1]");
      iterate_direction<0,-1>( left_image, dir_accumulated_costs );
      //write_image("effect_0_-1.tif", create_disparity_view( dir_accumulated_costs ) );
      inplace_sum_views( m_accum_buffer, dir_accumulated_costs );
    }
    
    std::cout << "4, ";
    {
      //Timer timer("\tCost Propagation [1,1]");
      iterate_direction<1,1>( left_image, dir_accumulated_costs );
      //write_image("effect_1_1.tif", create_disparity_view( dir_accumulated_costs ) );
      inplace_sum_views( m_accum_buffer, dir_accumulated_costs );
    }
    std::cout << "5, ";
    {
      //Timer timer("\tCost Propagation [-1,-1]");
      iterate_direction<-1,-1>( left_image, dir_accumulated_costs );
      //write_image("effect_-1_-1.tif", create_disparity_view( dir_accumulated_costs ) );
      inplace_sum_views( m_accum_buffer, dir_accumulated_costs );
    }
    std::cout << "6, ";
    {
      //Timer timer("\tCost Propagation [1,-1]");
      iterate_direction<1,-1>( left_image, dir_accumulated_costs );
      //write_image("effect_1_-1.tif", create_disparity_view( dir_accumulated_costs ) );
      inplace_sum_views( m_accum_buffer, dir_accumulated_costs );
    }
    std::cout << "7, ";
    {
      //Timer timer("\tCost Propagation [-1,1]");
      iterate_direction<-1,1>( left_image, dir_accumulated_costs );
      //write_image("effect_-1_1.tif", create_disparity_view( dir_accumulated_costs ) );
      inplace_sum_views( m_accum_buffer, dir_accumulated_costs );
    }
    std::cout << "8\n";
    
  }

  // Now that all the costs are calculated, fetch the best disparity for each pixel.
  return create_disparity_view( m_accum_buffer );
}


} // end namespace stereo
} // end namespace vw

