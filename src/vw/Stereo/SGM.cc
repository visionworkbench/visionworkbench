
#include <vw/Stereo/SGM.h>
#include <vw/Core/Debugging.h>
#include <vw/Cartography/GeoReferenceUtils.h>

namespace vw {

namespace stereo {


int32 get_disparity_dist(int32 d1, int32 d2) {

  //// This is the traditional way for 1D disparities
  //return(abs(d_p - d));
  
  // This method is for 2D disparities
  // - Currently using truncated euclidean distance!
  int32 dx1, dy1, dx2, dy2;
  disp_to_xy(d1, dx1, dy1);
  disp_to_xy(d2, dx2, dy2);
  int32 delX = dx2 - dx1;
  int32 delY = dy2 - dy1;
  //return sqrt(delX*delX + delY*delY);
  //return delX*delX + delY*delY;
  
  return abs(delX) + abs(delY); // How does this method compare?

}

AVector evaluate_path( AVector const& prior,
                       CVector const& local,
                       int path_intensity_gradient, bool debug ) {
                       
  // Two parameters controlling disparity smoothness, see the paper.
  const int32 PENALTY1 = 15;
  const int32 PENALTY2 = 100;
                       
  // Init the output costs to the local costs
  AVector curr_cost = local;
  
  for ( int32 d = 0; d < NUM_DISPS; d++ ) { // Loop through disparities in current pixel
    
    // For this disparity, find the prior disparity with the lowest cost increase.
    int32 lowest_combined_cost = std::numeric_limits<int32>::max();
    
    for ( int32 d_p = 0; d_p < NUM_DISPS; d_p++ ) { // Loop through disparities in prior pixel
    
      int32 prior_cost     = prior[d_p];                 // The cost to reach the disparity at the previous pixel
      int32 disparity_diff = get_disparity_dist(d, d_p); // Physical disparity distance from that disparity to this one
      
      // This penalty scheme is too simplistic for large possible disparity jumps which we can get with 2D disparity!
      // Get something reasonable here, don't want to spend too much time tweaking up front.
    
      // TODO!
      
      int32 penalty = 0;
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
      
     /*
     // TODO: Consider intensity diff?
     const int32 P_A = 20; // TODO: This needs to be related to the local cost sizes!
     int32 penalty  = disparity_diff * P_A;
      */
     
     
      lowest_combined_cost = std::min(lowest_combined_cost, prior_cost+penalty);

      if ((debug)){
        printf("d: %d, d_p %d, curr_cost = %d, lcc = %d, prior_cost = %d, penalty = %d\n", d, d_p, curr_cost[d], lowest_combined_cost, prior_cost, penalty);
      }
      

      //if ((d == 18) && (d_p == 18) && (lowest_combined_cost != 0)) {
      //  printf("curr_cost = %d, lcc = %d, prior_cost = %d, penalty = %d\n", curr_cost[d], lowest_combined_cost, prior_cost, penalty);
      //}
           
    } // End loop through prior disparity
    
    // The output cost = local cost + lowest combined cost
    curr_cost[d] += lowest_combined_cost;
  } // End loop through this disparity

  // Normalize by subtracting min of prior cost, this avoids overflow.
  return elem_diff(curr_cost, min(prior));
}



ImageView<Vector<uint8,2> >
create_disparity_view( ImageView<AVector> const& accumulated_costs ) {
  // Init output vector
  ImageView<Vector<uint8,2> > disparity( accumulated_costs.cols(),
                                        accumulated_costs.rows() );
  // For each element in the accumulated costs matrix, 
  //  select the disparity with the lowest accumulated cost.
  //Timer timer("\tCalculate Disparity Minimum");
  int32 dx, dy;
  for ( size_t j = 0; j < disparity.rows(); j++ ) {
    for ( size_t i = 0; i < disparity.cols(); i++ ) {
      
      // TODO: Handle type conversions!
      int32 disp = find_min_index(accumulated_costs(i,j));
      disp_to_xy(disp, dx, dy);
      int SCALE = 1; // For debug display!
      disparity(i,j) = SCALE*Vector<uint8,2>(dx, dy);      
      
      if ((i == 219) && (j==173))
        std::cout << "ACC costs: " << accumulated_costs(i,j) <<", "<< disp <<", "<< dx << ", "<< dy << std::endl;
    }
  }
  std::cout << "Done creating view\n";
  return disparity;
}

// TODO: Use our existing block matching code to do this!
//  - Do this later, cost propagation is the real bottleneck.
int16 get_cost(ImageView<uint8> const& left_image,
               ImageView<uint8> const& right_image,
               int left_x, int left_y, int right_x, int right_y, bool debug) {
  // Single pixel diff
  int16 cost = abs( int(left_image(left_x,left_y)) - int(right_image(right_x,right_y)) );
  return cost;
  
  /*
  // Block mean of abs dists
  const int half_kernel_size = (SGM_KERNEL_SIZE-1) / 2;
  int sum = 0, diff=0;
  for (int j=-half_kernel_size; j<=half_kernel_size; ++j) {
    for (int i=-half_kernel_size; i<=half_kernel_size; ++i) {
      sum += abs(int(left_image(left_x+i, left_y+j)) - int(right_image(right_x+i, right_y+j)));
    }
  }*/
  //int16 result = int16(sum / (SGM_KERNEL_SIZE*SGM_KERNEL_SIZE));;
  //printf("sum = %d, result = %d\n", sum, result);
  //return int16(sum);
}

ImageView<Vector<uint8,2> >
semi_global_matching_func( ImageView<uint8> const& left_image,
                           ImageView<uint8> const& right_image ) {

  std::cout << "Init!\n";

  // Calc one result for each left image pixel.
  Vector2i size( left_image.cols(), left_image.rows() );

  // Init storage for processing all costs. W*H*D. D= NUM_DISPS
  size_t num_output_pixels = size.x() * size.y();
  ImageView<CVector > costs( size.x(), size.y() );

  // Fill all values in "costs" with 255.
  CVector temporary;
  std::fill(&temporary[0], &temporary[0]+NUM_DISPS, 255u);
  std::fill(costs.data(),costs.data()+num_output_pixels, temporary);
  
  std::cout << "Done filling cost matrices\n";
  
  // For each pixel, for each disparity, compute the cost of that choice.
  // - OPT: Could use different disparity ranges for each pixel to reduce the amount of work done.
  {
    Timer timer("\tCost Calculation");

    const int half_kernel_size = (SGM_KERNEL_SIZE-1) / 2;

    size_t d=0;    
    for ( int dy = 0; dy < DISP_RANGE_Y; dy++ ) { // For each disparity
      for ( int dx = 0; dx < DISP_RANGE_X; dx++ ) {
  
        // Make sure we don't go out of bounds here due to the disparity shift.
        for ( int j = half_kernel_size; j < (size.y()-half_kernel_size)-dy; j++ ) { // For each row in left
          for ( int i = half_kernel_size; i < (size.x()-half_kernel_size)-dx; i++ ) { // For each column in left
          
            // Currently this computes only the SINGLE PIXEL difference.
            // Could improve results by converting to block matching around the center pixel.
          
            // TODO: Add or subtract the disparity?
            bool debug = ((i == 219) && (j == 173));
            costs(i,j)[d] = get_cost(left_image, right_image, i, j, i+dx,j+dy, debug);
            
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
  ImageView<uint8> cost_image( left_image.cols(), NUM_DISPS );
  for ( size_t i = 0; i < left_image.cols(); i++ ) {
    for ( size_t j = 0; j < NUM_DISPS; j++ ) {
      cost_image(i,j) = costs(i,185)[j];
    }
  }
  write_image("scanline_costs_block.tif",cost_image);
  std::cout << "Done writing line dump.\n";

  // Now that we have the local costs, accumulate the global costs.
  ImageView<AVector> accumulated_costs    ( left_image.cols(), left_image.rows() );
  ImageView<AVector> dir_accumulated_costs( left_image.cols(), left_image.rows() );
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
      iterate_direction<1,0>( left_image, costs, dir_accumulated_costs );
      write_image("effect_1_0.tif", create_disparity_view( dir_accumulated_costs ) );
      inplace_sum_views( accumulated_costs, dir_accumulated_costs );
    }
    
    std::cout << "1, ";
    {
      //Timer timer("\tCost Propagation [-1,0]");
      iterate_direction<-1,0>( left_image, costs, dir_accumulated_costs );
      //write_image("effect_-1_0.tif", create_disparity_view( dir_accumulated_costs ) );
      inplace_sum_views( accumulated_costs, dir_accumulated_costs );
    }
    std::cout << "2, ";
    {
      //Timer timer("\tCost Propagation [0,1]");
      iterate_direction<0,1>( left_image, costs, dir_accumulated_costs );
      //write_image("effect_0_1.tif", create_disparity_view( dir_accumulated_costs ) );
      inplace_sum_views( accumulated_costs, dir_accumulated_costs );
    }
    std::cout << "3, ";
    {
      //Timer timer("\tCost Propagation [0,-1]");
      iterate_direction<0,-1>( left_image, costs, dir_accumulated_costs );
      //write_image("effect_0_-1.tif", create_disparity_view( dir_accumulated_costs ) );
      inplace_sum_views( accumulated_costs, dir_accumulated_costs );
    }
    
    std::cout << "4, ";
    {
      //Timer timer("\tCost Propagation [1,1]");
      iterate_direction<1,1>( left_image, costs, dir_accumulated_costs );
      //write_image("effect_1_1.tif", create_disparity_view( dir_accumulated_costs ) );
      inplace_sum_views( accumulated_costs, dir_accumulated_costs );
    }
    std::cout << "5, ";
    {
      //Timer timer("\tCost Propagation [-1,-1]");
      iterate_direction<-1,-1>( left_image, costs, dir_accumulated_costs );
      //write_image("effect_-1_-1.tif", create_disparity_view( dir_accumulated_costs ) );
      inplace_sum_views( accumulated_costs, dir_accumulated_costs );
    }
    std::cout << "6, ";
    {
      //Timer timer("\tCost Propagation [1,-1]");
      iterate_direction<1,-1>( left_image, costs, dir_accumulated_costs );
      //write_image("effect_1_-1.tif", create_disparity_view( dir_accumulated_costs ) );
      inplace_sum_views( accumulated_costs, dir_accumulated_costs );
    }
    std::cout << "7, ";
    {
      //Timer timer("\tCost Propagation [-1,1]");
      iterate_direction<-1,1>( left_image, costs, dir_accumulated_costs );
      //write_image("effect_-1_1.tif", create_disparity_view( dir_accumulated_costs ) );
      inplace_sum_views( accumulated_costs, dir_accumulated_costs );
    }
    std::cout << "8\n";
    
  }

  // Now that all the costs are calculated, fetch the best disparity for each pixel.
  return create_disparity_view( accumulated_costs );
}


} // end namespace stereo
} // end namespace vw

