#ifndef __SEMI_GLOBAL_MATCHING_H__
#define __SEMI_GLBOAL_MATCHING_H__

#include <vw/Image/ImageView.h>
#include <vw/Image/PixelMask.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Stereo/Correlation.h>

#include <vw/InterestPoint/Detector.h> // TODO: REMOVE THIS!

namespace vw {


  // Registering the Pixel Disparity type for FileIO
  template<> struct PixelFormatID<Vector<uint8,2>     > { static const PixelFormatEnum value = VW_PIXEL_GENERIC_2_CHANNEL; };
  template<> struct PixelFormatID<Vector<int16,2>     > { static const PixelFormatEnum value = VW_PIXEL_GENERIC_2_CHANNEL; };

/* TODO LIST
  - Increase speed, it is way too slow!
  - Clean up interfaces/ROI handling.
  - Handle non-uint8 input images.
  - Integrate with disparity code as an option for the low-
    resolution disparity search.
    
NOTES:

 - Using a block size in the initial cost calculation does not seem to improve the results.
   Instead it gives splotchy output.  No real effect on speed.
 - All the required time is in the path iteration functions!
    
*/


namespace stereo {

class SemiGlobalMatcher {

public: // Definitions

  // The types are chosen to minimize storage costs
  typedef int16 DisparityType; ///< Contains the allowable dx, dy range.
  typedef uint8 CostType;      ///< Used to describe a single disparity cost.
  
  typedef uint16 AccumCostType; ///< Used to accumulate CostType values.

  typedef ImageView<PixelMask<Vector2i> > DisparityImage; // The usual VW disparity type

public: // Functions

  /// Set the parameters to be used for future SGM calls
  /// - Parameters that are not provided will be set to the best known default.
  /// - If kernel_size is 3 or 5, a census transform will be used (reccomended).
  ///   Otherwise a simple averaging over a block method will be used.
  void set_parameters(int min_disp_x, int min_disp_y,
                      int max_disp_x, int max_disp_y,
                      int kernel_size=5,
                      uint16 p1=0, uint16 p2=0);

  /// Compute SGM stereo on the images.
  /// - TODO: Make search_buffer a parameter
  DisparityImage
  semi_global_matching_func( ImageView<uint8> const& left_image,
                             ImageView<uint8> const& right_image,
                             DisparityImage const* prev_disparity=0,
                             int search_buffer = 2);

private: // Variables

    // The core parameters
    int m_min_disp_x, m_min_disp_y;
    int m_max_disp_x, m_max_disp_y;
    int m_kernel_size; ///< Must be odd. Use "1" for single pixel.

    int m_min_row, m_max_row;
    int m_min_col, m_max_col;
    int m_num_output_cols, m_num_output_rows;
    
    // Algorithm parameters
    AccumCostType m_p1;
    AccumCostType m_p2;
    
    // Derived parameters for convenience
    int m_num_disp_x, m_num_disp_y, m_num_disp;
    
    size_t m_buffer_step_size;
    boost::shared_array<CostType     > m_cost_buffer;
    boost::shared_array<AccumCostType> m_accum_buffer;
    
    /// Image containing the disparity bounds for each pixel.
    /// - Stored as min_col, min_row, max_col, max_row.
    ImageView<Vector4i> m_disp_bound_image;
    
    /// Lookup table of the adjacent disparities for each disparity
    /// - Handles outer boundaries by repitition.
    std::vector<DisparityType> m_adjacent_disp_lookup;

private: // Functions

  /// Populate the lookup table m_adjacent_disp_lookup
  void populate_adjacent_disp_lookup_table();

  /// Fill in m_disp_bound_image using image-wide contstants
  void populate_constant_disp_bound_image();

  /// Fill in m_disp_bound_image using a half resolution disparity image.
  void populate_disp_bound_image(DisparityImage const* prev_disparity,
                                 int search_buffer);

  /// Populates m_cost_buffer with all the disparity costs 
  void compute_disparity_costs(ImageView<uint8> const& left_image,
                               ImageView<uint8> const& right_image);
                               
  // The following functions are called from inside compute_disparity_costs()
  
  /// Compute mean of differences within a block of pixels.
  void fill_costs_block    (ImageView<uint8> const& left_image,
                            ImageView<uint8> const& right_image);
  // Compute a 8 bit hamming distance from a 5x5 census transform.
  void fill_costs_census3x3(ImageView<uint8> const& left_image,
                            ImageView<uint8> const& right_image);
  // Compute a 24 bit hamming distance from a 5x5 census transform.
  void fill_costs_census5x5(ImageView<uint8> const& left_image,
                            ImageView<uint8> const& right_image);

  /// Returns a block cost score at a given location
  CostType get_cost_block(ImageView<uint8> const& left_image,
                    ImageView<uint8> const& right_image,
                    int left_x, int left_y, int right_x, int right_y, bool debug);

  /// Return the index into one of the buffers for a given location
  /// - The data is stored row major interleaved format.
  size_t get_cost_index(int col, int row, DisparityType disp=0) const {
      return row*m_buffer_step_size + col*m_num_disp + disp;
  }

  /// Set provided buffer to the cost value at the selected pixel.
  void set_cost_vector(int col, int row,
                       boost::shared_array<AccumCostType> accum_vec) {
    size_t start_index = get_cost_index(col, row, 0);
    for (int d=0; d<m_num_disp; ++d)
      accum_vec[start_index+d] = m_cost_buffer[start_index+d];
  }
  
  CostType* get_cost_vector(int col, int row) {
    size_t start_index = get_cost_index(col, row);
    return &(m_cost_buffer[start_index]);
  };
  
  AccumCostType* get_accum_vector(boost::shared_array<AccumCostType> accum_vec,
                                  int col, int row) {
    size_t start_index = get_cost_index(col, row);
    return &(accum_vec[start_index]);
  };


  /// Goes across all the viterbi diagrams and extracts out the minimum vector.
  DisparityImage
  create_disparity_view( boost::shared_array<AccumCostType> const accumulated_costs );

  // TODO: Use a lookup table?
  /// Converts from a linear disparity index to the dx, dy values it represents.
  /// - This function is too slow to use inside the inner loop!
  void disp_to_xy(DisparityType disp, DisparityType &dx, DisparityType &dy) {
    dy = (disp / m_num_disp_x) + m_min_disp_y; // 2D implementation
    dx = (disp % m_num_disp_x) + m_min_disp_x;
  } 
  
  /// v1 += v2.
  void inplace_sum_views( boost::shared_array<AccumCostType>       v1,
                          boost::shared_array<AccumCostType> const v2) {
    size_t count = m_num_output_cols * m_num_output_rows * m_num_disp;
    for (size_t i=0; i<count; ++i)
      v1[i] = v1[i] + v2[i];
  }


  /// Get the index if the smallest element in a vector
  DisparityType find_min_index( const AccumCostType* const vec ) {
    AccumCostType value = std::numeric_limits<AccumCostType>::max();
    int min_index = 0;
    for (int i=0; i<m_num_disp; ++i) {
      if (vec[i] < value) {
        value = vec[i];
        min_index = i;
      }
    }
    return min_index;
  }
  
  /// Get the min disparity of an AccumCost vector
  AccumCostType get_accum_vector_min(const AccumCostType* const vec){
    AccumCostType value = std::numeric_limits<AccumCostType>::max();
    for (int i=0; i<m_num_disp; ++i) {
      if (vec[i] < value)
        value = vec[i];
    }
    return value;
  }

  // Print out a disparity vector
  template <typename T>
  void print_disparity_vector(T* const vec){
    std::cout << "V: ";
    for (int i=0; i<m_num_disp; ++i)
      std::cout << vec[i] << " ";
    std::cout << std::endl;
  }
  
  /// Get the pixel diff along a line at a specified output location.
  int get_path_pixel_diff(ImageView<uint8> const& left_image,
                          int col, int row, int dir_x, int dir_y) const {
    // Take the offset between the output location and the input pixel coordinates.
    return abs(left_image(col-m_min_col,         row-m_min_row) - 
               left_image((col-dir_x)-m_min_col, (row-dir_y)-m_min_row));
  }
  

  /// Create an updated cost accumulation vector for the next pixel along an SGM evaluation path.
  /// - For each disparity in the current pixel, add that disparity's cost with the "cheapest"
  ///   prior pixel disparity.
  void evaluate_path( int col, int row, int col_p, int row_p,
                      AccumCostType* const prior, // Accumulated costs leading up to this pixel
                      CostType     * const local, // The disparity costs of the current pixel
                      AccumCostType*       output,
                      int path_intensity_gradient, bool debug=false ); // This variable is the magnitude of intensity change to this pixel
/*
  // TODO: Move to cc file!
  /// Perform all eight path accumulations in two passes through the image
  void two_pass_path_accumulation(ImageView<uint8> const& left_image) {

    // Instantiate two single-row buffers that will be used to temporarily store
    //  accumulated cost info until it is no longer needed.
    // - Within each buffer, data is indexed in order [col][pass][disparity]
    const size_t NUM_PATHS_IN_PASS = 4;
    const size_t buffer_pixel_size = m_num_disp*NUM_PATHS_IN_PASS;
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

    // First pass, raster top left to bottom right.
    memset(top_buffer, 0, buffer_size_bytes);
    
    for (int row=0; row<m_num_output_rows; ++row) {
    
      // Init the bottom buffer to zero
      memset(bot_buffer, 0, buffer_size_bytes);
    
      for (int col=0; col<m_num_output_cols; ++col) {
      
        // Set some pointers for this pixel
        CostType     * const local_cost_ptr   = get_cost_vector(col, row);
        AccumCostType*       output_accum_ptr = bot_buffer + col*buffer_pixel_size;
        bool debug = false;
                
        // Top left
        if ((row > 0) && (col > 0) {
          // Fill in the accumulated value in the bottom buffer
          int pixel_diff = get_path_pixel_diff(left_image, col, row, 1, 1);
          AccumCostType* const prior_accum_ptr = top_buffer + (col-1)*buffer_pixel_size;
          evaluate_path( col, row, col-1, row-1,
                         prior_accum_ptr, local_cost_ptr, output_accum_ptr, 
                         pixel_diff, debug );
        }
        output_accum_ptr += m_num_disp; // Move to the next path accumulation location
        
        // Top
        if (row > 0) {
          int pixel_diff = get_path_pixel_diff(left_image, col, row, 0, 1);
          AccumCostType* const prior_accum_ptr = top_buffer + col*buffer_pixel_size;
          evaluate_path( col, row, col, row-1,
                         prior_accum_ptr, local_cost_ptr, output_accum_ptr, 
                         pixel_diff, debug );
        }       
        output_accum_ptr += m_num_disp;
        
        // Top right
        if ((row > 0) && (col < last_column) {
          int pixel_diff = get_path_pixel_diff(left_image, col, row, -1, 1);
          AccumCostType* const prior_accum_ptr = top_buffer + (col+1)*buffer_pixel_size;
          evaluate_path( col, row, col+1, row-1,
                         prior_accum_ptr, local_cost_ptr, output_accum_ptr, 
                         pixel_diff, debug );
        }
        output_accum_ptr += m_num_disp;
        
        // Left
        if (col > 0) {
          int pixel_diff = get_path_pixel_diff(left_image, col, row, 1, 0);
          AccumCostType* const prior_accum_ptr = output_accum_ptr - buffer_pixel_size;
          evaluate_path( col, row, col-1, row,
                         prior_accum_ptr, local_cost_ptr, output_accum_ptr, 
                         pixel_diff, debug );
        }
      
      } // End col loop
      
      // Sum up the contents of the top row of the buffer into m_accum_buffer
      size_t buffer_index = 0;
      for (int col=0; col<m_num_output_cols; ++col) {
        for (int pass=0; pass<NUM_PATHS_IN_PASS; ++pass) {
          size_t out_index = get_cost_index(col, row);
          for (int d=0; d<m_num_disp; ++d) {
            m_accum_buffer[out_index] += top_buffer[buffer_index];
            ++out_index;
            ++buffer_index;
          }
        }
      } // Done adding the temp buffer to m_accum_buffer
      
      std::swap(top_buffer, bot_buffer); // Swap the buffers
      
    } // End row loop
    
    std::cout << "DEBUG - first accum pass done!\n";
    return;
    
    // Second pass, raster bottom left to top right.
    // - Note that the roles of the top and bottom buffers are reversed here
    memset(bot_buffer, 0, buffer_size_bytes);

    for (int row = last_row; row >= 0; --row) {
    
      // Init the top buffer to zero
      memset(top_buffer, 0, buffer_size_bytes);
    
      for (int col = last_column; col >= 0; --col) {
      
        // Set some pointers for this pixel
        CostType     * const local_cost_ptr   = get_cost_vector(col, row);
        AccumCostType*       output_accum_ptr = top_buffer + col*buffer_pixel_size;
        bool debug = false;
                
        // Bottom right
        if ((row < last_row) && (col < last_column) {
          // Fill in the accumulated value in the bottom buffer
          int pixel_diff = get_path_pixel_diff(left_image, col, row, -1, -1);
          AccumCostType* const prior_accum_ptr = bot_buffer + (col+1)*buffer_pixel_size;
          evaluate_path( col, row, col+1, row+1,
                         prior_accum_ptr, local_cost_ptr, output_accum_ptr, 
                         pixel_diff, debug );
        }
        output_accum_ptr -= m_num_disp; // Move to the next path accumulation location
        
        // Bottom
        if (row < last_row) {
          int pixel_diff = get_path_pixel_diff(left_image, col, row, 0, -1);
          AccumCostType* const prior_accum_ptr = bot_buffer + col*buffer_pixel_size;
          evaluate_path( col, row, col, row+1,
                         prior_accum_ptr, local_cost_ptr, output_accum_ptr, 
                         pixel_diff, debug );
        }       
        output_accum_ptr -= m_num_disp;
        
        // Bottom left
        if ((row < last_row) && (col > 0) {
          int pixel_diff = get_path_pixel_diff(left_image, col, row, 1, -1);
          AccumCostType* const prior_accum_ptr = bot_buffer + (col-1)*buffer_pixel_size;
          evaluate_path( col, row, col-1, row+1,
                         prior_accum_ptr, local_cost_ptr, output_accum_ptr, 
                         pixel_diff, debug );
        }
        output_accum_ptr -= m_num_disp;
        
        // Right
        if (col < last_column) {
          int pixel_diff = get_path_pixel_diff(left_image, col, row, -1, 0);
          AccumCostType* const prior_accum_ptr = output_accum_ptr + buffer_pixel_size;
          evaluate_path( col, row, col+1, row,
                         prior_accum_ptr, local_cost_ptr, output_accum_ptr, 
                         pixel_diff, debug );
        }
      
      } // End col loop
      
      // Sum up the contents of the bottom row of the buffer into m_accum_buffer
      size_t buffer_index = 0;
      for (int col=0; col<m_num_output_cols; ++col) {
        for (int pass=0; pass<NUM_PATHS_IN_PASS; ++pass) {
          size_t out_index = get_cost_index(col, row);
          for (int d=0; d<m_num_disp; ++d) {
            m_accum_buffer[out_index] += bot_buffer[buffer_index];
            ++out_index;
            ++buffer_index;
          }
        }
      } // Done adding the temp buffer to m_accum_buffer
      
      std::swap(top_buffer, bot_buffer); // Swap the buffers
      
    } // End row loop

    std::cout << "DEBUG - second accum pass done!\n";

    // Done with both passes!
  }
*/

  /// Compute the accumulated costs in a pixel direction from the local costs at each pixel.
  /// - TODO: This implementation seems inefficient!
  template <int DIRX, int DIRY>
  void iterate_direction( ImageView<uint8   > const& left_image,
                          boost::shared_array<AccumCostType>      & accumulated_costs ) {


    // Clear the accumulation values before we write to them.
    size_t num_cost_elements = m_num_output_cols*m_num_output_rows*m_num_disp;
    for (size_t i=0; i<num_cost_elements; ++i)
      accumulated_costs[i] = 0;


    // Walk along the edges in a clockwise fashion
    if ( DIRX > 0 ) {
      // LEFT MOST EDGE
      // Init the edge pixels with just the cost (no accumulation yet)
      for ( int32 j = 0; j < m_num_output_rows; j++ ) {
        set_cost_vector(0, j, accumulated_costs);
      }
      // Loop across to the opposite edge
      for ( int32 i = 1; i < m_num_output_cols; i++ ) {
        // Loop through the pixels in this column, limiting the range according
        //  to the iteration direction progress.
        int32 jstart = std::max( 0,      0      + DIRY * i );
        int32 jstop  = std::min( m_num_output_rows, m_num_output_rows + DIRY * i );
        for ( int32 j = jstart; j < jstop; j++ ) {
        
          bool debug = false;//((i >= 195) && (i <= 196) && (j==203));
        
          int pixel_diff = get_path_pixel_diff(left_image, i, j, DIRX, DIRY);
          evaluate_path( i, j, i-DIRX,j-DIRY,
                         get_accum_vector(accumulated_costs, i-DIRX,j-DIRY), // Previous pixel
                         get_cost_vector(i,j),
                         get_accum_vector(accumulated_costs, i,j), 
                         pixel_diff, debug );         // Current pixel
        }
      }
    } 
    if ( DIRY > 0 ) {
      // TOP MOST EDGE
      // Process every pixel along this edge only if DIRX == 0. Otherwise skip the top left most pixel
      for ( int32 i = (DIRX <= 0 ? 0 : 1 ); i < m_num_output_cols; i++ ) {
        set_cost_vector(i, 0, accumulated_costs);
      }
      for ( int32 j = 1; j < m_num_output_rows; j++ ) {
        int32 istart = std::max( (DIRX <= 0 ? 0 : 1),
                                 (DIRX <= 0 ? 0 : 1) + DIRX * j );
        int32 istop  = std::min( m_num_output_cols, m_num_output_cols + DIRX * j );
        for ( int32 i = istart; i < istop; i++ ) {
          int pixel_diff = get_path_pixel_diff(left_image, i, j, DIRX, DIRY);
          evaluate_path( i, j, i-DIRX,j-DIRY,
                         get_accum_vector(accumulated_costs, i-DIRX,j-DIRY), // Previous pixel
                         get_cost_vector(i,j), 
                         get_accum_vector(accumulated_costs, i,j),
                         pixel_diff );         // Current pixel
        }
      }
    } 
    if ( DIRX < 0 ) {
      // RIGHT MOST EDGE
      // Process every pixel along this edge only if DIRY == 0. Otherwise skip the top right most pixel
      for ( int32 j = (DIRY <= 0 ? 0 : 1); j < m_num_output_rows; j++ ) {
        set_cost_vector(m_num_output_cols-1, j, accumulated_costs);
      }
      for ( int32 i = m_num_output_cols-2; i >= 0; i-- ) {
        int32 jstart = std::max( (DIRY <= 0 ? 0 : 1),
                                 (DIRY <= 0 ? 0 : 1) - DIRY * (i - m_num_output_cols + 1) );
        int32 jstop  = std::min( m_num_output_rows, m_num_output_rows - DIRY * (i - m_num_output_cols + 1) );
        for ( int32 j = jstart; j < jstop; j++ ) {
          int pixel_diff = get_path_pixel_diff(left_image, i, j, DIRX, DIRY);
          evaluate_path( i, j, i-DIRX,j-DIRY,
                         get_accum_vector(accumulated_costs, i-DIRX,j-DIRY), // Previous pixel
                         get_cost_vector(i,j), 
                         get_accum_vector(accumulated_costs, i,j),
                         pixel_diff );         // Current pixel
        }
      }
    } 
    if ( DIRY < 0 ) {
      // BOTTOM MOST EDGE
      // Process every pixel along this edge only if DIRX == 0. Otherwise skip the bottom left and bottom right pixel
      for ( int32 i = (DIRX <= 0 ? 0 : 1);
            i < (DIRX >= 0 ? m_num_output_cols : m_num_output_cols-1); i++ ) {
        set_cost_vector(i,m_num_output_rows-1, accumulated_costs);
      }
      for ( int32 j = m_num_output_rows-2; j >= 0; j-- ) {
        int32 istart = std::max( (DIRX <= 0 ? 0 : 1),
                                 (DIRX <= 0 ? 0 : 1) - DIRX * (j - m_num_output_rows + 1) );
        int32 istop  = std::min( (DIRX >= 0 ? m_num_output_cols : m_num_output_cols - 1),
                                 (DIRX >= 0 ? m_num_output_cols : m_num_output_cols - 1) - DIRX * (j - m_num_output_rows + 1) );
        for ( int32 i = istart; i < istop; i++ ) {
          int pixel_diff = get_path_pixel_diff(left_image, i, j, DIRX, DIRY);
          evaluate_path( i, j, i-DIRX,j-DIRY,
                         get_accum_vector(accumulated_costs, i-DIRX,j-DIRY), // Previous pixel
                         get_cost_vector(i,j), 
                         get_accum_vector(accumulated_costs, i,j),
                         pixel_diff );         // Current pixel
        }
      }
    }
  } // End function iterate_direction

}; // end class SemiGlobalMatcher


/// Wrapper function for SGM that handles ROIs.
/// - Merge with the function in Correlation.h!
/// - This function only searches positive disparities. The input images need to be
///   already cropped so that this makes sense.
template <class ImageT1, class ImageT2>
  ImageView<PixelMask<Vector2i> >
  calc_disparity_sgm(//CostFunctionType cost_type,
                 ImageViewBase<ImageT1> const& left_in,
                 ImageViewBase<ImageT2> const& right_in,
                 BBox2i                 const& left_region,   // Valid region in the left image
                 Vector2i               const& search_volume, // Max disparity to search in right image
                 Vector2i               const& kernel_size,  // Only really takes an N by N kernel!
                 SemiGlobalMatcher::DisparityImage  const* prev_disparity=0){ 

    
    // Sanity check the input:
    VW_DEBUG_ASSERT( kernel_size[0] % 2 == 1 && kernel_size[1] % 2 == 1,
                     ArgumentErr() << "best_of_search_convolution: Kernel input not sized with odd values." );
    VW_DEBUG_ASSERT( kernel_size[0] <= left_region.width() &&
                     kernel_size[1] <= left_region.height(),
                     ArgumentErr() << "best_of_search_convolution: Kernel size too large of active region." );
    VW_DEBUG_ASSERT( search_volume[0] > 0 && search_volume[1] > 0,
                     ArgumentErr() << "best_of_search_convolution: Search volume must be greater than 0." );
    VW_DEBUG_ASSERT( left_region.min().x() >= 0 &&  left_region.min().y() >= 0 &&
                     left_region.max().x() <= left_in.impl().cols() &&
                     left_region.max().y() <= left_in.impl().rows(),
                     ArgumentErr() << "best_of_search_convolution: Region not inside left image." );

    // Rasterize input so that we can do a lot of processing on it.
    BBox2i right_region = left_region;
    right_region.max() += search_volume - Vector2i(1,1); // Why the -1?
    
    std::cout << "calc_disparity_sgm: left  region  = " << left_region   << std::endl;
    std::cout << "calc_disparity_sgm: right region  = " << right_region  << std::endl;
    std::cout << "calc_disparity_sgm: search_volume = " << search_volume << std::endl;
    
    ImageView<typename ImageT1::pixel_type> left_crop ( crop(left_in.impl(),  left_region) );
    ImageView<typename ImageT2::pixel_type> right_crop( crop(right_in.impl(), right_region) );

    // TODO: Make scaling optional
    // Convert the input image to uint8 with 2%-98% intensity scaling.
    ImageView<PixelGray<vw::uint8> > left, right;
    ip::percentile_scale_convert(left_crop,  left,  0.02, 0.98);
    ip::percentile_scale_convert(right_crop, right, 0.02, 0.98);
    
    //write_image("final_sgm_left.tif", left);
    //write_image("final_sgm_right.tif", right);
    
    
    // TODO: Support different cost types?
    SemiGlobalMatcher matcher;
    matcher.set_parameters(0, 0, search_volume[0], search_volume[1], kernel_size[0]);
    return matcher.semi_global_matching_func(left, right, prev_disparity);
    
  } // End function calc_disparity


} // end namespace stereo
} // end namespace vw

#endif //__SEMI_GLOBAL_MATCHING__
