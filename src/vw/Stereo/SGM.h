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
  typedef int   DisparityType; ///< Contains the allowable dx, dy range.
  typedef uint8 CostType;      ///< Used to describe a single disparity cost.
  
  typedef uint16 AccumCostType; ///< Used to accumulate CostType values.

  typedef ImageView<PixelMask<Vector2i> > DisparityImage; // The usual VW disparity type

public: // Functions

  SemiGlobalMatcher() {}

  /// Set set_parameters for details
  SemiGlobalMatcher(int min_disp_x, int min_disp_y,
                    int max_disp_x, int max_disp_y,
                    int kernel_size=5,
                    uint16 p1=0, uint16 p2=0) {
    set_parameters(min_disp_x, min_disp_y, max_disp_x, max_disp_y, kernel_size, p1, p2);
  }

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
                             ImageView<uint8> const* left_image_mask,
                             ImageView<uint8> const* right_image_mask,
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
    
    //size_t m_buffer_step_size;
    boost::shared_array<CostType     > m_cost_buffer;
    boost::shared_array<AccumCostType> m_accum_buffer;
    
    /// Image containing the disparity bounds for each pixel.
    /// - Stored as min_col, min_row, max_col, max_row.
    ImageView<Vector4i> m_disp_bound_image;
    
    /// Lookup table of the adjacent disparities for each disparity
    /// - For each disparity index, store the disparity indices of the 
    ///   eight adjacent disparities.
    /// - Handles outer boundaries by repitition.
    /// - This vector stores a table of size m_num_disp*8.
    std::vector<DisparityType> m_adjacent_disp_lookup;

/* 
Proposed changes for sparse accum buffer:
- Each pixel only has space in the buffer to contain all the values 
  specified in m_disp_bound_image.  Will be either m_num_disp or something
  potentially much smaller.
- Entries for each pixel are contiguous.
- A new data structure is needed to store the starting index in m_accum_buffer
  for each pixel.
- Issue: The prior disparity vector may be small, and won't contain all the indices
         specified in m_adjacent_disp_lookup.  
       : Proposed solution = Build full size array from sparse?
                           = Bounds checking in the function?
       
*/

    /// For each output pixel, store the starting index in m_cost_buffer/m_accum_buffer
    ImageView<size_t> m_buffer_starts;

private: // Functions

  /// Populate the lookup table m_adjacent_disp_lookup
  void populate_adjacent_disp_lookup_table();

  /// Fill in m_disp_bound_image using image-wide contstants
  void populate_constant_disp_bound_image();

  /// Fill in m_disp_bound_image using some image information.
  /// - Returns false if there is no valid image data.
  /// - The left and right image masks contain a nonzero value if the pixel is valid.
  ///   No search is performed at masked pixels.
  /// - prev_disparity is a half resolution disparity image.
  bool populate_disp_bound_image(ImageView<uint8> const* left_image_mask,
                                 ImageView<uint8> const* right_image_mask,
                                 DisparityImage const* prev_disparity,
                                 int search_buffer);

  /// Fills m_buffer_starts and allocates m_cost_buffer and m_accum_buffer
  void allocate_large_buffers();
  
  /// Return a bad accumulation value used to fill locations we don't visit
  AccumCostType get_bad_accum_val() const { return std::numeric_limits<CostType>::max() + m_p2; }

  /// Returns the number of disparities searched for a given pixel.
  /// - This gets called a lot, may need to speed it up!
  int get_num_disparities(int col, int row) const {
    const Vector4i bounds = m_disp_bound_image(col,row);
    return (bounds[2] - bounds[0] + 1) * (bounds[3] - bounds[1] + 1);
  }

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
  //size_t get_cost_index(int col, int row, DisparityType disp=0) const {
  //    return row*m_buffer_step_size + col*m_num_disp + disp;
  //}

  /// Set provided buffer to the cost value at the selected pixel.
  //void set_cost_vector(int col, int row,
  //                     boost::shared_array<AccumCostType> accum_vec) {
  //  size_t start_index = get_cost_index(col, row, 0);
  //  for (int d=0; d<m_num_disp; ++d)
  //    accum_vec[start_index+d] = m_cost_buffer[start_index+d];
  //}
  
  CostType* get_cost_vector(int col, int row) {
    size_t start_index = m_buffer_starts(col, row);
    return m_cost_buffer.get() + start_index;
  };
  
  AccumCostType* get_accum_vector(int col, int row) {
    size_t start_index = m_buffer_starts(col, row);
    return m_accum_buffer.get() + start_index;
  };


  /// Generate the output disparity view from the accumulated costs.
  DisparityImage create_disparity_view();

/*
  // TODO: Use a lookup table?
  /// Converts from a linear disparity index to the dx, dy values it represents.
  /// - This function is too slow to use inside the inner loop!
  void disp_to_xy(DisparityType disp, DisparityType &dx, DisparityType &dy) {
    dy = (disp / m_num_disp_x) + m_min_disp_y; // 2D implementation
    dx = (disp % m_num_disp_x) + m_min_disp_x;
  }
*/
  /// Given the dx and dy positions of a pixel, return the 
  ///  full size disparity index.
  DisparityType xy_to_disp(DisparityType dx, DisparityType dy) {
    return (dy-m_min_disp_y)*m_num_disp_x + (dx-m_min_disp_x);
  }


  
  /// Get the value and index of the smallest element in a vector
  AccumCostType get_accum_vector_min(int col, int row,
                                     DisparityType &dx, DisparityType &dy){
  
    AccumCostType* vec = get_accum_vector(col, row);
    //CostType* vec = get_cost_vector(col, row); // DEBUG!!!
    const int num_disp = get_num_disparities(col, row);
    
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
  
/*
  // Print out a disparity vector
  template <typename T>
  void print_disparity_vector(T* const vec){
    std::cout << "V: ";
    for (int i=0; i<m_num_disp; ++i)
      std::cout << vec[i] << " ";
    std::cout << std::endl;
  }
*/
  
  /// Get the pixel diff along a line at a specified output location.
  int get_path_pixel_diff(ImageView<uint8> const& left_image,
                          int col, int row, int dir_x, int dir_y) const {
    // Take the offset between the output location and the input pixel coordinates.
    int a = left_image(col        +m_min_col, row        +m_min_row);
    int b = left_image((col-dir_x)+m_min_col, (row-dir_y)+m_min_row);
    return abs(a - b);
  }

  /// Create an updated cost accumulation vector for the next pixel along an SGM evaluation path.
  /// - For each disparity in the current pixel, add that disparity's cost with the "cheapest"
  ///   prior pixel disparity.
  /// - Returns the minimum disparity score
  void evaluate_path( int col, int row, int col_p, int row_p,
                      AccumCostType* const prior,      // Accumulated costs leading up to this pixel, truncated
                      AccumCostType*       full_prior_buffer, // Buffer to store all accumulated costs
                      CostType     * const local, // The disparity costs of the current pixel
                      AccumCostType*       output,
                      int path_intensity_gradient, bool debug=false ); // This variable is the magnitude of intensity change to this pixel

  /// Perform all eight path accumulations in two passes through the image
  void two_trip_path_accumulation(ImageView<uint8> const& left_image);


  /// Allow this helper class to access private members
  friend class MultiAccumRowBuffer;

}; // end class SemiGlobalMatcher


/// Wrapper function for SGM that handles ROIs.
/// - Merge with the function in Correlation.h!
/// - This function only searches positive disparities. The input images need to be
///   already cropped so that this makes sense.
/// - This function could be made more flexible by accepting other varieties of mask images!
template <class ImageT1, class ImageT2>
  ImageView<PixelMask<Vector2i> >
  calc_disparity_sgm(ImageViewBase<ImageT1> const& left_in,
                     ImageViewBase<ImageT2> const& right_in,
                     BBox2i                 const& left_region,   // Valid region in the left image
                     Vector2i               const& search_volume, // Max disparity to search in right image
                     Vector2i               const& kernel_size,  // Only really takes an N by N kernel!
                     ImageView<uint8>       const* left_mask_ptr=0,  
                     ImageView<uint8>       const* right_mask_ptr=0,
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

    Vector2i search_volume_inclusive = search_volume - Vector2i(1,1);

    // Rasterize input so that we can do a lot of processing on it.
    BBox2i right_region = left_region;
    right_region.max() += search_volume_inclusive;
    
    std::cout << "calc_disparity_sgm: left  region  = " << left_region   << std::endl;
    std::cout << "calc_disparity_sgm: right region  = " << right_region  << std::endl;
    std::cout << "calc_disparity_sgm: search_volume = " << search_volume_inclusive << std::endl;
    
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
    matcher.set_parameters(0, 0, search_volume_inclusive[0], search_volume_inclusive[1], kernel_size[0]);
    return matcher.semi_global_matching_func(left, right, left_mask_ptr, right_mask_ptr, prev_disparity);
    
  } // End function calc_disparity


} // end namespace stereo
} // end namespace vw

#endif //__SEMI_GLOBAL_MATCHING__
