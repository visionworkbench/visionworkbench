#ifndef __SEMI_GLOBAL_MATCHING_H__
#define __SEMI_GLOBAL_MATCHING_H__

#include <vw/Image/ImageView.h>
#include <vw/Image/PixelMask.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Stereo/Correlation.h>
#include <vw/Image/CensusTransform.h>
#include <vw/Image/Algorithms.h>

#include <boost/smart_ptr/shared_ptr.hpp>

//#include <vw/InterestPoint/Detector.h> // TODO: REMOVE THIS!

#if defined(VW_ENABLE_SSE) && (VW_ENABLE_SSE==1)
  #include <emmintrin.h>
  #include <smmintrin.h> // SSE4.1
#endif

namespace vw {

namespace stereo {

/**
A 2D implentation of the popular Semi-Global Matching (SGM) algorithm.  This 
implementation has the following features:
- 2D search using the passed in search range.
- Uses the popular Census cost function in a variable kernel size.
- Accepts an input low-resolution disparity image as a search seed.
  By using this input, the search range of each pixel is individually
  computed to minimize the run time.
- The large memory buffers required by the algorithm are compressed to contain
  only the individual search range for every pixel.  When combined with an
  input low-resolution disparity image, this can massively reduce the amount
  of memory required.
- SSE instructions are used to increase speed but currently they only provide
  a small improvement.
  
Even with the included optimizations this algorithm is slow and requires huge
amounts of memory to operate on large images.  Be careful not to exceed your
available memory when using it!

Future improvements:
- Implement an option in our pyramid correlation to short-circuit the lowest
  levels of the pyramid, enabling a fast computation of a low-resolution stereo output.
- Speed up the census cost function calculations.
- Optimize the algorithm parameters for our common use cases.
- Create a sub-pixel disparity step that can be used as an alternative
  to our existing sub-pixel algorithms.
- Experiment with filtering operations (such as described in the original SGM paper)
  that play well with the results of SGM.  This would be a way to avoid low-confidence
  detections
- Try to find algorithmic improvements.
- Try to further optimize the speed of the expensive accumulation step.
- Generate a confidence score for each pixel's disparity result.
- Make sure everything works with negative disparity search ranges.
*/

class SemiGlobalMatcher {

public: // Definitions

  // The types are chosen to minimize storage costs
  typedef int    DisparityType; ///< Contains the allowable dx, dy range.
  typedef uint8  CostType;      ///< Used to describe a single disparity cost.
  typedef uint16 AccumCostType; ///< Used to accumulate CostType values.
  
  typedef ImageView<PixelMask<Vector2i> > DisparityImage; // The usual VW disparity type

public: // Functions

  SemiGlobalMatcher() {} ///< Default constructor
  ~SemiGlobalMatcher() {} ///< Destructor

  /// Set set_parameters for details
  SemiGlobalMatcher(CostFunctionType cost_type,
                    bool use_mgm,
                    int min_disp_x, int min_disp_y,
                    int max_disp_x, int max_disp_y,
                    int kernel_size=5,
                    uint16 p1=0, uint16 p2=0,
                    int ternary_census_threshold=5) {
    set_parameters(cost_type, use_mgm, min_disp_x, min_disp_y, max_disp_x, max_disp_y, 
                   kernel_size, p1, p2);
  }

  /// Set the parameters to be used for future SGM calls
  /// - Parameters that are not provided will be set to the best known default.
  /// - If kernel_size is 3 or 5, a census transform will be used (reccomended).
  ///   Otherwise a simple averaging over a block method will be used.
  /// - p1 and p2 are algorithm constants very similar to those from the original SGM algorithm.
  ///   If not provided, they well be set to defaults according to the kernel size.
  void set_parameters(CostFunctionType cost_type,
                      bool use_mgm,
                      int min_disp_x, int min_disp_y,
                      int max_disp_x, int max_disp_y,
                      int kernel_size=5,
                      uint16 p1=0, uint16 p2=0,
                      int ternary_census_threshold=5);

  /// Compute SGM stereo on the images.
  /// The masks and disparity inputs are used to improve the searched disparity range.
  /// - The search buffer value is very important, it defines the radius around each
  ///   estimated disparity value that we will search.  A value of 2 means a 5x5 search 
  ///   region.  A larger region directly affects the speed and memory usage of SGM.
  DisparityImage
  semi_global_matching_func( ImageView<uint8> const& left_image,
                             ImageView<uint8> const& right_image,
                             ImageView<uint8> const* left_image_mask=0,
                             ImageView<uint8> const* right_image_mask=0,
                             DisparityImage const* prev_disparity=0,
                             int search_buffer = 4); // TODO: Restore this?

  /// Create a subpixel leves disparity image using parabola interpolation
  ImageView<PixelMask<Vector2f> > create_disparity_view_subpixel(DisparityImage const& integer_disparity);

private: // Variables

    // The core parameters
    int  m_min_disp_x, m_min_disp_y;
    int  m_max_disp_x, m_max_disp_y;
    int  m_kernel_size;              ///< Must be odd. Use "1" for single pixel.
    int  m_ternary_census_threshold; ///< Used only with the ternary census option
    bool m_use_mgm;
    CostFunctionType m_cost_type;

    int m_min_row, m_max_row;
    int m_min_col, m_max_col;
    int m_num_output_cols, m_num_output_rows;
    
    // Algorithm parameters
    AccumCostType m_p1;
    AccumCostType m_p2;
    
    // Derived parameters for convenience
    int m_num_disp_x, m_num_disp_y, m_num_disp;
    
    // The two main memory buffers that must be allocated.
    boost::shared_array<CostType     > m_cost_buffer;
    boost::shared_array<AccumCostType> m_accum_buffer;
    
    /// Image containing the inclusive disparity bounds for each pixel.
    /// - Stored as min_col, min_row, max_col, max_row.
    ImageView<Vector4i> m_disp_bound_image;
    
    /// Lookup table of the adjacent disparities for each disparity
    /// - For each disparity index, store the disparity indices of the 
    ///   eight adjacent disparities.
    /// - Handles outer boundaries by repetition.
    /// - This vector stores a table of size m_num_disp*8.
    /// - This allows us to avoid performing any bounds checking in the algorithm core.
    std::vector<DisparityType> m_adjacent_disp_lookup;

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
                                 DisparityImage   const* prev_disparity,
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
  // The following functions use two census transform cost function options
  void fill_costs_census3x3(ImageView<uint8> const& left_image, ImageView<uint8> const& right_image);
  void fill_costs_census5x5(ImageView<uint8> const& left_image, ImageView<uint8> const& right_image);
  void fill_costs_census7x7(ImageView<uint8> const& left_image, ImageView<uint8> const& right_image);
  void fill_costs_census9x9(ImageView<uint8> const& left_image, ImageView<uint8> const& right_image);

  /// Used to finish computing the census-based disparity costs in the above functions.
  template <typename T>
  void get_hamming_distance_costs(ImageView<T> const& left_binary_image,
                                  ImageView<T> const& right_binary_image);
                         
  /// Returns a block cost score at a given location
  CostType get_cost_block(ImageView<uint8> const& left_image,
                    ImageView<uint8> const& right_image,
                    int left_x, int left_y, int right_x, int right_y, bool debug);
  
  /// Get a pointer to a cost vector
  CostType* get_cost_vector(int col, int row) {
    size_t start_index = m_buffer_starts(col, row);
    return m_cost_buffer.get() + start_index;
  };
  
  /// Get a pointer to an accumulated cost vector
  AccumCostType* get_accum_vector(int col, int row) {
    size_t start_index = m_buffer_starts(col, row);
    return m_accum_buffer.get() + start_index;
  };

  /// Generate the output disparity view from the accumulated costs.
  DisparityImage create_disparity_view();


  /// Given the dx and dy positions of a pixel, return the 
  ///  full size disparity index.
  DisparityType xy_to_disp(DisparityType dx, DisparityType dy) {
    return (dy-m_min_disp_y)*m_num_disp_x + (dx-m_min_disp_x);
  }
  
  /// Get the value and index of the smallest element in an accumulation vector
  AccumCostType get_accum_vector_min(int col, int row,
                                     DisparityType &dx, DisparityType &dy);

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
                      AccumCostType* const prior,             // Accumulated costs leading up to this pixel, truncated
                      AccumCostType*       full_prior_buffer, // Buffer to store all accumulated costs
                      CostType     * const local,             // The disparity costs of the current pixel
                      AccumCostType*       output,
                      int path_intensity_gradient, bool debug=false ); // The magnitude of intensity change to this pixel

  /// Perform all eight path accumulations in two passes through the image
  void two_trip_path_accumulation(ImageView<uint8> const& left_image);
  
  /// Perform a smoother path accumulation using the MGM algorithm.
  /// - This method requires four passes and takes longer.
  void smooth_path_accumulation(ImageView<uint8> const& left_image);

  /// Allow this helper class to access private members.
  /// - This class can be found in SGM.cc.
  friend class MultiAccumRowBuffer;

  /// Converts from a linear disparity index to the dx, dy values it represents.
  /// - This function is too slow to use inside the inner loop!
  void disp_to_xy(DisparityType disp, DisparityType &dx, DisparityType &dy) {
    dy = (disp / m_num_disp_x) + m_min_disp_y; // 2D implementation
    dx = (disp % m_num_disp_x) + m_min_disp_x;
  }
  

  //// Print out a disparity vector
  //template <typename T>
  //void print_disparity_vector(T* const vec){
  //  std::cout << "V: ";
  //  for (int i=0; i<m_num_disp; ++i)
  //    std::cout << vec[i] << " ";
  //  std::cout << std::endl;
  //}

  // The following functions are inlined for speed
#if defined(VW_ENABLE_SSE) && (VW_ENABLE_SSE==1)
  /// Use SSE instructions to simultaneousy compute the scores for up to 8 disparities in evaluate_path()
  inline void compute_path_internals_sse(uint16* dL, uint16* d0, uint16* d1, uint16* d2, uint16* d3,
                                         uint16* d4, uint16* d5, uint16* d6, uint16* d7, uint16* d8,
                                         __m128i& _dJ, __m128i& _dP, __m128i& _dp1, uint16* dRes,
                                         int sse_index, int &output_index, AccumCostType* output);
#endif

  /// Non-sse backup for compute_path_internals_sse
  inline void compute_path_internals(uint16* dL, uint16* d0, uint16* d1, uint16* d2, uint16* d3,
                                     uint16* d4, uint16* d5, uint16* d6, uint16* d7, uint16* d8,
                                     AccumCostType dJ, AccumCostType dP, AccumCostType dp1, uint16* dRes,
                                     int sse_index, int &output_index, AccumCostType* output);


}; // end class SemiGlobalMatcher

/// Wrapper function for SGM that handles ROIs.
/// - This call is set up to be easily compatible with our existing calls in the
///   PyramidCorrelationView class.
/// - This function only searches positive disparities. The input images need to be
///   already cropped so that this makes sense.
/// - This function could be made more flexible by accepting other varieties of mask images.
/// - TODO: Merge with the function in Correlation.h?
template <class ImageT1, class ImageT2>
ImageView<PixelMask<Vector2i> >
calc_disparity_sgm(CostFunctionType cost_type,
                   ImageViewBase<ImageT1> const& left_in,
                   ImageViewBase<ImageT2> const& right_in,
                   BBox2i                 const& left_region,   // Valid region in the left image
                   Vector2i               const& search_volume, // Max disparity to search in right image
                   Vector2i               const& kernel_size,  // Only really takes an N by N kernel!
                   bool                   const  use_mgm,
                   boost::shared_ptr<SemiGlobalMatcher> &matcher_ptr,
                   ImageView<uint8>       const* left_mask_ptr=0,  
                   ImageView<uint8>       const* right_mask_ptr=0,
                   SemiGlobalMatcher::DisparityImage  const* prev_disparity=0);


//#################################################################################################
// Function definitions

#if defined(VW_ENABLE_SSE) && (VW_ENABLE_SSE==1)
void SemiGlobalMatcher::compute_path_internals_sse(uint16* dL, uint16* d0, uint16* d1, uint16* d2, uint16* d3,
                                                   uint16* d4, uint16* d5, uint16* d6, uint16* d7, uint16* d8,
                                                   __m128i& _dJ, __m128i& _dP, __m128i& _dp1, uint16* dRes,
                                                   int sse_index, int &output_index,
                                                   AccumCostType*       output) {
  // Load data from arrays into SSE registers
  __m128i _dL = _mm_load_si128( (__m128i*) dL );
  __m128i _d0 = _mm_load_si128( (__m128i*) d0 );
  __m128i _d1 = _mm_load_si128( (__m128i*) d1 );
  __m128i _d2 = _mm_load_si128( (__m128i*) d2 );
  __m128i _d3 = _mm_load_si128( (__m128i*) d3 );
  __m128i _d4 = _mm_load_si128( (__m128i*) d4 );
  __m128i _d5 = _mm_load_si128( (__m128i*) d5 );
  __m128i _d6 = _mm_load_si128( (__m128i*) d6 );
  __m128i _d7 = _mm_load_si128( (__m128i*) d7 );
  __m128i _d8 = _mm_load_si128( (__m128i*) d8 );

  // Operation = min( min(d1...d8)+dp1, d0, dJ) + dL - dP
  
  // Start computing the min
  __m128i _min12   = _mm_min_epu16(_d1, _d2);
  __m128i _min34   = _mm_min_epu16(_d3, _d4);
  __m128i _min56   = _mm_min_epu16(_d5, _d6);
  __m128i _min78   = _mm_min_epu16(_d7, _d8);
  // Keep computing the min
  __m128i _min1234 = _mm_min_epu16(_min12, _min34);
  __m128i _min5678 = _mm_min_epu16(_min56, _min78);
  // Finish computing the min
  __m128i _minAdj = _mm_min_epu16(_min1234, _min5678);
  __m128i _minO   = _mm_min_epu16(_d0, _dJ);
 
  // Perform the required computations
  __m128i _result = _mm_adds_epu16(_minAdj, _dp1);
  _result = _mm_min_epu16(_result, _minO);
  _result = _mm_adds_epu16(_result, _dL);
  _result = _mm_subs_epu16(_result, _dP);

  // Fetch results from the output register
  _mm_store_si128( (__m128i*) dRes, _result );

  // Copy the valid results from the register.
  for (int i=0; i<sse_index; ++i){
    output[output_index++] = dRes[i];
  }
} // end function compute_path_internals_sse
#endif

void SemiGlobalMatcher::compute_path_internals(uint16* dL, uint16* d0, uint16* d1, uint16* d2, uint16* d3,
                                               uint16* d4, uint16* d5, uint16* d6, uint16* d7, uint16* d8,
                                               AccumCostType dJ, AccumCostType dP, AccumCostType dp1, uint16* dRes,
                                               int sse_index, int &output_index,
                                               AccumCostType*       output) {
  //std::cout << "dL[i], d0[i], d1[i], d2[i], d3[i], d4[i], d5[i], d6[i], d7[i], d8[i]\n";      
  //for (int i=0; i<sse_index; ++i) {
  //  printf("%d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n", 
  //      dL[i], d0[i], d1[i], d2[i], d3[i], d4[i], d5[i], d6[i], d7[i], d8[i]);
  //}

  // Operation = min( min(d1...d8)+dp1, d0, dJ) + dL - dP

  for (int i=0; i<sse_index; ++i){
  
    uint16 minAdj = std::min(d1[i], d2[i]);
    minAdj = std::min(minAdj, d3[i]);
    minAdj = std::min(minAdj, d4[i]);
    minAdj = std::min(minAdj, d5[i]);
    minAdj = std::min(minAdj, d6[i]);
    minAdj = std::min(minAdj, d7[i]);
    minAdj = std::min(minAdj, d8[i]);
    minAdj += dp1;
    
    uint16 minVal = std::min(minAdj, d0[i]);
           minVal = std::min(minVal, dJ);
    output[output_index++] = minVal + (dL[i] - dP);
  }
} // end function compute_path_internals


// From the census transformed input images, compute the cost of each disparity value.
template <typename T>
void SemiGlobalMatcher::get_hamming_distance_costs(ImageView<T> const& left_binary_image,
                                                   ImageView<T> const& right_binary_image) {

  const int half_kernel = (m_kernel_size - 1) / 2;

  // Now compute the disparity costs for each pixel.
  // Make sure we don't go out of bounds here due to the disparity shift and kernel.
  size_t cost_index = 0;
  for ( int r = m_min_row; r <= m_max_row; r++ ) { // For each row in left
    int output_row = r - m_min_row;
    int binary_row = r - half_kernel;
    for ( int c = m_min_col; c <= m_max_col; c++ ) { // For each column in left
      int output_col = c - m_min_col;
      int binary_col = c - half_kernel;
      
      Vector4i pixel_disp_bounds = m_disp_bound_image(output_col, output_row);
    
      for ( int dy = pixel_disp_bounds[1]; dy <= pixel_disp_bounds[3]; dy++ ) { // For each disparity
        for ( int dx = pixel_disp_bounds[0]; dx <= pixel_disp_bounds[2]; dx++ ) {
          
          CostType cost = hamming_distance(left_binary_image (binary_col   , binary_row   ), 
                                           right_binary_image(binary_col+dx, binary_row+dy) );
          m_cost_buffer[cost_index] = cost;
          ++cost_index;
        }    
      } // End disparity loops   
    } // End x loop
  }// End y loop 
                      
}


//TODO: Move this function!
/// Converts a single channel image into a uint8 image with percentile based intensity scaling.
template <class ViewT>
void u8_convert(ImageViewBase<ViewT> const& input_image, ImageView<PixelGray<vw::uint8> > &output_image, int num_bins=256) {
  // First get the min and max values
  double min_val, max_val;
  find_image_min_max(input_image, min_val, max_val);

  // Scale the image using the computed values and convert to uint8
  output_image = pixel_cast<vw::uint8>(normalize( clamp(input_image, min_val, max_val),
					    min_val, max_val, 0.0, 255.0 ));
}



template <class ImageT1, class ImageT2>
ImageView<PixelMask<Vector2i> >
calc_disparity_sgm(CostFunctionType cost_type,
                   ImageViewBase<ImageT1> const& left_in,
                   ImageViewBase<ImageT2> const& right_in,
                   BBox2i                 const& left_region,   // Valid region in the left image
                   Vector2i               const& search_volume, // Max disparity to search in right image
                   Vector2i               const& kernel_size,  // Only really takes an N by N kernel!
                   bool                   const  use_mgm,
                   boost::shared_ptr<SemiGlobalMatcher> &matcher_ptr,
                   ImageView<uint8>       const* left_mask_ptr,  
                   ImageView<uint8>       const* right_mask_ptr,
                   SemiGlobalMatcher::DisparityImage  const* prev_disparity){ 

    
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
    
    vw_out(VerboseDebugMessage, "stereo") << "calc_disparity_sgm: left  region  = " << left_region   << std::endl;
    vw_out(VerboseDebugMessage, "stereo") << "calc_disparity_sgm: right region  = " << right_region  << std::endl;
    vw_out(VerboseDebugMessage, "stereo") << "calc_disparity_sgm: search_volume_inclusive = " << search_volume_inclusive << std::endl;
    
    // TODO: Ignore masked values when computing this!
    // Convert the input image to uint8
    ImageView<PixelGray<vw::uint8> > left, right;
    //ip::percentile_scale_convert(crop(left_in.impl(),  left_region),  left,  0.00, 1.00); // Any stretching seems to cause problems!
    //ip::percentile_scale_convert(crop(right_in.impl(), right_region), right, 0.00, 1.00);
    u8_convert(crop(left_in.impl(),  left_region),  left);
    u8_convert(crop(right_in.impl(), right_region), right);    
    
    //write_image("final_left.tif", left);
    //write_image("final_right.tif", right);
    
    matcher_ptr.reset(new SemiGlobalMatcher(cost_type, use_mgm, 0, 0, 
                      search_volume_inclusive[0], search_volume_inclusive[1], kernel_size[0]));
    return matcher_ptr->semi_global_matching_func(left, right, left_mask_ptr, right_mask_ptr, prev_disparity);
    
  } // End function calc_disparity


} // end namespace stereo
} // end namespace vw

#endif //__SEMI_GLOBAL_MATCHING__
