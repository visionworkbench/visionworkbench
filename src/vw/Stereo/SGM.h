// __BEGIN_LICENSE__
//  Copyright (c) 2009-2025, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NGT platform is licensed under the Apache License, Version 2.0 (the
//  "License"); you may not use this file except in compliance with the
//  License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__

#ifndef __SEMI_GLOBAL_MATCHING_H__
#define __SEMI_GLOBAL_MATCHING_H__

#include <vw/Image/ImageView.h>
#include <vw/Image/PixelMask.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Image/CensusTransform.h>
#include <vw/Image/Algorithms.h>
#include <vw/Stereo/CostFunctions.h>

#include <boost/smart_ptr/shared_ptr.hpp>

#if defined(VW_ENABLE_SSE) && (VW_ENABLE_SSE==1)
  #include <emmintrin.h>
  #include <smmintrin.h> // SSE4.1
#endif

namespace vw {

namespace stereo {

/**
A 2D implementation of the popular Semi-Global Matching (SGM) algorithm. This
implementation has the following features:
- 2D search using the passed in search range.
- Uses the popular Census cost function in a variable kernel size.
- Accepts an input low-resolution disparity image as a search seed.
  By using this input, the search range of each pixel is individually
  computed to minimize the run time.
- The large memory buffers required by the algorithm are compressed to contain
  only the individual search range for every pixel. When combined with an
  input low-resolution disparity image, this can massively reduce the amount
  of memory required.
- SSE instructions are used to increase speed but currently they only provide
  a small improvement.

Even with the included optimizations this algorithm is slow and requires huge
amounts of memory to operate on large images. Be careful not to exceed your
available memory when using it!

Future improvements:
- Implement an option in our pyramid correlation to short-circuit the lowest
  levels of the pyramid, enabling a fast computation of a low-resolution stereo output.
- Speed up the census cost function calculations.
- Optimize the algorithm parameters for our common use cases.
- Create a sub-pixel disparity step that can be used as an alternative
  to our existing sub-pixel algorithms.
- Experiment with filtering operations (such as described in the original SGM paper)
  which play well with the results of SGM. This would be a way to avoid low-confidence
  detection.
- Try to find algorithmic improvements.
- Try to further optimize the speed of the expensive accumulation step.
- Generate a confidence score for each pixel's disparity result.
- Make sure everything works with negative disparity search ranges. This never
  comes up when called from CorrelationView, but would make the class more flexible.
*/

class SemiGlobalMatcher {

public: // Definitions

  // The types are chosen to minimize storage costs
  typedef int    DisparityType; ///< Contains the allowable dx, dy range.
  typedef uint8  CostType;      ///< Used to describe a single disparity cost.
  typedef uint16 AccumCostType; ///< Used to accumulate CostType values.

  // For converting buffer sizes to MB
  static constexpr size_t BYTES_PER_MB = 1024 * 1024;
  static constexpr double MainBufToMB
    = double(sizeof(CostType) + sizeof(AccumCostType)) / BYTES_PER_MB;
  static constexpr double SmallBufToMB = double(sizeof(AccumCostType)) / BYTES_PER_MB;

  typedef ImageView<PixelMask<Vector2i>> DisparityImage; // The usual VW disparity type

  /// Available subpixel options
  enum SgmSubpixelMode {SUBPIXEL_NONE     = 0, // Skip subpixel processing
                        SUBPIXEL_PARABOLA = 1, // The only combined 2D option
                        SUBPIXEL_LINEAR   = 2,
                        SUBPIXEL_POLY4    = 3,
                        SUBPIXEL_COSINE   = 4,
                        SUBPIXEL_LC_BLEND = 5  // Probably the best option
                        };

public: // Functions

  SemiGlobalMatcher() {} ///< Default constructor
  ~SemiGlobalMatcher() {} ///< Destructor

  /// Set set_parameters for details
  SemiGlobalMatcher(CostFunctionType cost_type,
                    bool use_mgm,
                    int min_disp_x, int min_disp_y,
                    int max_disp_x, int max_disp_y,
                    int kernel_size=5,
                    SgmSubpixelMode subpixel=SUBPIXEL_LC_BLEND,
                    Vector2i search_buffer=Vector2i(2,2),
                    size_t memory_limit_mb=6000,
                    uint16 p1=0, uint16 p2=0,
                    int ternary_census_threshold=5) {
    set_parameters(cost_type, use_mgm, min_disp_x, min_disp_y, max_disp_x, max_disp_y,
                   kernel_size, subpixel, search_buffer, memory_limit_mb,
                   p1, p2, ternary_census_threshold);
  }

  /// Set the parameters to be used for future SGM calls
  /// - Parameters that are not provided will be set to the best known default.
  /// - p1 and p2 are algorithm constants very similar to those from the original SGM algorithm.
  ///   If not provided, they well be set to defaults according to the kernel size and cost type.
  /// - The search buffer value is important, it defines the radius around each
  ///   estimated disparity value that we will search.  A value of 2 means a 5x5 search
  ///   region.  A larger region directly affects the speed and memory usage of SGM.
  /// - memory_limit_mb is the maximum amount of memory that the algorithm is allowed to allocate
  ///   for its large buffers (total memory usage can go slightly over this).  The program will
  ///   attempt a more conservative search range if needed to get under this target and if that
  ///   fails then the program will throw an exception.
  void set_parameters(CostFunctionType cost_type,
                      bool use_mgm,
                      int min_disp_x, int min_disp_y,
                      int max_disp_x, int max_disp_y,
                      int kernel_size=5,
                      SgmSubpixelMode subpixel=SUBPIXEL_LC_BLEND,
                      Vector2i search_buffer=Vector2i(2,2),
                      size_t memory_limit_mb=6000,
                      uint16 p1=0, uint16 p2=0,
                      int ternary_census_threshold=5);

  /// Compute SGM stereo on the images.
  /// The masks and disparity inputs are used to improve the searched disparity range.
  DisparityImage
  semi_global_matching_func(ImageView<uint8> const& left_image,
                            ImageView<uint8> const& right_image,
                            ImageView<uint8> const* left_image_mask,
                            ImageView<uint8> const* right_image_mask,
                            DisparityImage const* prev_disparity);

  /// Create a subpixel disparity image using parabola interpolation
  ImageView<PixelMask<Vector2f> > create_disparity_view_subpixel(DisparityImage const& integer_disparity);

private: // Variables

    // The core parameters
    int  m_min_disp_x, m_min_disp_y;
    int  m_max_disp_x, m_max_disp_y;
    int  m_kernel_size;              ///< Must be odd. Use "1" for single pixel.
    int  m_ternary_census_threshold; ///< Used only with the ternary census option
    bool m_use_mgm;
    CostFunctionType m_cost_type;
    SgmSubpixelMode  m_subpixel_type;
    Vector2i m_search_buffer;
    size_t m_memory_limit_mb; ///< Maximum memory usage allowed in main buffers

    int m_min_row, m_max_row;
    int m_min_col, m_max_col;
    int m_num_output_cols, m_num_output_rows;

    // Algorithm parameters
    AccumCostType m_p1;
    AccumCostType m_p2;

    // Derived parameters for convenience
    int m_num_disp_x, m_num_disp_y, m_num_disp;

    // The two main memory buffers that must be allocated.
    boost::shared_array<CostType> m_cost_buffer;
    boost::shared_array<AccumCostType> m_accum_buffer;
    size_t m_main_buf_size;

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

  /// Fill in m_disp_bound_image using image-wide constants
  void populate_constant_disp_bound_image();

  /// Fill in m_disp_bound_image using some image information.
  /// - Returns false if there is no valid image data.
  /// - The left and right image masks contain a nonzero value if the pixel is valid.
  ///   No search is performed at masked pixels.
  /// - prev_disparity is a half resolution disparity image.  This input is optional,
  ///   pass in a null pointer to ignore it.  If provided, it will be used to limit
  ///   the disparity range searched at each pixel.
  /// - The range of disparities searched controls the run time and memory usage of the
  ///   algorithm so this is an important function!
  bool populate_disp_bound_image(ImageView<uint8> const* left_image_mask,
                                 ImageView<uint8> const* right_image_mask,
                                 DisparityImage   const* prev_disparity);

  /// Reduce the search range of full-search-range pixel by looking at nearby
  ///  pixels with a smaller search range.
  /// - conserve_memory controls how aggressive the function is in finding possible
  ///   search areas for uncertain pixels.  Level 0 uses full search ranges.
  /// - Returns false if there is a problem processing the image.
  bool constrain_disp_bound_image(ImageView<uint8> const &full_search_image,
                                  DisparityImage const* prev_disparity,
                                  double percent_trusted, double percent_masked, double area,
                                  int conserve_memory=0);

  /// Return the number of elements in each of the large buffers.
  /// - Also perform a check to make sure our memory usage falls within the user specified limit.
  size_t calc_main_buf_size();

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

  /// Compute the mean and STD of a small image patch.
  /// - Does not perform bounds checking.
  void compute_patch_mean_std(ImageView<uint8> const& image, int x, int y,
                              double &mean, double &std) const;

  /// Returns a block cost score at a given location.
  /// - Supports non-census cost modes.
  CostType get_cost_block(ImageView<uint8> const& left_image,
                    ImageView<uint8> const& right_image,
                    double mean_left, double std_left,
                    int left_x, int left_y, int right_x, int right_y, bool debug) const;

  /// Get a pointer to a cost vector
  CostType * get_cost_vector(int col, int row);

  /// Get a pointer to an accumulated cost vector
  AccumCostType* get_accum_vector(int col, int row);

  /// Generate the output disparity view from the accumulated costs.
  DisparityImage create_disparity_view();

  /// Select the best disparity index in the accumulation vector.
  /// - If needed, applies smoothing to the values in order to yield a single minimum value.
  int select_best_disparity(AccumCostType * accum_vec,
                            Vector4i const& bounds,
                            int &min_index,
                            std::vector<AccumCostType> & buffer,
                            bool debug);

  /// Get the pixel diff along a line at a specified output location.
  int get_path_pixel_diff(ImageView<uint8> const& left_image,
                          int col, int row, int dir_x, int dir_y) const;

  /// Create an updated cost accumulation vector for the next pixel along an SGM evaluation path.
  /// - For each disparity in the current pixel, add that disparity's cost with the "cheapest"
  ///   prior pixel disparity.
  /// - Returns the minimum disparity score
  void evaluate_path(int col, int row, int col_p, int row_p,
                      AccumCostType* const prior,             // Accumulated costs leading up to this pixel, truncated
                      AccumCostType*       full_prior_buffer, // Buffer to store all accumulated costs
                      CostType     * const local,             // The disparity costs of the current pixel
                      AccumCostType*       output,
                      int path_intensity_gradient, bool debug=false); // The magnitude of intensity change to this pixel

  /// Perform all eight path accumulations in two passes through the image
  void two_trip_path_accumulation(ImageView<uint8> const& left_image);

  /// Perform a smoother path accumulation using the MGM algorithm.
  /// - This method requires four passes and takes longer.
  void smooth_path_accumulation(ImageView<uint8> const& left_image);

  void accum_mgm_multithread(ImageView<uint8> const& left_image);

  /// Multi-threaded version of the normal SGM accumulation method;
  void accum_sgm_multithread(ImageView<uint8> const& left_image);

  /// Allow this helper class to access private members.
  /// - Some of these classes can be found in SGMAssist.h.
  friend class MultiAccumRowBuffer;
  friend class SingleAccumRowBuffer;
  friend class OneLineBuffer;
  friend class PixelPassTask;
  friend class SmoothPathAccumTask;

  /// Given the dx and dy positions of a pixel, return the full size disparity
  /// index. - Note that the big storage vectors do not store the entire
  /// disparity range for each pixel.
  DisparityType xy_to_disp(DisparityType dx, DisparityType dy) const;

  /// Converts from a linear disparity index to the dx, dy values it represents.
  /// - This function is too slow to use inside the inner loop!
  void disp_to_xy(DisparityType disp, DisparityType &dx, DisparityType &dy) const;

  /// Convert a pixel's minimum disparity index to dx, dy.
  void disp_index_to_xy(int min_index, int col, int row,
                        DisparityType &dx, DisparityType &dy) const;

  /// Given disparity cost and adjacent costs, compute subpixel offset.
  double compute_subpixel_offset(AccumCostType prev, AccumCostType center, AccumCostType next,
                                 bool left_bound = false, bool right_bound = false,
                                 bool debug = false);

  /// Crude two-element subpixel estimation.
  double two_value_subpixel(AccumCostType primary, AccumCostType other);

  ///// Function to help with developing subpixel functions
  //double compute_subpixel_ratio(AccumCostType prev, AccumCostType center, AccumCostType next);

}; // end class SemiGlobalMatcher

/// Wrapper function for SGM that handles ROIs.
/// - This call is set up to be easily compatible with our existing calls in the
///   PyramidCorrelationView class.
/// - This function only searches positive disparities. The input images need to be
///   already cropped so that this makes sense.
/// - This function could be made more flexible by accepting other varieties of mask images.
ImageView<PixelMask<Vector2i>>
calc_disparity_sgm(
  CostFunctionType cost_type,
  ImageView<PixelGray<float>> const& left_in,
  ImageView<PixelGray<float>> const& right_in,
  BBox2i                 const& left_region,   // Valid region in the left image
  Vector2i               const& search_volume, // Max disparity to search in right image
  Vector2i               const& kernel_size,   // The kernel dimensions are always equal
  bool                   const  use_mgm,
  SemiGlobalMatcher::SgmSubpixelMode const& subpixel_mode,
  Vector2i               const  search_buffer, // Search buffer applied around prev_disparity
  size_t                 const memory_limit_mb,
  boost::shared_ptr<SemiGlobalMatcher> &matcher_ptr,
  ImageView<uint8>       const* left_mask_ptr=0,
  ImageView<uint8>       const* right_mask_ptr=0,
  SemiGlobalMatcher::DisparityImage const* prev_disparity=0);

} // end namespace stereo
} // end namespace vw

#endif //__SEMI_GLOBAL_MATCHING__
