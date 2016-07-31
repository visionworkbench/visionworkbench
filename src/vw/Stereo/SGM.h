#ifndef __SEMI_GLOBAL_MATCHING_H__
#define __SEMI_GLBOAL_MATCHING_H__

#include <vw/Image/ImageView.h>
#include <vw/Math/Vector.h>
#include <vw/FileIO.h>

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
  typedef int16 CostType;      ///< Used to describe a single disparity cost.
  typedef int32 AccumCostType; ///< Used to accumulate CostType values.

  typedef ImageView<Vector<DisparityType,2> > DisparityImage;

public: // Functions

  // Set the parameters to be used for future SGM calls
  void setParameters(int min_disp_x, int min_disp_y,
                     int max_disp_x, int max_disp_y,
                     int kernel_size);

  /// Invokes a 8 path version of SGM
  DisparityImage
  semi_global_matching_func( ImageView<uint8> const& left_image,
                             ImageView<uint8> const& right_image );

private: // Variables

    // The core parameters
    int m_min_disp_x, m_min_disp_y;
    int m_max_disp_x, m_max_disp_y;
    int m_kernel_size; ///< Must be odd. Use "1" for single pixel.

    int m_num_output_cols, m_num_output_rows;
    
    // Derived parameters for convenience
    int m_num_disp_x, m_num_disp_y, m_num_disp;
    
    size_t m_buffer_step_size;
    boost::shared_array<CostType     > m_cost_buffer;
    boost::shared_array<AccumCostType> m_accum_buffer;
    

private: // Functions

  /// Return the index into one of the buffers for a given location
  /// - The data is stored row major interleaved format.
  size_t get_cost_index(int col, int row, DisparityType disp=0) const {
      return row*m_buffer_step_size + col*m_num_disp + disp;
  }

  /// Add the value of the cost buffer to the accumulated cost vector at a pixel.
  void add_cost_vector(int col, int row,
                       boost::shared_array<AccumCostType> accum_vec) {
    size_t start_index = get_cost_index(col, row, 0);
    for (int d=0; d<m_num_disp; ++d)
      accum_vec[start_index+d] += m_cost_buffer[start_index+d];
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

  /// Converts from a linear disparity index to the dx, dy values it represents.
  void disp_to_xy(DisparityType disp, DisparityType &dx, DisparityType &dy) {
    dy = (disp / m_num_disp_x) + m_min_disp_y; // 2D implementation
    dx = (disp % m_num_disp_x) + m_min_disp_x;
  }

  /// Compute the physical distance between the disparity values of two adjacent pixels.
  AccumCostType get_disparity_dist(DisparityType d1, DisparityType d2);
  
  
  
  
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

  /// Returns a cost score at a given location
  CostType get_cost(ImageView<uint8> const& left_image,
                 ImageView<uint8> const& right_image,
                 int left_x, int left_y, int right_x, int right_y, bool debug);


  // Print out a disparity vector
  template <typename T>
  void print_disparity_vector(T* const vec){
    std::cout << "V: ";
    for (int i=0; i<m_num_disp; ++i)
      std::cout << vec[i] << " ";
    std::cout << std::endl;
  }
  
  

  /// Create an updated cost accumulation vector for the next pixel along an SGM evaluation path.
  /// - For each disparity in the current pixel, add that disparity's cost with the "cheapest"
  ///   prior pixel disparity.
  void evaluate_path( AccumCostType* const prior, // Accumulated costs leading up to this pixel
                      CostType     * const local, // The disparity costs of the current pixel
                      AccumCostType*       output,
                      int path_intensity_gradient, bool debug=false ); // This variable is the magnitude of intensity change to this pixel

  /// Compute the accumulated costs in a pixel direction from the local costs at each pixel.
  /// - TODO: This implementation seems inefficient!
  template <int DIRX, int DIRY>
  void iterate_direction( ImageView<uint8   > const& left_image,
                          boost::shared_array<AccumCostType>      & accumulated_costs ) {

    // Zero out the output data   
    size_t num_cost_elements = m_num_output_cols*m_num_output_rows*m_num_disp;
    for (size_t i=0; i<num_cost_elements; ++i)
      accumulated_costs[i] = 0;
    

    // Walk along the edges in a clockwise fashion
    if ( DIRX > 0 ) {
      // LEFT MOST EDGE
      // Init the edge pixels with just the cost (no accumulation yet)
      for ( int32 j = 0; j < m_num_output_rows; j++ ) {
        add_cost_vector(0, j, accumulated_costs);
      }
      //std::cout << "L costs: " << costs(0,185) << std::endl;

      // Loop across to the opposite edge
      for ( int32 i = 1; i < m_num_output_cols; i++ ) {
        // Loop through the pixels in this column, limiting the range according
        //  to the iteration direction progress.
        int32 jstart = std::max( 0,      0      + DIRY * i );
        int32 jstop  = std::min( m_num_output_rows, m_num_output_rows + DIRY * i );
        for ( int32 j = jstart; j < jstop; j++ ) {
          int pixel_diff = abs(left_image(i,j)-left_image(i-DIRX,j-DIRY));
          evaluate_path( get_accum_vector(accumulated_costs, i-DIRX,j-DIRY), // Previous pixel
                         get_cost_vector(i,j),
                         get_accum_vector(accumulated_costs, i,j), 
                         pixel_diff, false );         // Current pixel
                                                  
          //if (j == 185) {
          //  std::cout << "i: "<<i<< " costs: " << costs(i,j) << "\n      accum: " << accumulated_costs(i-DIRX,j-DIRY)  
          //                                                   << "\n      accum: " << accumulated_costs(i,j) << std::endl;
          //}
                                                  
        }
      }
    } 
    if ( DIRY > 0 ) {
      // TOP MOST EDGE
      // Process every pixel along this edge only if DIRX == 0. Otherwise skip the top left most pixel
      for ( int32 i = (DIRX <= 0 ? 0 : 1 ); i < m_num_output_cols; i++ ) {
        add_cost_vector(i, 0, accumulated_costs);
      }
      for ( int32 j = 1; j < m_num_output_rows; j++ ) {
        int32 istart = std::max( (DIRX <= 0 ? 0 : 1),
                                 (DIRX <= 0 ? 0 : 1) + DIRX * j );
        int32 istop  = std::min( m_num_output_cols, m_num_output_cols + DIRX * j );
        for ( int32 i = istart; i < istop; i++ ) {
          int pixel_diff         = abs(left_image(i,j)-left_image(i-DIRX,j-DIRY));
          evaluate_path( get_accum_vector(accumulated_costs, i-DIRX,j-DIRY), // Previous pixel
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
        add_cost_vector(m_num_output_cols-1, j, accumulated_costs);
      }
      for ( int32 i = m_num_output_cols-2; i >= 0; i-- ) {
        int32 jstart = std::max( (DIRY <= 0 ? 0 : 1),
                                 (DIRY <= 0 ? 0 : 1) - DIRY * (i - m_num_output_cols + 1) );
        int32 jstop  = std::min( m_num_output_rows, m_num_output_rows - DIRY * (i - m_num_output_cols + 1) );
        for ( int32 j = jstart; j < jstop; j++ ) {
          int pixel_diff         = abs(left_image(i,j)-left_image(i-DIRX,j-DIRY));
          evaluate_path( get_accum_vector(accumulated_costs, i-DIRX,j-DIRY), // Previous pixel
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
        add_cost_vector(i,m_num_output_rows-1, accumulated_costs);
      }
      for ( int32 j = m_num_output_rows-2; j >= 0; j-- ) {
        int32 istart = std::max( (DIRX <= 0 ? 0 : 1),
                                 (DIRX <= 0 ? 0 : 1) - DIRX * (j - m_num_output_rows + 1) );
        int32 istop  = std::min( (DIRX >= 0 ? m_num_output_cols : m_num_output_cols - 1),
                                 (DIRX >= 0 ? m_num_output_cols : m_num_output_cols - 1) - DIRX * (j - m_num_output_rows + 1) );
        for ( int32 i = istart; i < istop; i++ ) {
          int pixel_diff         = abs(left_image(i,j)-left_image(i-DIRX,j-DIRY));
          evaluate_path( get_accum_vector(accumulated_costs, i-DIRX,j-DIRY), // Previous pixel
                         get_cost_vector(i,j), 
                         get_accum_vector(accumulated_costs, i,j),
                         pixel_diff );         // Current pixel
        }
      }
    }
  } // End function iterate_direction



}; // end class SemiGlobalMatcher


/*

  /// SGM view class.
  /// - Currently only accepts uint8 images!
  /// - Do we really need a view for this?
  template <class Image1T, class Image2T>
  class SemiGlobalMatchingView : public ImageViewBase<SemiGlobalMatchingView<Image1T,Image2T> > {
    Image1T m_left_image;
    Image2T m_right_image;
  public:
    typedef uint8 pixel_type;
    typedef uint8 result_type;
    typedef ProceduralPixelAccessor<SemiGlobalMatchingView> pixel_accessor;

    SemiGlobalMatchingView( ImageViewBase<Image1T> const& left,
                            ImageViewBase<Image2T> const& right ) :
      m_left_image(left.impl()), m_right_image(right.impl()) {}

    inline int32 cols  () const { return m_left_image.cols(); }
    inline int32 rows  () const { return m_left_image.rows(); }
    inline int32 planes() const { return 1; }

    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }
    inline pixel_type operator()( int32 i, int32 j, int32 p = 0) const {
      vw_throw( NoImplErr() << "CorrelationView::operator()(....) has not been implemented." );
      return pixel_type();
    }

    // Block rasterization section that does actual work
    typedef CropView<ImageView<pixel_type> > prerasterize_type;
    inline prerasterize_type prerasterize(BBox2i const& bbox) const {
      // Rasterize the left image in the desired bbox
      ImageView<PixelGray<uint8> > left = crop( edge_extend(m_left_image), bbox );
      
      // Figure out the associated bbox in the right image
      // - This is the maximum possible match size.
      // - TODO: If we switch to block matching, need to add the width of the matching kernel.
      BBox2i rbbox = bbox;
      rbbox.max() += Vector2i(DISP_RANGE_X,DISP_RANGE_Y);
      
      // Rasterize the needed section of the right image
      ImageView<PixelGray<uint8> > right = crop( edge_extend(m_right_image), rbbox );
      
      // Call the SGM function, then do the crop trick to fake the whole resolution image.
      return prerasterize_type( semi_global_matching_func( left, right ),
                                -bbox.min().x(), -bbox.min().y(), cols(), rows() );
    }

    template <class DestT>
    inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
      vw::rasterize(prerasterize(bbox), dest, bbox);
    }
  };

  template <class Image1T, class Image2T>
  SemiGlobalMatchingView<Image1T,Image2T>
  semi_global_matching( ImageViewBase<Image1T> const& left,
                        ImageViewBase<Image2T> const& right ) {
    typedef SemiGlobalMatchingView<Image1T,Image2T> result_type;
    return result_type( left.impl(), right.impl() );
  }*/

} // end namespace stereo
} // end namespace vw

#endif //__SEMI_GLOBAL_MATCHING__
