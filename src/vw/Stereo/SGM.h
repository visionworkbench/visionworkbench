#ifndef __SEMI_GLOBAL_MATCHING_H__
#define __SEMI_GLBOAL_MATCHING_H__

#include <vw/Image/ImageView.h>
#include <vw/Math/Vector.h>
#include <vw/FileIO.h>

namespace vw {


  // Registering the Pixel Disparity type for FileIO
  template<> struct PixelFormatID<Vector<uint8,2>     > { static const PixelFormatEnum value = VW_PIXEL_GENERIC_2_CHANNEL; };



namespace stereo {

// TODO: Move everything into a namespace!

// TODO: Make into parameters!
#define DISP_RANGE_X 20
#define DISP_RANGE_Y 3
#define NUM_DISPS    60
#define SGM_KERNEL_SIZE 1

  // Cost Vector type, one cost for each disparity.
  typedef Vector<int16,NUM_DISPS> CVector;
  // Accumulation Vector type, used to accumulate CVector values.
  typedef Vector<int32,NUM_DISPS> AVector;




  /// Converts from a linear disparity index to the dx, dy values it represents.
  void disp_to_xy(int32 disp, int32 &dx, int32 &dy) {
    dy = disp / DISP_RANGE_X; // 2D implementation
    dx = disp % DISP_RANGE_X;
    //dx = disp; // This is the regular, 1D implementation!
    //dy = 0;
  }

  /// Compute the physical distance between the disparity values of two adjacent pixels.
  int32 get_disparity_dist(int32 d1, int32 d2);
  

  /// Create an updated cost accumulation vector for the next pixel along an SGM evaluation path.
  /// - For each disparity in the current pixel, add that disparity's cost with the "cheapest"
  ///   prior pixel disparity.
  AVector evaluate_path( AVector const& prior, // Accumulated costs leading up to this pixel
                         CVector const& local, // The disparity costs of the current pixel
                         int path_intensity_gradient, bool debug=false ); // This variable is the magnitude of intensity change to this pixel

  /// Compute the accumulated costs in a pixel direction from the local costs at each pixel.
  /// - TODO: This implementation seems inefficient!
  template <int DIRX, int DIRY>
  void iterate_direction( ImageView<uint8   > const& left_image,
                          ImageView<CVector > const& costs, // Costs for each disparity at each pixel
                          ImageView<AVector >      & accumulated_costs ) {
    const int32 WIDTH  = costs.cols();
    const int32 HEIGHT = costs.rows();

    // Zero out the output data
    std::fill(accumulated_costs.data(), accumulated_costs.data()+HEIGHT*WIDTH, AVector());

    // Walk along the edges in a clockwise fashion
    if ( DIRX > 0 ) {
      // LEFT MOST EDGE
      // Init the edge pixels with just the cost (no accumulation yet)
      for ( int32 j = 0; j < HEIGHT; j++ ) {
        accumulated_costs(0,j) += costs(0,j);
      }
      //std::cout << "L costs: " << costs(0,185) << std::endl;

      // Loop across to the opposite edge
      for ( int32 i = 1; i < WIDTH; i++ ) {
        // Loop through the pixels in this column, limiting the range according
        //  to the iteration direction progress.
        int32 jstart = std::max( 0,      0      + DIRY * i );
        int32 jstop  = std::min( HEIGHT, HEIGHT + DIRY * i );
        for ( int32 j = jstart; j < jstop; j++ ) {
          int pixel_diff         = abs(left_image(i,j)-left_image(i-DIRX,j-DIRY));
          accumulated_costs(i,j) = evaluate_path( accumulated_costs(i-DIRX,j-DIRY), // Previous pixel
                                                  costs(i,j), pixel_diff, false );         // Current pixel
                                                  
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
      for ( int32 i = (DIRX <= 0 ? 0 : 1 ); i < WIDTH; i++ ) {
        accumulated_costs(i,0) += costs(i,0); // TODO: Replace some contional logic by replacing += with = ?
      }
      for ( int32 j = 1; j < HEIGHT; j++ ) {
        int32 istart = std::max( (DIRX <= 0 ? 0 : 1),
                                 (DIRX <= 0 ? 0 : 1) + DIRX * j );
        int32 istop  = std::min( WIDTH, WIDTH + DIRX * j );
        for ( int32 i = istart; i < istop; i++ ) {
          int pixel_diff         = abs(left_image(i,j)-left_image(i-DIRX,j-DIRY));
          accumulated_costs(i,j) = evaluate_path( accumulated_costs(i-DIRX,j-DIRY), // Previous pixel
                                                  costs(i,j), pixel_diff );         // Current pixel
        }
      }
    } 
    if ( DIRX < 0 ) {
      // RIGHT MOST EDGE
      // Process every pixel along this edge only if DIRY == 0. Otherwise skip the top right most pixel
      for ( int32 j = (DIRY <= 0 ? 0 : 1); j < HEIGHT; j++ ) {
        accumulated_costs(WIDTH-1,j) += costs(WIDTH-1,j);
      }
      for ( int32 i = WIDTH-2; i >= 0; i-- ) {
        int32 jstart = std::max( (DIRY <= 0 ? 0 : 1),
                                 (DIRY <= 0 ? 0 : 1) - DIRY * (i - WIDTH + 1) );
        int32 jstop  = std::min( HEIGHT, HEIGHT - DIRY * (i - WIDTH + 1) );
        for ( int32 j = jstart; j < jstop; j++ ) {
          int pixel_diff         = abs(left_image(i,j)-left_image(i-DIRX,j-DIRY));
          accumulated_costs(i,j) = evaluate_path( accumulated_costs(i-DIRX,j-DIRY), // Previous pixel
                                                  costs(i,j), pixel_diff );         // Current pixel
        }
      }
    } 
    if ( DIRY < 0 ) {
      // BOTTOM MOST EDGE
      // Process every pixel along this edge only if DIRX == 0. Otherwise skip the bottom left and bottom right pixel
      for ( int32 i = (DIRX <= 0 ? 0 : 1);
            i < (DIRX >= 0 ? WIDTH : WIDTH-1); i++ ) {
        accumulated_costs(i,HEIGHT-1) += costs(i,HEIGHT-1);
      }
      for ( int32 j = HEIGHT-2; j >= 0; j-- ) {
        int32 istart = std::max( (DIRX <= 0 ? 0 : 1),
                                 (DIRX <= 0 ? 0 : 1) - DIRX * (j - HEIGHT + 1) );
        int32 istop  = std::min( (DIRX >= 0 ? WIDTH : WIDTH - 1),
                                 (DIRX >= 0 ? WIDTH : WIDTH - 1) - DIRX * (j - HEIGHT + 1) );
        for ( int32 i = istart; i < istop; i++ ) {
          int pixel_diff         = abs(left_image(i,j)-left_image(i-DIRX,j-DIRY));
          accumulated_costs(i,j) = evaluate_path( accumulated_costs(i-DIRX,j-DIRY), // Previous pixel
                                                  costs(i,j), pixel_diff );         // Current pixel
        }
      }
    }
  } // End function iterate_direction

  /// im1 += im2.
  // ADD two image views of vector type. Vectors are not from pixel math base
  //
  // This only works if these are ImageViews!
  template <class PixelT>
  void inplace_sum_views( ImageView<PixelT>      & im1,
                          ImageView<PixelT> const& im2 ) {
          PixelT * im1_ptr = im1.data();
    const PixelT * im2_ptr = im2.data();
    while ( im1_ptr != im1.data() + im1.cols()*im1.rows() ) {
      *im1_ptr += *im2_ptr;
      im1_ptr++;
      im2_ptr++;
    }
  }

  // TODO: Typedef uint8 as the disparity unit!

  /// Get the index if the smallest element in a vector
  inline uint8 find_min_index( AVector const& v ) {
    return std::distance(v.begin(),
                         std::min_element(v.begin(), v.end())); // Gets the min iterator
  }

  /// Returns a cost score at a given location
  int16 get_cost(ImageView<uint8> const& left_image,
                 ImageView<uint8> const& right_image,
                 int left_x, int left_y, int right_x, int right_y, bool debug);

  /// Goes across all the viterbi diagrams and extracts out the minimum vector.
  ImageView<Vector<uint8,2> >
  create_disparity_view( ImageView<AVector> const& accumulated_costs );

  /// Invokes a 8 path version of SGM
  ImageView<Vector<uint8,2> >
  semi_global_matching_func( ImageView<uint8> const& left_image,
                             ImageView<uint8> const& right_image );



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
    inline pixel_type operator()( int32 /*i*/, int32 /*j*/, int32 /*p*/ = 0) const {
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
  }

} // end namespace stereo
} // end namespace vw

#endif //__SEMI_GLOBAL_MATCHING__
