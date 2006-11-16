#ifndef __VW_STEREO_DISPARITY_MAP_H__
#define __VW_STEREO_DISPARITY_MAP_H__

#include <vw/Image/PixelTypes.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/Filter.h>
#include <vw/Image/Statistics.h>

#include <vw/FileIO.h>

// For the PixelDisparity math.
#include <boost/operators.hpp>

namespace vw { 
  
  /// The disparity pixel type for vision workbench has two channels of
  /// disparity (horizontal, vertical) and one alpha channel for recording
  /// which pixels are valid.
  ///
  /// An alpha value of 0 denotes a bad (missing) pixel, and an alpha
  /// channel of 1 denotes a good pixel.
  template <class ChannelT>
  class PixelDisparity : 
    boost::additive< PixelDisparity<ChannelT> >,
    boost::multiplicative< PixelDisparity<ChannelT>, ChannelT >
  {
  private:
    ChannelT m_ch[3];
  public:
    PixelDisparity() { m_ch[0]=m_ch[1]=0; m_ch[2]=1; }
    PixelDisparity( ChannelT const& scalar ) { m_ch[0]=m_ch[1]=scalar; m_ch[2] = 0; }
    PixelDisparity( ChannelT const& h, ChannelT const& v ) { m_ch[0]=h; m_ch[1]=v; m_ch[2]=0; }
  
    template <class OtherT> explicit PixelDisparity( PixelDisparity<OtherT> const& other ) { m_ch[0]=other[0]; m_ch[1]=other[1]; m_ch[2]=other[2];}

    inline ChannelT& operator[](int i) { return m_ch[i]; }
    inline ChannelT const& operator[](int i) const { return m_ch[i]; }
    inline ChannelT& h() { return m_ch[0]; }
    inline ChannelT const& h() const { return m_ch[0]; }
    inline ChannelT& v() { return m_ch[1]; }
    inline ChannelT const& v() const { return m_ch[1]; }
    inline ChannelT& missing() { return m_ch[2]; }
    inline ChannelT missing() const { return m_ch[2]; }
    inline double magnitude() const { return sqrtf((double)m_ch[0]*m_ch[0] + (double)m_ch[1]*m_ch[1]); }
    inline double magnitude_squared() const { return (double)m_ch[0]*m_ch[0] + (double)m_ch[1]*m_ch[1]; }

    inline PixelDisparity& operator-() { if (!missing()) {m_ch[0]=-m_ch[0]; m_ch[1]=-m_ch[1];} return *this; }
    inline PixelDisparity& operator+=( PixelDisparity const& p ) { if (missing() || p.missing()) {m_ch[0]=m_ch[1]=0; m_ch[2] = 1;} 
      else { m_ch[0]+=p[0]; m_ch[1]+=p[1]; } return *this; }
    inline PixelDisparity& operator-=( PixelDisparity const& p ) { if (missing() || p.missing()) {m_ch[0]=m_ch[1]=0; m_ch[2] = 1;} 
      else { m_ch[0]-=p[0]; m_ch[1]-=p[1]; } return *this; }
    inline PixelDisparity& operator*=( ChannelT s ) { if (!missing()) { m_ch[0]*=s; m_ch[1]*=s; } return *this; }
    inline PixelDisparity& operator/=( ChannelT s ) { if (!missing()) { m_ch[0]/=s; m_ch[1]/=s; } return *this; }
  };

  /// \cond INTERNAL
  VW_DECLARE_PIXEL_TYPE(PixelDisparity,3);
  /// \endcond

  template <class ChannelT>
  std::ostream& operator<<( std::ostream& os, PixelDisparity<ChannelT> const& pix ) {
    if (pix.missing()) 
      return os << "Disparity( MISSING )";
    else 
      return os << "Disparity(" << pix.h() << "," << pix.v() << ")";
  }
} // namespace vw

namespace vw {
namespace stereo {
namespace disparity {
  
  /// Masks all of the black pixels along the edges of an image.  This
  /// algorithm "eats away" at the pixels on all four sides of the
  /// image; masking pixels until it encounters a non-black pixel.
  /// 
  /// You can supply an optional buffer argument that will mask some
  /// of the good data bordering the image as well.  The width of this
  /// border is set using the additional_border_width argument (in
  /// units of pixels).
  template <class ViewT>
  vw::ImageView<bool> generate_mask(vw::ImageViewBase<ViewT> const& input_image,
                                    unsigned int additional_border_width = 0) {

    const double BLACK_PIXEL = 0;
    vw::ImageView<bool> mask(input_image.impl().cols(), input_image.impl().rows(), 1);
    unsigned int i, j;

    vw::fill(mask, true);

    for (i = 0; i < input_image.impl().cols(); i++) {
      // Search from the left side of the image for black pixels 
      j = 0;
      while ( j < input_image.impl().rows() && input_image.impl()(i,j)[0] == BLACK_PIXEL ) {
        mask(i,j) = false;
        j++;
      }

      // Mask additional "buffer" pixels
      int j_start = j;
      while ( j < input_image.impl().rows() && j - j_start < additional_border_width ) {
        mask(i,j) = false;
        j++;
      }
    
      // Search from the right side of the image for black pixels 
      j = input_image.impl().rows() - 1;
      while ( j > 0 && input_image.impl()(i,j)[0] == BLACK_PIXEL ) {
        mask(i,j) = false;
        j--;
      }

      // Mask additional "buffer" pixels
      int j_end = j;
      while ( j > 0 && j_end - j < additional_border_width) {
        mask(i,j) = false;
        j--;
      }
    }

    for (j = 0; j < input_image.impl().rows(); j++) {
      // Search from the top side of the image for black pixels 
      i = 0;
      while ( i < input_image.impl().cols() && input_image.impl()(i,j)[0] == BLACK_PIXEL ) {
        mask(i,j) = false;
        i++;
      }

      // Mask additional "buffer" pixels
      int i_start = i;
      while ( i < input_image.impl().cols() && i - i_start < additional_border_width) {
        mask(i,j) = false;
        i++;
      }
    
      // Search from the bottom side of the image for black pixels 
      i = input_image.impl().cols() - 1;
      while ( i > 0 && input_image.impl()(i,j)[0] == BLACK_PIXEL) {
        mask(i,j) = false;
        i--;
      }

      // Mask additional "buffer" pixels
      int i_end = i;
      while ( i > 0 && i_end - i < additional_border_width) {
        mask(i,j) = false;
        i--;
      }
    }

    return mask;
  }

  /// Apply a binary mask to the disparity map (see also \ref disparity::generate_mask())
  template <class PixelT>
  void mask(ImageView<PixelDisparity<PixelT> > &disparity_map, 
            ImageView<bool> const& left_mask,
            ImageView<bool> const& right_mask) {
    
    VW_ASSERT(disparity_map.cols() == left_mask.cols() && disparity_map.rows() == left_mask.rows() &&
              disparity_map.cols() == right_mask.cols() && disparity_map.rows() == right_mask.rows(),
              ArgumentErr() << "disparity::mask() : Mask images and disparity map image are not the same size.\n");
    
    for (unsigned i = 0; i < disparity_map.cols() ; i++) 
      for (unsigned j = 0; j < disparity_map.rows() ; j++)
        if ( disparity_map(i,j).missing() || !left_mask(i,j) || !right_mask((int)(i+disparity_map(i,j).h()), (int)(j+disparity_map(i,j).v())) ) 
          disparity_map(i,j) = PixelDisparity<PixelT>();  // Set to missing pixel value
  }
  
  /// Remove outliers from a disparity map image
  template <class PixelT>
  void remove_outliers(ImageView<PixelDisparity<PixelT> > &disparity_map,
                       int half_h_kernel, int half_v_kernel,
                       int min_matches, double threshold, bool verbose = false) {
    
    unsigned int width = disparity_map.cols();
    unsigned int height = disparity_map.rows();
    int	matched;
    int	total;
    unsigned int	x, y, xk, yk;
    ImageView<bool> disparity_mask(width, height);
    fill(disparity_mask, false);  
    double rejection_threshold = (double)min_matches/100.0;
    
    if(verbose) {
      printf("Removing low confidence pixels");
      fflush(stdout);
    }
      
    for(y = 0; y < height; y++){
      for(x = 0; x < width; x++){
        // if valid pixel 
        if ( !disparity_map(x,y).missing() ) {
          // walk the kernel 
          matched = 0;
          total = 0;
          for(yk = y - half_v_kernel; yk <= y + half_v_kernel; yk++){
            if(yk >= 0 && yk < height){
              for(xk = x - half_h_kernel; xk <= x + half_h_kernel; xk++){
                if(xk >=0 && xk < width){
                  if(fabs(disparity_map(x,y).h()-disparity_map(xk, yk).h()) <= threshold &&
                     fabs(disparity_map(x,y).v()-disparity_map(xk, yk).v()) <= threshold) {
                    matched++;
                  }
                  total++;
                }
              }
            }
          } // end walk kernel 

          if(total != 0){
            if( ((double)matched/(double)total) >= rejection_threshold){
              disparity_mask(x,y) = true;
            }
          }
        } // end if valid pixel 
      }
    }
    total = 0;
    for(y = 0; y < height; y++){
      for(x = 0; x < width; x++){
        if( !disparity_mask(x,y) ) {
          disparity_map(x,y) = PixelDisparity<PixelT>(); // Reset to missing pixel 
          total++;
        }
      }
    }
    if(verbose) {
      printf("\r        %d/%d low confidence pixels removed (%0.2f%%)\n", total, width*height, (double)total/(width*height));
      fflush(stdout);
    }
  }


  /// Clean up a disparity map.
  ///
  /// You supply the half dimensions of the kernel window.  
  ///
  /// Next, you supply the percentage of the pixels within the kernel
  /// that must "match" the center pixel if that pixel is to be
  /// considered an inlier. (given in units of percent [0..100]).
  ///
  /// Finally, you supply the threshold that determines whether a
  /// pixel is considered "close" to its neightbors (in units of
  /// pixels).
  template <class PixelT>
  inline void clean_up(ImageView<PixelDisparity<PixelT> > &disparity_map,
                       int v_half_kernel, int h_half_kernel,
                       int min_matches, double threshold,
                       bool verbose = false) {
    
    // Remove outliers using user specified parameters
    remove_outliers(disparity_map, 
                    v_half_kernel, h_half_kernel, 
                    min_matches, threshold, verbose);
    
    // Remove outliers using a heuristic that isolates single pixel
    // outliers.
    remove_outliers(disparity_map, 1, 1, 75, 0.5, verbose);
  }


  template <class PixelT>
  inline void disparity_debug_images(ImageView<PixelDisparity<PixelT> > const& disparity_map,
                                     ImageView<PixelGray<float> > &horizontal,
                                     ImageView<PixelGray<float> > &vertical) {
    double min_h_disp, min_v_disp, max_h_disp, max_v_disp;
    get_disparity_range(disparity_map, min_h_disp, max_h_disp, min_v_disp, max_v_disp);
    horizontal = clamp(select_channel(disparity_map,0), min_h_disp, max_h_disp);
    vertical = clamp(select_channel(disparity_map,1), min_v_disp, max_v_disp);
  }


  template <class PixelT>
  inline ImageView<PixelRGB<float> > rgb_missing_pixel_image(ImageView<PixelDisparity<PixelT> > const& disparity_map) {
    ImageView<PixelRGB<float> > mask(disparity_map.cols(), disparity_map.rows());

    for (unsigned i = 0; i < mask.cols(); i++) {
      for (unsigned j = 0; j < mask.rows(); j++) {
        if ( !disparity_map(i,j).missing() ) {
          mask(i,j).r() = 0.8;
          mask(i,j).g() = 0.8;
          mask(i,j).b() = 0.8;
        } else {
          mask(i,j).r() = 1.0;
          mask(i,j).g() = 0.0;
          mask(i,j).b() = 0.0;
        }
      }
    }
    return mask;
  }

  template <class PixelT>
  inline ImageView<PixelGray<float> > missing_pixel_image(ImageView<PixelDisparity<PixelT> > const& disparity_map) {
    ImageView<PixelGray<float> > mask(disparity_map.cols(), disparity_map.rows());

    for (unsigned i = 0; i < mask.cols(); i++) {
      for (unsigned j = 0; j < mask.rows(); j++) {
        if ( !disparity_map(i,j).missing() ) {
          mask(i,j).v() = 1.0;
        } else {
          mask(i,j).v() = 0.0;
        }
      }
    }
    return mask;
  }


  template <class PixelT>
  void get_disparity_range(ImageView<PixelDisparity<PixelT> > const& disparity_map, 
                           double &min_horz_disp, double& max_horz_disp, 
                           double &min_vert_disp, double& max_vert_disp,
                           bool verbose = false) {
    
    max_horz_disp = -1e100;
    min_horz_disp = 1e100;
    max_vert_disp = -1e100;
    min_vert_disp = 1e100;

    // Find the max/min disparity values
    int missing = 0;
    for (unsigned i = 0; i < disparity_map.cols(); i++) {
      for (unsigned j = 0; j < disparity_map.rows(); j++) {
        if ( !disparity_map(i,j).missing() ) {
          max_horz_disp = disparity_map(i,j).h() > max_horz_disp ? disparity_map(i,j).h() : max_horz_disp;
          min_horz_disp = disparity_map(i,j).h() < min_horz_disp ? disparity_map(i,j).h() : min_horz_disp;
          max_vert_disp = disparity_map(i,j).v() > max_vert_disp ? disparity_map(i,j).v() : max_vert_disp;
          min_vert_disp = disparity_map(i,j).v() < min_vert_disp ? disparity_map(i,j).v() : min_vert_disp;
        } else {
          missing++;
        }

      }
    }
    
    if (verbose) {
      printf("Disparity range -- Horizontal: [%f, %f]   Vertical: [%f, %f]  (%d missing)\n", 
             min_horz_disp, max_horz_disp, min_vert_disp, max_vert_disp, missing);
    }
  }

  template <class ChannelT>
  void sparse_disparity_filter(ImageView<PixelDisparity<ChannelT> > &disparity_map, float blur_stddev, float rejection_threshold) {
    std::cout << "\tIsolating and rejecting large areas of very low confidence..." << std::flush;
    
    ImageView<PixelGray<float> > test_image = disparity::missing_pixel_image(disparity_map);
    //    write_image("test-a.png", test_image);
    ImageView<float> blurred_image = gaussian_filter(select_channel(test_image,0), blur_stddev);
    //    write_image("test-b.png", blurred_image);
    ImageView<float> threshold_image = threshold(blurred_image, rejection_threshold);
    //    write_image("test-c.png", threshold_image);
    
    for (int i = 0; i < disparity_map.cols(); i++) 
      for (int j = 0; j < disparity_map.rows(); j++) 
        if (threshold_image(i,j) == 0)
          disparity_map(i,j) = PixelDisparity<ChannelT>(); // Set to missing pixel
    
    std::cout << " done.\n";
  }

  template <class ChannelT, class ImagePixelT>
  void low_contrast_filter(ImageView<PixelDisparity<ChannelT> > &disparity_map, 
                           ImageView<ImagePixelT> const& left_image,
                           ImageView<ImagePixelT> const& right_image,
                           int kernel_width, int kernel_height,
                           float rejection_threshold) {
    std::cout << "\tIsolating and rejecting large areas of low contrast..." << std::flush;
    
    ImageView<float> left_contrast_image(left_image.cols(), left_image.rows());
    ImageView<float> right_contrast_image(right_image.cols(), right_image.rows());

    // First, compute the standard deviation in each image patch
    for (int j = kernel_height/2; j < left_image.rows()-kernel_height/2; ++j) {
      for (int i = kernel_width/2; i < left_image.cols()-kernel_width/2; ++i) {
        typename CompoundChannelType<ImagePixelT>::type std_dev;
        int crop_i = i - kernel_width/2;
        int crop_j = j - kernel_height/2;
        left_contrast_image(i,j) = stddev_channel_value(crop(channel_cast<float>(left_image), crop_i, crop_j, kernel_width, kernel_height));
        right_contrast_image(i,j) = stddev_channel_value(crop(channel_cast<float>(right_image), crop_i, crop_j, kernel_width, kernel_height));
      }
    }

    // For debugging
    //     write_image("left-contrast.png", normalize(left_contrast_image));
    //     write_image("right-contrast.png", normalize(right_contrast_image));
    
    // Reject pixels that have a standard deviation below the supplied threshold
    for (int j = 0; j < disparity_map.rows(); ++j) {
      for (int i = 0; i < disparity_map.cols(); ++i) {
        if (left_contrast_image(i,j) < rejection_threshold ||
            right_contrast_image(int(i+disparity_map(i,j).h()), int(j+disparity_map(i,j).v())) < rejection_threshold) {
          disparity_map(i,j) = PixelDisparity<ChannelT>();
        }
      }
    }
    std::cout << " done.\n";
  }
  
  
} // namespace disparity
  
}}    // namespace vw::stereo


#endif //  _VWDISPARITYMAP_H_

