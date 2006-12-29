#ifndef __VW_STEREO_PYRAMID_CORRELATOR_H__
#define __VW_STEREO_PYRAMID_CORRELATOR_H__

#include <vw/Image/ImageView.h>
#include <vw/Image/Algorithms.h>
#include <vw/Math/BBox.h>
#include <vw/Stereo/Correlate.h>
#include <vw/Stereo/BlockCorrelator.h>

namespace vw { 
namespace stereo {
  
  class PyramidCorrelator {
    
    int m_kernel_width, m_kernel_height;
    int m_min_h, m_max_h, m_min_v, m_max_v;
    int m_verbose;
    double m_cross_corr_threshold;
    double m_slog_sigma;
    bool m_do_horizontal_calibration;
    bool m_do_vertical_calibration;
    bool m_debug_mode;
    std::string m_debug_prefix;

  public:

    PyramidCorrelator(int minH,	/* left bound disparity search window*/
                      int maxH,	/* right bound disparity search window*/
                      int minV,	/* bottom bound disparity search window */ 
                      int maxV,	/* top bound disparity search window */
                      int kernWidth,	/* size of the kernel */
                      int kernHeight,       
                      int verbose,
                      double cross_corr_threshold,
                      double slog_sigma) : m_min_h(minH), m_max_h(maxH), 
                                           m_min_v(minV), m_max_v(maxV),
                                           m_kernel_width(kernWidth), m_kernel_height(kernHeight),
                                           m_verbose(verbose), m_cross_corr_threshold(cross_corr_threshold),
                                           m_slog_sigma(slog_sigma) {
      m_debug_mode = false;
    }
  
    /// Enable debugging output.  The pyramid correlator will write
    /// images to disk using the provided prefix.
    void enable_debug_mode(std::string const& debug_prefix) {
      m_debug_mode = true;
      m_debug_prefix = debug_prefix;
    }
    void disable_debug_mode() { m_debug_mode = false; }

    // Run the correlator
    template <class PixelT>
    vw::BBox2i operator()(vw::ImageView<PixelT>& image0,
                          vw::ImageView<PixelT>& image1, 
                          bool do_horizontal_search = true,
                          bool do_vertical_search = true) {
      m_do_horizontal_calibration = do_horizontal_search;
      m_do_vertical_calibration = do_vertical_search;

      // Check to make sure that image0 and image1 have equal dimensions 
      if ((image0.cols() != image1.cols()) ||
          (image0.rows() != image1.rows())) {
        vw_throw( ArgumentErr() << "Primary and secondary image dimensions do not agree!" );
      }
        
      // Check to make sure that the images are single channel/single plane
      if (!(image0.channels() == 1 && image0.planes() == 1 &&
            image1.channels() == 1 && image1.planes() == 1)) {
        vw_throw( ArgumentErr() << "Both images must be single channel/single plane images!" );
      }

      ImageView<typename PixelChannelType<PixelT>::type> l_img = channels_to_planes(image0);
      ImageView<typename PixelChannelType<PixelT>::type> r_img = channels_to_planes(image1);        
      
      // First we hone in on the correct search range using a pyramid of images
      this->pyramid_search( l_img, r_img);
      
      return BBox2i(Vector2i(m_min_h, m_min_v), Vector2i(m_max_h, m_max_v));
    }
    
  private:
    // Reduce the image size by a factor of two by averaging the pixels
    template <class PixelT>
    ImageView<PixelT> subsample_and_average(ImageView<PixelT> const& img) {
  
      ImageView<PixelT> outImg(img.cols()/2, img.rows()/2,img.planes());		
      unsigned int i, j, p;
      
      for (p = 0; p < outImg.planes() ; p++) {
        for (i = 0; i < outImg.cols(); i++) {
          for (j = 0; j < outImg.rows(); j++) {  
            outImg(i,j,p) = 0.0f;
            outImg(i,j,p) += img(2*i     , 2*j    ,p);
            outImg(i,j,p) += img(2*i + 1 , 2*j    ,p);
            outImg(i,j,p) += img(2*i     , 2*j + 1,p);
            outImg(i,j,p) += img(2*i + 1 , 2*j + 1,p);
            outImg(i,j,p) /= 4;
          }
        }
      }
      return outImg;
    }

    // This is the core routine for pyramid correlation.  It builds a
    // pyramid of slog'd, subsampled images, and then runs the
    // correlator starting with the smallest image pair.  As each
    // level is completed, the disparity map from that level is used
    // to form a guess about the search range at the next level.  We
    // then pop up a level and refine our estimate of the search
    // range.
    template <class PixelT>
    void pyramid_search(ImageView<PixelT> const& left_image, ImageView<PixelT> const& right_image) {
      const int max_dimension = 256;
      int width = left_image.cols();
      int height = left_image.rows();

      std::vector<ImageView<PixelGray<uint8> > > left_pyramid;
      std::vector<ImageView<PixelGray<uint8> > > right_pyramid;

      ImageView<PixelT> current_left = subsample_and_average(left_image);
      ImageView<PixelT> current_right = subsample_and_average(right_image);

      // Create a pyramid of SLOG images for us to run through the correlator
      std::cout << "Creating image pyramid... " << std::flush;
      int num_levels = 0;
      while (current_left.cols()*2 > max_dimension && current_left.rows()*2 > max_dimension) {
        left_pyramid.push_back(vw::channel_cast<vw::uint8>(vw::threshold(vw::laplacian_filter(vw::gaussian_filter(current_left,m_slog_sigma)), 0.0)) );
        right_pyramid.push_back(vw::channel_cast<vw::uint8>(vw::threshold(vw::laplacian_filter(vw::gaussian_filter(current_right,m_slog_sigma)), 0.0)) );
        current_left = subsample_and_average(current_left);
        current_right = subsample_and_average(current_right);
        num_levels++;
      }
      std::cout << num_levels << " levels.\n";
      if (num_levels == 0) {
        // vw_throw( CorrelatorErr() << "PyramidCorrelator failed.  Image pyramid contained zero levels.!" );
        std::cout <<"ZERO LEVELS!!!\n" << std::flush;
        exit(0);
      }
      
      // Now, run the correlator on each level
      for (int i = num_levels-1; i >= 0; --i) {
        int min_h = int(ceil(m_min_h / pow(2,i+1)));
        int max_h = int(ceil(m_max_h / pow(2,i+1)));
        int min_v = int(ceil(m_min_v / pow(2,i+1)));
        int max_v = int(ceil(m_max_v / pow(2,i+1)));
        
        // Print out some useful feedback for the user
        if(m_verbose) {
          std::cout << "\nPyramid search: level " << i << "  (" << left_pyramid[i].cols()<< " " << left_pyramid[i].rows() << ")\n";
          std::cout << "\tRange -- H [" << min_h << ", " << max_h 
                    << "   V [" << min_v << ", " << max_v << "]\n";
        }
  
        // Create an instance of the optimized correlator.  We run the
        // correlator on one of the subsampled, slog'd image pairs and
        // the resulting disparity map is used to narrow down the
        // search range.
        BlockCorrelator corr(min_h, max_h, min_v, max_v, 
                             m_kernel_width, m_kernel_height, 
                             true, m_cross_corr_threshold, 2048,  // 2048 is block size
                             false, false);
        ImageView<PixelDisparity<float> > disparity_map = corr(left_pyramid[i], right_pyramid[i], true);

        // For debugging 
        std::ostringstream current_level;
        current_level << i;
        double min_h_disp, min_v_disp, max_h_disp, max_v_disp;
        if (m_debug_mode) {
          disparity::get_disparity_range(disparity_map, min_h_disp, max_h_disp, min_v_disp, max_v_disp,true);
          write_image( m_debug_prefix+"-H-raw-"+current_level.str()+".jpg", normalize(clamp(select_channel(disparity_map,0), min_h_disp, max_h_disp)));
          write_image( m_debug_prefix+"-V-raw-"+current_level.str()+".jpg", normalize(clamp(select_channel(disparity_map,1), min_v_disp, max_v_disp)));
        }

        // Now, clean up the disparity map by rejecting outliers 
        unsigned rm_half_kernel = 5 / (int)pow(2,i+1);
        unsigned rm_min_matches_percent = 60 / (i+1);
        double rm_threshold = 3;
  
        for(int nn=0; nn < 1; nn++) {  // Run the clean up routine three times
          disparity::clean_up(disparity_map,
                              rm_half_kernel, 
                              rm_half_kernel,
                              rm_min_matches_percent,
                              rm_threshold,
                              true);
        }
        std::cout << "\tRemoving solitary pixels [20x20 window, 20% threshold]\n";
        disparity::remove_outliers(disparity_map, 20, 20, 20, 200, true);

        int mask_buffer = std::max(m_kernel_width, m_kernel_height);
        disparity::mask(disparity_map,
                        disparity::generate_mask(left_pyramid[i], mask_buffer),
                        disparity::generate_mask(right_pyramid[i], mask_buffer));

        // For debugging
        if (m_debug_mode) {
          disparity::get_disparity_range(disparity_map, min_h_disp, max_h_disp, min_v_disp, max_v_disp,true);
          write_image( m_debug_prefix+"-H-"+current_level.str()+".jpg", normalize(clamp(select_channel(disparity_map,0), min_h_disp, max_h_disp)));
          write_image( m_debug_prefix+"-V-"+current_level.str()+".jpg", normalize(clamp(select_channel(disparity_map,1), min_v_disp, max_v_disp)));
        } 

        // The disparity map of the subsampled problem is used to
        // determine the best guess for search range at the next
        // level.
        try {
          double new_h_min, new_v_min, new_h_max, new_v_max;
          disparity::get_disparity_range(disparity_map, new_h_min, new_h_max, new_v_min, new_v_max,true);
          
          if (m_do_horizontal_calibration) {
            m_min_h = int(floor( new_h_min * pow(2,i+1))) - 4;
            m_max_h = int(ceil(  new_h_max * pow(2,i+1))) + 4;
          }
          
          if (m_do_vertical_calibration) {
            m_min_v = int(floor( new_v_min * pow(2,i+1))) - 4;
            m_max_v = int(ceil(  new_v_max * pow(2,i+1))) + 4;
          }
          
          std::cout << "\n\tNew disparity range  --  H: [" << m_min_h << ", " << m_max_h
                    << "]  V: [" << m_min_v << ", " << m_max_v << "]\n";
          
        } catch (ArgumentErr &e) { // Couldn't adjust disparity range
          std::cout << "INSUFFICIENT MATHCHES\n" << std::flush;
          exit(0);
          // vw_throw( CorrelatorErr() << "Pyramid Align Failed.  Correlation return insufficient results to compute a match." );
        }
        
      }

    }
  };

}} // namespace vw::stereo 

#endif // __PYRAMID_CORRELATOR_H__ 

