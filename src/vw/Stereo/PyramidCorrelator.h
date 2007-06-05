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
    template <class ViewT>
    vw::BBox2i operator()(vw::ImageViewBase<ViewT>& image0, vw::ImageViewBase<ViewT>& image1, 
                          bool do_horizontal_search = true, bool do_vertical_search = true) {
      m_do_horizontal_calibration = do_horizontal_search;
      m_do_vertical_calibration = do_vertical_search;

      // Check to make sure that image0 and image1 have equal dimensions 
      if ((image0.impl().cols() != image1.impl().cols()) || (image0.impl().rows() != image1.impl().rows())) {
        vw_throw( ArgumentErr() << "Primary and secondary image dimensions do not agree!" );
      }
        
      // Check to make sure that the images are single channel/single plane
      if (!(image0.impl().channels() == 1 && image0.impl().planes() == 1 &&
            image1.impl().channels() == 1 && image1.impl().planes() == 1)) {
        vw_throw( ArgumentErr() << "Both images must be single channel/single plane images!" );
      }

      ImageView<typename PixelChannelType<typename ViewT::pixel_type>::type> l_img = channels_to_planes(image0.impl());
      ImageView<typename PixelChannelType<typename ViewT::pixel_type>::type> r_img = channels_to_planes(image1.impl());        
      
      // First we hone in on the correct search range using a pyramid of images
      this->pyramid_search( l_img, r_img);
      
      return BBox2i(Vector2i(m_min_h, m_min_v), Vector2i(m_max_h, m_max_v));
    }
    
  private:
    // Reduce the image size by a factor of two by averaging the pixels
    template <class ViewT>
    ImageView<typename ViewT::pixel_type> subsample_and_average(ImageViewBase<ViewT> const& img) {
  
      ImageView<typename ViewT::pixel_type> outImg(img.impl().cols()/2, img.impl().rows()/2, img.impl().planes());		
      int32 i, j, p;
      
      for (p = 0; p < outImg.planes() ; p++) {
        for (i = 0; i < outImg.cols(); i++) {
          for (j = 0; j < outImg.rows(); j++) {  
            outImg(i,j,p) = 0.0f;
            outImg(i,j,p) += img.impl()(2*i     , 2*j    ,p);
            outImg(i,j,p) += img.impl()(2*i + 1 , 2*j    ,p);
            outImg(i,j,p) += img.impl()(2*i     , 2*j + 1,p);
            outImg(i,j,p) += img.impl()(2*i + 1 , 2*j + 1,p);
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
    template <class ViewT>
    void pyramid_search(ImageViewBase<ViewT> const& left_image, ImageViewBase<ViewT> const& right_image) {
      const int max_dimension = 256;
      int width = left_image.impl().cols();
      int height = left_image.impl().rows();

      std::vector<ImageView<PixelGray<uint8> > > left_pyramid;
      std::vector<ImageView<PixelGray<uint8> > > right_pyramid;

      ImageView<typename ViewT::pixel_type> current_left = subsample_and_average(left_image);
      ImageView<typename ViewT::pixel_type> current_right = subsample_and_average(right_image);

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
          BBox2i disp_range = disparity::get_disparity_range(disparity_map);
          write_image( m_debug_prefix+"-H-raw-"+current_level.str()+".jpg", normalize(clamp(select_channel(disparity_map,0), disp_range.min().x(), disp_range.max().x())));
          write_image( m_debug_prefix+"-V-raw-"+current_level.str()+".jpg", normalize(clamp(select_channel(disparity_map,1), disp_range.min().y(), disp_range.max().y())));
        }

        // Now, clean up the disparity map by rejecting outliers 
        int32 rm_half_kernel = 5;
        int32 rm_min_matches_percent = 60;
        double rm_threshold = 3;
  
        ImageViewRef<PixelDisparity<float> > processed_disparity_map_ref = disparity::clean_up(disparity_map,
                                                                                               rm_half_kernel, 
                                                                                               rm_half_kernel,
                                                                                               rm_threshold,
                                                                                               rm_min_matches_percent/100.0);
        //        processed_disparity_map_ref = disparity::remove_outliers(processed_disparity_map_ref, 20, 20, 2, 0.2);

        int mask_buffer = std::max(m_kernel_width, m_kernel_height);
        processed_disparity_map_ref = disparity::mask(processed_disparity_map_ref,
                                                      disparity::generate_mask(left_pyramid[i], mask_buffer),
                                                      disparity::generate_mask(right_pyramid[i], mask_buffer));
        
        ImageView<PixelDisparity<float> > processed_disparity_map = processed_disparity_map_ref;
        
        // For debugging
        if (m_debug_mode) {
          BBox2i disp_range = disparity::get_disparity_range(processed_disparity_map);
          write_image( m_debug_prefix+"-H-"+current_level.str()+".jpg", normalize(clamp(select_channel(processed_disparity_map,0), disp_range.min().x(), disp_range.max().x())));
          write_image( m_debug_prefix+"-V-"+current_level.str()+".jpg", normalize(clamp(select_channel(processed_disparity_map,1), disp_range.min().y(), disp_range.max().y())));
        } 

        // The disparity map of the subsampled problem is used to
        // determine the best guess for search range at the next
        // level.
        try {

          BBox2i disp_range = disparity::get_disparity_range(processed_disparity_map);
          
          if (m_do_horizontal_calibration) {
            m_min_h = int(floor( disp_range.min().x() * pow(2,i+1))) - 4;
            m_max_h = int(ceil(  disp_range.max().x() * pow(2,i+1))) + 4;
          }
          
          if (m_do_vertical_calibration) {
            m_min_v = int(floor( disp_range.min().y() * pow(2,i+1))) - 4;
            m_max_v = int(ceil(  disp_range.max().y() * pow(2,i+1))) + 4;
          }
          
          std::cout << "\n\tNew disparity range  --  H: [" << disp_range.min().x() << ", " << disp_range.max().x()
                    << "]  V: [" << disp_range.min().y() << ", " << disp_range.max().y() << "]\n";
          
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

