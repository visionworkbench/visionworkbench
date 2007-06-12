#ifndef __VW_STEREO_MULTIRESOLUTION_CORRELATOR__
#define __VW_STEREO_MULTIRESOLUTION_CORRELATOR__
#include <vw/Stereo/Correlate.h>
#include <vw/Stereo/OptimizedCorrelator.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Image/Transform.h>
#include <vw/Core/Debugging.h>

namespace vw {
namespace stereo {

  // Reduce the image size by a factor of two by averaging the pixels
  template <class PixelT>
  ImageView<PixelT> subsample_by_two(ImageView<PixelT> &img) {
    
    ImageView<PixelT> outImg(img.cols()/2, img.rows()/2,img.planes());		
    int32 i, j, p;
    
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
  
  class MultiresolutionCorrelator {
    
    int m_lKernWidth, m_lKernHeight;
    int m_lMinH, m_lMaxH, m_lMinV, m_lMaxV;
    int m_verbose;
    double m_crossCorrThreshold;
    float m_corrscore_rejection_threshold;
    int m_useHorizSubpixel;
    int m_useVertSubpixel;
    double m_slog_width;
    int m_min_dimension;
    std::string m_debug_prefix;

  public:
    MultiresolutionCorrelator(int minH,	/* left bound disparity search window*/
                              int maxH,	/* right bound disparity search window*/
                              int minV,	/* bottom bound disparity search window */ 
                              int maxV,	/* top bound disparity search window */
                              int kernWidth,	/* size of the kernel */
                              int kernHeight,       
                              int verbose,
                              double crosscorrThreshold,
                              float corrscore_rejection_threshold,
                              double slog_width,
                              int useSubpixelH,
                              int useSubpixelV) {
      m_lKernWidth = kernWidth;
      m_lKernHeight = kernHeight;
      m_lMinH = minH;
      m_lMaxH = maxH;
      m_lMinV = minV;
      m_lMaxV = maxV;  
      m_verbose = verbose;
      m_slog_width = slog_width;
    
      m_crossCorrThreshold = crosscorrThreshold;
      m_corrscore_rejection_threshold = corrscore_rejection_threshold;
      m_useHorizSubpixel = useSubpixelH;
      m_useVertSubpixel = useSubpixelV;
      m_min_dimension = 512;

      m_debug_prefix = "";
    }

    void set_min_dimension(int value) { m_min_dimension = value; }
    void set_debug_mode(std::string const& debug_file_prefix) { m_debug_prefix = debug_file_prefix; }

    template <class PixelT>
    ImageView<PixelDisparity<float> > correlate(vw::ImageView<PixelT>& left_image,
                                                vw::ImageView<PixelT>& right_image,
                                                bool swap) {  

      set_debug_level(InfoMessage+1);
    
      // Compute the number of image pyramid levels
      int dimension = std::min(left_image.cols(), left_image.rows());
      int levels = 1;
      while (dimension > m_min_dimension) { levels++; dimension /= 2; }
      std::cout << "\tBuilding image pyramid with " << levels << " levels.\n\n";

      // Build the image pyramid
      std::vector<ImageView<float> > left_pyramid(levels), right_pyramid(levels);
      std::vector<ImageView<uint8> > left_slog_pyramid(levels), right_slog_pyramid(levels);
      std::vector<ImageView<PixelDisparity<float> > > disparity_pyramid(levels);
      left_pyramid[0] = channels_to_planes(left_image);
      right_pyramid[0] = channels_to_planes(right_image);
      left_slog_pyramid[0] = vw::channel_cast<uint8>(vw::threshold(vw::laplacian_filter(vw::gaussian_filter(left_image,m_slog_width)), 0.0));
      right_slog_pyramid[0] = vw::channel_cast<uint8>(vw::threshold(vw::laplacian_filter(vw::gaussian_filter(right_image,m_slog_width)), 0.0));
    
      // Apply the SLOG filter to each level.  
      for (int n = 1; n < levels; ++n) {
        left_pyramid[n] = subsample_by_two(left_pyramid[n-1]);
        right_pyramid[n] = subsample_by_two(right_pyramid[n-1]);
        left_slog_pyramid[n] = vw::channel_cast<uint8>(vw::threshold(vw::laplacian_filter(vw::gaussian_filter(left_pyramid[n],m_slog_width)), 0.0));
        right_slog_pyramid[n] = vw::channel_cast<uint8>(vw::threshold(vw::laplacian_filter(vw::gaussian_filter(right_pyramid[n],m_slog_width)), 0.0));
      }

      int h_min, h_max, v_min, v_max;
      if (swap) {
        h_min = int(-m_lMaxH / pow(2,levels-1));
        h_max = int(-m_lMinH / pow(2,levels-1));
        v_min = int(-m_lMaxV / pow(2,levels-1));
        v_max = int(-m_lMinV / pow(2,levels-1));
      } else {
        h_min = int(m_lMinH / pow(2,levels-1));
        h_max = int(m_lMaxH / pow(2,levels-1));
        v_min = int(m_lMinV / pow(2,levels-1));
        v_max = int(m_lMaxV / pow(2,levels-1));
      }

      int h_kern = m_lKernWidth;
      int v_kern = m_lKernHeight;

      // Get things started with a full correlation.
      OptimizedCorrelator correlator(h_min, h_max,
                                     v_min, v_max,
                                     h_kern, v_kern,
                                     true,          // verbose
                                     m_crossCorrThreshold,
                                     m_corrscore_rejection_threshold,
                                     false, false); // no subpixel for now
      disparity_pyramid[levels-1] = correlator(left_slog_pyramid[levels-1], 
                                               right_slog_pyramid[levels-1], 
                                               true);

      // Print Out the range of disparity values
      BBox2 disp_range = disparity::get_disparity_range(disparity_pyramid[levels-1], true);

      // Print out the disparity map at the lowest resolution
      if (m_debug_prefix.size() > 0) {
        write_image( m_debug_prefix + "-DH-0.jpg", normalize(clamp(select_channel(disparity_pyramid[levels-1],0), disp_range.min().x(), disp_range.max().x() )));
        write_image( m_debug_prefix + "-DV-0.jpg", normalize(clamp(select_channel(disparity_pyramid[levels-1],1), disp_range.min().y(), disp_range.max().y() )));
      }

      //     // Now, clean up the disparity map by rejecting outliers 
      //     int32 rm_half_kernel = 5 / (int32)powf(2,levels) + 1;
      //     int32 rm_min_matches_percent = 60 / levels;
      //     double rm_threshold = 3;
      //     for(int nn=0; nn < 1; nn++) {  // Run the clean up routine three times
      //       disparity::clean_up(disparity_pyramid[levels-1],
      //                           rm_half_kernel, 
      //                           rm_half_kernel,
      //                           rm_min_matches_percent,
      //                           rm_threshold,
      //                           true);
      //     }
      //     std::cout << "\tRemoving solitary pixels [20x20 window, 20% threshold]\n";
      //     disparity::remove_outliers(disparity_pyramid[levels-1], 20, 20, 20, 200, true);
    
      // Clear out any regions with a large number of outliers
      // 
      // These settings should move from being hard coded to being user configurable
      //
      // FIXME: This should be re-enabled to use the newer disparity map
      //      filtering architecture that uses PerPixelAccessorViews.
      //      disparity::sparse_disparity_filter(disparity_pyramid[levels-1],
      //      100, 0.1);

      // Print out the disparity map at the lowest resolution
      if (m_debug_prefix.size() > 0) {
        BBox2 disp_range = disparity::get_disparity_range(disparity_pyramid[levels-1]);
        write_image( m_debug_prefix+"-DH-0-filt.jpg", normalize(clamp(select_channel(disparity_pyramid[levels-1],0), disp_range.min().x(), disp_range.max().x() )));
        write_image( m_debug_prefix+"-DV-0-filt.jpg", normalize(clamp(select_channel(disparity_pyramid[levels-1],1), disp_range.min().y(), disp_range.max().y() )));
      }

      // Refined the disparity map by searching in the local region where the last good disparity value was found
      std::cout << "\tRefining disparity map.\n";
      for (int n = levels - 2; n >=0; --n) {

        // Resample the image and multiply the disparity range by 2 to
        // reach the next level of detail.
        disparity_pyramid[n] = 2 * resample(disparity_pyramid[n+1], 2,
                                            (int)left_pyramid[n].cols(), (int)left_pyramid[n].rows(), 
                                            ZeroEdgeExtension(), NearestPixelInterpolation());
      
        // Save a copy of the pre-refinement horiz. disparity for later comparison
        
        ImageView<double> orig;
        if (m_debug_prefix.size() > 0) 
          orig = select_channel(disparity_pyramid[n],0);

        int min_i_search_range, max_i_search_range;
        int min_j_search_range, max_j_search_range;
        if (m_lMaxV-m_lMinV == 0) {
          min_j_search_range = 0; 
          max_j_search_range = 0;
        } else {
          min_j_search_range = -2;
          max_j_search_range = 2;
        }

        if (m_lMaxH-m_lMinH == 0) {
          min_i_search_range = 0;
          max_i_search_range = 0;
        } else {
          min_i_search_range = -2;
          max_i_search_range = 2;
        }
      
        // Refinement: adjust disparity ranges at successively higher resolutions...
        for (int j=v_kern; j<disparity_pyramid[n].rows()-v_kern; ++j) {
          if (j%33 == 0)
            std::cout << "\tLevel " << n << ": " << float(j)/float(disparity_pyramid[n].rows()-2*v_kern)*100.0 << "%               \r" << std::flush;
          for (int i=h_kern; i<disparity_pyramid[n].cols()-h_kern; ++i) {
          
            if (!disparity_pyramid[n](i,j).missing()) {
              double best_score = 1e99;
              double best_i_disp = 0;
              double best_j_disp = 0;

              for (int j_disp = min_j_search_range; j_disp <= max_j_search_range; ++j_disp) {
                for (int i_disp = min_i_search_range; i_disp <= max_i_search_range; ++i_disp) {
                  double current_score = compute_soad(&(left_slog_pyramid[n](0,0)), &(right_slog_pyramid[n](0,0)),
                                                      j,i,
                                                      int(i_disp+disparity_pyramid[n](i,j).h()), 
                                                      int(j_disp+disparity_pyramid[n](i,j).v()),
                                                      h_kern, v_kern,
                                                      left_pyramid[n].cols(), 
                                                      left_pyramid[n].rows());
                  if (current_score < best_score) {
                    best_score = current_score;
                    best_i_disp = i_disp;
                    best_j_disp = j_disp;
                  }
                }              
              }
              disparity_pyramid[n](i,j).h() += best_i_disp;
              disparity_pyramid[n](i,j).v() += best_j_disp;
            }
          }
        }
        std::cout << "\tLevel " << n << ": Done.         \n";

        // Debugging output at each level
        if (m_debug_prefix.size() > 0) {
          std::ostringstream current_level;
          current_level << n;
          BBox2 disp_range = disparity::get_disparity_range(disparity_pyramid[n]);
          write_image( m_debug_prefix+"-refined-H-" + current_level.str() + "-filt.jpg", normalize(clamp(select_channel(disparity_pyramid[n],0), disp_range.min().x(), disp_range.max().x() )));
          write_image( m_debug_prefix+"-refined-V-" + current_level.str() + "-filt.jpg", normalize(clamp(select_channel(disparity_pyramid[n],1), disp_range.min().y(), disp_range.max().y() )));

          std::ostringstream ostream3;
          ostream3 << m_debug_prefix << "-error-" << n << ".jpg";
          write_image(ostream3.str(), normalize(clamp(select_channel(disparity_pyramid[n],0)-orig, -100, 100)));
        }

      }
      return disparity_pyramid[0];
    }

    template <class ViewT>
    ImageView<PixelDisparity<float> > operator()(vw::ImageViewBase<ViewT>& left_image,
                                                 vw::ImageViewBase<ViewT>& right_image) {

      VW_ASSERT(left_image.impl().cols() == right_image.impl().cols() &&
                left_image.impl().rows() == right_image.impl().rows(),
                ArgumentErr() << "multiresolution_correlator: input image dimensions do not agree.\n");

      VW_ASSERT(left_image.impl().channels() == 1 && left_image.impl().planes()==1 &&
                right_image.impl().channels() == 1 && right_image.impl().planes()==1,
                ArgumentErr() << "multiresolution_correlator: does not support multi-channel, multi-plane images.\n");

      int width = left_image.impl().cols();
      int height = left_image.impl().rows();

      //Run the correlator and record how long it takes to run.
      double begin__ = Time();

      std::cout << "Multiresolution Correlator\n";
      ImageView<typename PixelChannelType<typename ViewT::pixel_type>::type> l_image = channels_to_planes(left_image);
      ImageView<typename PixelChannelType<typename ViewT::pixel_type>::type> r_image = channels_to_planes(right_image);        

      std::cout << "\nRunning left-to-right correlation... \n";
      ImageView<PixelDisparity<float> > resultL2R = this->correlate(l_image, r_image, false);
      std::cout << "\nRunning right-to-left correlation... \n";
      ImageView<PixelDisparity<float> > resultR2L = this->correlate(r_image, l_image, true);

      // Cross check the left and right disparity maps
      std::cout << "Cross-correlation consistency check.\n";
      cross_corr_consistency_check(resultL2R, resultR2L,m_crossCorrThreshold);

      // Do subpixel correlation
      ImageView<uint8> l_subpix_image = vw::channel_cast<uint8>(vw::threshold(vw::laplacian_filter(vw::gaussian_filter(l_image,m_slog_width)), 0.0));
      ImageView<uint8> r_subpix_image = vw::channel_cast<uint8>(vw::threshold(vw::laplacian_filter(vw::gaussian_filter(r_image,m_slog_width)), 0.0));
      subpixel_correlation(resultL2R, l_subpix_image, r_subpix_image, m_lKernWidth, m_lKernHeight, m_useHorizSubpixel, m_useVertSubpixel);

      int matched = 0;
      int total = 0;
      int nn = 0;
      for (int j = 0; j < resultL2R.rows(); j++) {
        for (int i = 0; i < resultL2R.cols(); i++) {
          total++;
          if (!(resultL2R(i,j).missing())) {
            matched++;
          }
          nn++;
        }
      } 

      double lapse__ = Time() - begin__;
      if (m_verbose) {
        std::cout << "\tTotal correlation + subpixel took " << lapse__ << " sec";
        double nTries = 2.0 * (m_lMaxV - m_lMinV + 1) * (m_lMaxH - m_lMinH + 1);
        double rate = nTries * width * height / lapse__ / 1.0e6;
        std::cout << "\t(" << rate << " M disparities/second)" << std::endl;
      }

      double score = (100.0 * matched) / total;
      printf("\tCorrelation rate: %6.4f\n\n", score);

      return resultL2R;
    }
  };



}} // namespace vw::stereo

#endif // __VW_STEREO_REFERENCE_CORRELATOR__
