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
  
  template <class ViewT>
  double sum_of_pixels(const ImageViewBase<ViewT> &view) { 
    const ViewT& m_view = view.impl();

    typedef typename ViewT::pixel_accessor pixel_accessor;
    typedef typename ViewT::pixel_type pixel_type;
    typedef typename PixelChannelType<typename ViewT::pixel_type>::type channel_type;

    unsigned num_channels = PixelNumChannels<typename ViewT::pixel_type>::value;

    double accum = 0;

    pixel_accessor plane_iter = m_view.origin();
    for (unsigned p = 0; p < m_view.planes(); p++, plane_iter.next_plane()) { 
      pixel_accessor col_iter = plane_iter;
      for (unsigned i = 0; i < m_view.cols(); i++, col_iter.next_col()) {
        pixel_accessor row_iter = col_iter;
        for (unsigned j = 0; j < m_view.rows(); j++, row_iter.next_row()) {
					pixel_type pix = *row_iter;
          for (unsigned channel = 0; channel < num_channels; channel++) {
            channel_type channel_value = compound_select_channel<channel_type>(pix,channel);
            accum += *row_iter;
          }
        }
      }  
    }
    return accum;
  }
  
  inline double compute_soad_local(const ImageView<double> &Limg,
                                   const ImageView<double> &Rimg,
                                   int left_center_x, 
                                   int left_center_y,
                                   int right_center_x,
                                   int right_center_y,
                                   int half_width_x, 
                                   int half_width_y) {

    return sum_of_pixels(abs(crop(edge_extend(Limg), 
                                  left_center_x - half_width_x, 
                                  left_center_y - half_width_y, 
                                  half_width_x, 
                                  half_width_y) - 
                             crop(edge_extend(Rimg), 
                                  right_center_x - half_width_x, 
                                  right_center_y - half_width_y, 
                                  half_width_x, 
                                  half_width_y)));
  }

class MultiresolutionCorrelator {

  int m_lKernWidth, m_lKernHeight;
  int m_lMinH, m_lMaxH, m_lMinV, m_lMaxV;
  int m_verbose;
  double m_crossCorrThreshold;
  int m_useHorizSubpixel;
  int m_useVertSubpixel;
  double m_slog_width;

public:
  MultiresolutionCorrelator(int minH,	/* left bound disparity search window*/
                            int maxH,	/* right bound disparity search window*/
                            int minV,	/* bottom bound disparity search window */ 
                            int maxV,	/* top bound disparity search window */
                            int kernWidth,	/* size of the kernel */
                            int kernHeight,       
                            int verbose,
                            double crosscorrThreshold,
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
    m_useHorizSubpixel = useSubpixelH;
    m_useVertSubpixel = useSubpixelV;
  }

//   template <class PixelT>
//   ImageView<PixelDisparity<float> > correlate(vw::ImageView<PixelT>& left_image,
//                                               vw::ImageView<PixelT>& right_image,
//                                               bool swap) {  
    
//     ImageView<PixelDisparity<float> > result(left_image.cols(), left_image.rows());

//     for (int i = 0; i < left_image.cols(); i++) {
//       for (int j = 0; j < left_image.rows(); j++) {
//         if (swap) {
//           result(i,j) = compute_disparity(left_image, right_image, i, j, 
//                                           m_lKernWidth, m_lKernHeight,
//                                           -m_lMaxH, -m_lMinH, -m_lMaxV, -m_lMinV);
//         } else {
//           result(i,j) = compute_disparity(left_image, right_image, i, j, 
//                                           m_lKernWidth, m_lKernHeight,
//                                           m_lMinH, m_lMaxH, m_lMinV, m_lMaxV);
//         }
//       }
//     }
//     return result;
//   }

  template <class PixelT>
  ImageView<PixelDisparity<float> > correlate(vw::ImageView<PixelT>& left_image,
                                              vw::ImageView<PixelT>& right_image,
                                              bool swap) {  

    set_debug_level(InfoMessage+1);
    
    // Compute the number of image pyramid levels
    const int min_dimension = 256;
    int dimension = std::min(left_image.cols(), left_image.rows());
    int levels = 1;
    while (dimension > min_dimension) { levels++; dimension /= 2; }
    std::cout << "Building image pyramid with " << levels << " levels. \n";

    // Build the image pyramid
    std::vector<ImageView<float> > left_pyramid(levels), right_pyramid(levels);
    std::vector<ImageView<PixelDisparity<float> > > disparity_pyramid(levels);
    left_pyramid[0] = channels_to_planes(left_image);
    right_pyramid[0] = channels_to_planes(right_image);
    
    // Apply the SLOG filter to each level.  Should we change the SLOG width here? Probably...
    const float slog_width = m_slog_width;
    for (int n = 1; n < levels; ++n) {
      std::cout << "\tProcessing level " << n << "        \n";
      left_pyramid[n] = subsample_by_two(left_pyramid[n-1]);
      right_pyramid[n] = subsample_by_two(right_pyramid[n-1]);
      std::ostringstream current_level;
      current_level << n;
      std::ostringstream previous_level;
      previous_level << n-1;
      write_image("test/level" + current_level.str() + "-L.jpg", left_pyramid[n]);
      write_image("test/level" + current_level.str() + "-R.jpg", right_pyramid[n]);
      left_pyramid[n-1] = vw::threshold(vw::laplacian_filter(vw::gaussian_filter(left_pyramid[n-1],slog_width)), 0.0);
      right_pyramid[n-1] = vw::threshold(vw::laplacian_filter(vw::gaussian_filter(right_pyramid[n-1],slog_width)), 0.0);
      write_image("test/slog" + previous_level.str() + "-L.jpg", left_pyramid[n-1]);
      write_image("test/slog" + previous_level.str() + "-R.jpg", right_pyramid[n-1]);
    }
    left_pyramid[levels-1] = vw::threshold(vw::laplacian_filter(vw::gaussian_filter(left_pyramid[levels-1],slog_width)), 0.0);
    right_pyramid[levels-1] = vw::threshold(vw::laplacian_filter(vw::gaussian_filter(right_pyramid[levels-1],slog_width)), 0.0);

    int h_min = m_lMinH / (levels);
    int h_max = m_lMaxH / (levels);
    int v_min = m_lMinV / (levels);
    int v_max = m_lMaxV / (levels);
    int h_kern = m_lKernWidth;
    int v_kern = m_lKernHeight;
    
    // Get things started with a full correlation.
    OptimizedCorrelator correlator(h_min, h_max,
                                  v_min, v_max,
                                  h_kern, v_kern,
                                  true,          // verbose
                                  m_crossCorrThreshold,
                                  false, false); // no subpixel for now
    disparity_pyramid[levels-1] = correlator(left_pyramid[levels-1], 
                                             right_pyramid[levels-1], 
                                             true);
    
    // Print out the range of disparity values
    double min_h_disp, min_v_disp, max_h_disp, max_v_disp;
    disparity::get_disparity_range(disparity_pyramid[levels-1], min_h_disp, max_h_disp, min_v_disp, max_v_disp,true);

    // Print out the disparity map at the lowest resolution
    write_image( "test/DH-0.jpg", normalize(clamp(select_channel(disparity_pyramid[levels-1],0), min_h_disp, max_h_disp)));
    write_image( "test/DV-0.jpg", normalize(clamp(select_channel(disparity_pyramid[levels-1],1), min_v_disp, max_v_disp)));

    // Now, clean up the disparity map by rejecting outliers 
    unsigned rm_half_kernel = 5 / (int)powf(2,levels) + 1;
    unsigned rm_min_matches_percent = 60 / levels;
    double rm_threshold = 3;
    for(int nn=0; nn < 1; nn++) {  // Run the clean up routine three times
      disparity::clean_up(disparity_pyramid[levels-1],
                          rm_half_kernel, 
                          rm_half_kernel,
                          rm_min_matches_percent,
                          rm_threshold,
                          true);
    }
    std::cout << "\tRemoving solitary pixels [20x20 window, 20% threshold]\n";
    disparity::remove_outliers(disparity_pyramid[levels-1], 20, 20, 20, 200, true);

    // These settings should move from being hard coded to being user configurable
    disparity::sparse_disparity_filter(disparity_pyramid[levels-1], 20, 0.1);

    // Print out the disparity map at the lowest resolution
    write_image( "test/DH-0-filt.jpg", normalize(clamp(select_channel(disparity_pyramid[levels-1],0), min_h_disp, max_h_disp)));
    write_image( "test/DV-0-filt.jpg", normalize(clamp(select_channel(disparity_pyramid[levels-1],1), min_v_disp, max_v_disp)));

    // Print out the range of disparity values
    disparity::get_disparity_range(disparity_pyramid[levels-1], min_h_disp, max_h_disp, min_v_disp, max_v_disp,true);

    // Refined the disparity map by searching in the local region where the last good disparity value was found
    std::cout << "Refining disparity map.\n";
    for (int n = levels - 2; n >=0; --n) {
      std::cout << "Level " << n << "\n";
      // FIXME: Remove the clamp hack here that takes care of missing pixel values
      disparity_pyramid[n] = resample(disparity_pyramid[n+1], 2,
                                      (int)left_pyramid[n].cols(),
                                      (int)left_pyramid[n].rows(), 
                                      ZeroEdgeExtend(),
                                      NearestPixelInterpolation());
      
      // Print out the disparity map at the lowest resolution
      std::ostringstream current_level;
      current_level << n;
      disparity::get_disparity_range(disparity_pyramid[n], min_h_disp, max_h_disp, min_v_disp, max_v_disp,true);
      write_image( "test/disp-resample-H-" + current_level.str() + "-filt.jpg", normalize(clamp(select_channel(disparity_pyramid[n],0), min_h_disp, max_h_disp)));
      write_image( "test/disp-resample-V-" + current_level.str() + "-filt.jpg", normalize(clamp(select_channel(disparity_pyramid[n],1), min_v_disp, max_v_disp)));

//       ImageView<double> orig = copy(disparity_pyramid[n]);
    
      double best_score, best_j_disp, best_i_disp;
      int h_kern_halfwidth  = h_kern / 2;
      int v_kern_halfwidth  = v_kern / 2;
      for (int i=h_kern; i<disparity_pyramid[n].cols()-h_kern; ++i) {
        for (int j=v_kern; j<disparity_pyramid[n].rows()-v_kern; ++j) {
          
          if (!disparity_pyramid[n](i,j).missing()) {
            disparity_pyramid[n](i,j).h() *= 2;
            //            disparity_pyramid[n](i,j).v() *= 2;
            best_score = 1e99;
            for (int i_disp = -5; i_disp < 5; ++i_disp) {
              //              for (int j_disp = -2; j_disp < 2; ++j_disp) {
                //              std::cout << "\t" << i << " " << j << " " << disparity_pyramid[n](i,j,0) << " " << disparity_pyramid[n](i,j,1) << " " << (int)(i+i_disp+disparity_pyramid[n](i,j,0)) << " " << (int)(j+j_disp+disparity_pyramid[levels-2](i,j,1)) << "\n";
//                 double current_score = compute_soad_local(left_pyramid[n],
//                                                           right_pyramid[n],
//                                                           i, j,
//                                                           (int)(i+i_disp+disparity_pyramid[n](i,j).h()),
//                                                           int(j),
// //                                                           (int)(j+j_disp+disparity_pyramid[n](i,j).v()),
//                                                           h_kern_halfwidth, v_kern_halfwidth);
                double current_score = compute_soad(&(left_pyramid[n](0,0)), &(right_pyramid[n](0,0)),
                                                    j,i,
                                                    int(i_disp+disparity_pyramid[n](i,j).h()), 
                                                    j,
                                                    //                                                    int(j_disp+disparity_pyramid[n](i,j).v()),
                                                    h_kern, v_kern,
                                                    left_pyramid[n].cols(), 
                                                    left_pyramid[n].rows());
                //  std::cout << "Score: (" << i_disp << ", " << j_disp << ") : " << current_score << "   (best: " << best_score << ")\n";
                if (current_score < best_score) {
                  best_score = current_score;
                  best_i_disp = i_disp;
                  //                  best_j_disp = j_disp;
                }
                
                //              }
            }
            // std::cout << "Best:  " << best_i_disp << ", " << best_j_disp << ".\n";
            disparity_pyramid[n](i,j).h() += best_i_disp;
            //            disparity_pyramid[n](i,j).v() += best_j_disp;
          }
        }
      }

      disparity::get_disparity_range(disparity_pyramid[n], min_h_disp, max_h_disp, min_v_disp, max_v_disp,true);
      write_image( "test/disp-refined-H-" + current_level.str() + "-filt.jpg", normalize(clamp(select_channel(disparity_pyramid[n],0), min_h_disp, max_h_disp)));
      write_image( "test/disp-refined-V-" + current_level.str() + "-filt.jpg", normalize(clamp(select_channel(disparity_pyramid[n],1), min_v_disp, max_v_disp)));

//       ostringstream ostream3;
//       ostream3 << "disparity-diff-" << n << "*.jpg";
//       write_image(ostream3.str(), normalize(clamp(disparity_pyramid[n]-orig, -100, 100)));
    }
    return disparity_pyramid[0];
  }

  template <class PixelT>
  ImageView<PixelDisparity<float> > operator()(vw::ImageView<PixelT>& left_image,
                                               vw::ImageView<PixelT>& right_image, 
                                               bool bit_image = false) {

    VW_ASSERT(left_image.cols() == right_image.cols() &&
              left_image.rows() == right_image.rows(),
              ArgumentErr() << "subpixel_correlation: input image dimensions do not agree.\n");

    int width = left_image.cols();
    int height = left_image.rows();

    //Run the correlator and record how long it takes to run.
    double begin__ = Time();

    // Ask the worker threads for the actual results of the disparity correlation
    std::cout << "Multiresolution Correlator\n";
    ImageView<typename PixelChannelType<PixelT>::type> l_image = channels_to_planes(left_image);
    ImageView<typename PixelChannelType<PixelT>::type> r_image = channels_to_planes(right_image);        

//     std::cout << "\tRunning left-to-right correlation... " << std::flush;
    ImageView<PixelDisparity<float> > resultL2R = this->correlate(l_image, r_image,false);
//     std::cout << "done.\n\tRunning right-to-left correlation... " << std::flush;
//     ImageView<PixelDisparity<float> > resultR2L = this->correlate(r_image, l_image,true);
//     std::cout << "done.\n";

//     // Cross check the left and right disparity maps
//     std::cout << "Cross-correlation consistency check.\n";
//     cross_corr_consistency_check(resultL2R, resultR2L,m_crossCorrThreshold);

//     // Do subpixel correlation
//     std::cout << "Subpixel.\n";
//     subpixel_correlation(resultL2R, l_image, r_image, m_lKernWidth, m_lKernHeight, m_useHorizSubpixel, m_useVertSubpixel);
//     std::cout << "Done Subpixel.\n";

//     int matched = 0;
//     int total = 0;
//     int nn = 0;
//     for (int j = 0; j < resultL2R.rows(); j++) {
//       for (int i = 0; i < resultL2R.cols(); i++) {
//         total++;
//         if (!(resultL2R(i,j).missing())) {
//           matched++;
//         }
//         nn++;
//       }
//     } 

//     double lapse__ = Time() - begin__;
//     if (m_verbose) {
//       std::cout << "\tTotal correlation + subpixel took " << lapse__ << " sec";
//       double nTries = 2.0 * (m_lMaxV - m_lMinV + 1) * (m_lMaxH - m_lMinH + 1);
//       double rate = nTries * width * height / lapse__ / 1.0e6;
//       std::cout << "\t(" << rate << " M disparities/second)" << std::endl;
//     }

//     double score = (100.0 * matched) / total;
//     printf("\tCorrelation rate: %6.4f\n\n", score);

    return resultL2R;
  }
};



}} // namespace vw::stereo

#endif // __VW_STEREO_REFERENCE_CORRELATOR__
