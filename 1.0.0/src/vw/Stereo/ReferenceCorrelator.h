#ifndef __VW_STEREO_REFERENCE_CORRELATOR__
#define  __VW_STEREO_REFERENCE_CORRELATOR__
#include <vw/Stereo/Correlate.h>

namespace vw {
namespace stereo {

class ReferenceCorrelator {
  
  int m_lKernWidth, m_lKernHeight;
  int m_lMinH, m_lMaxH, m_lMinV, m_lMaxV;
  int m_verbose;
  double m_crossCorrThreshold;
  int m_useHorizSubpixel;
  int m_useVertSubpixel;

public:
  ReferenceCorrelator(int minH,	/* left bound disparity search window*/
                      int maxH,	/* right bound disparity search window*/
                      int minV,	/* bottom bound disparity search window */ 
                      int maxV,	/* top bound disparity search window */
                      int kernWidth,	/* size of the kernel */
                      int kernHeight,       
                      int verbose,
                      double crosscorrThreshold,
                      int useSubpixelH,
                      int useSubpixelV) {
    m_lKernWidth = kernWidth;
    m_lKernHeight = kernHeight;
    m_lMinH = minH;
    m_lMaxH = maxH;
    m_lMinV = minV;
    m_lMaxV = maxV;  
    m_verbose = verbose;
    
    m_crossCorrThreshold = crosscorrThreshold;
    m_useHorizSubpixel = useSubpixelH;
    m_useVertSubpixel = useSubpixelV;
  }

  template <class PixelT>
  ImageView<PixelDisparity<float> > correlate(vw::ImageView<PixelT>& left_image,
                                              vw::ImageView<PixelT>& right_image,
                                              bool swap) {  
    
    ImageView<PixelDisparity<float> > result(left_image.cols(), left_image.rows());

    for (int i = 0; i < left_image.cols(); i++) {
      for (int j = 0; j < left_image.rows(); j++) {
        if (swap) {
          result(i,j) = compute_disparity(left_image, right_image, i, j, 
                                          m_lKernWidth, m_lKernHeight,
                                          -m_lMaxH, -m_lMinH, -m_lMaxV, -m_lMinV);
        } else {
          result(i,j) = compute_disparity(left_image, right_image, i, j, 
                                          m_lKernWidth, m_lKernHeight,
                                          m_lMinH, m_lMaxH, m_lMinV, m_lMaxV);
        }
      }
    }
    return result;
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
    std::cout << "Reference Correlator\n";
    ImageView<typename PixelChannelType<PixelT>::type> l_image = channels_to_planes(left_image);
    ImageView<typename PixelChannelType<PixelT>::type> r_image = channels_to_planes(right_image);        

    std::cout << "\tRunning left-to-right correlation... " << std::flush;
    ImageView<PixelDisparity<float> > resultL2R = this->correlate(l_image, r_image,false);
    std::cout << "done.\n\tRunning right-to-left correlation... " << std::flush;
    ImageView<PixelDisparity<float> > resultR2L = this->correlate(r_image, l_image,true);
    std::cout << "done.\n";

    // Cross check the left and right disparity maps
    std::cout << "Cross-correlation consistency check.\n";
    cross_corr_consistency_check(resultL2R, resultR2L,m_crossCorrThreshold);

    // Do subpixel correlation
    std::cout << "Subpixel.\n";
    subpixel_correlation(resultL2R, l_image, r_image, m_lKernWidth, m_lKernHeight, m_useHorizSubpixel, m_useVertSubpixel);
    std::cout << "Done Subpixel.\n";

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
