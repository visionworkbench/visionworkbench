// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_STEREO_REFERENCE_CORRELATOR__
#define  __VW_STEREO_REFERENCE_CORRELATOR__

#include <vw/Core/Stopwatch.h>
#include <vw/Stereo/Correlate.h>
#include <vw/Stereo/CorrelateResearch.h>

namespace vw {
namespace stereo {

class ReferenceCorrelator {

  int m_lKernWidth, m_lKernHeight;
  int m_lMinH, m_lMaxH, m_lMinV, m_lMaxV;
  int m_verbose;
  double m_crossCorrThreshold;
  int m_useHorizSubpixel;
  int m_useVertSubpixel;
  int m_do_affine_subpixel;

public:
  ReferenceCorrelator(int minH, /* left bound disparity search window*/
                      int maxH, /* right bound disparity search window*/
                      int minV, /* bottom bound disparity search window */
                      int maxV, /* top bound disparity search window */
                      int kernWidth,    /* size of the kernel */
                      int kernHeight,
                      int verbose,
                      double crosscorrThreshold,
                      int useSubpixelH,
                      int useSubpixelV,
                      bool do_affine_subpixel) {
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
    m_do_affine_subpixel = do_affine_subpixel;
  }

  template <class PixelT>
  ImageView<PixelMask<Vector2f> > correlate(ImageView<PixelT>& left_image,
                                            ImageView<PixelT>& right_image,
                                            bool swap, bool /*use_bit_image*/) {
    ImageView<PixelMask<Vector2f> > result(left_image.cols(), left_image.rows());

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

  template <class ViewT, class PreProcFilterT>
  ImageView<PixelMask<Vector2f> > operator()(ImageViewBase<ViewT> const& image0,
                                             ImageViewBase<ViewT> const& image1,
                                             PreProcFilterT const& preproc_filter) {


    // Check to make sure that image0 and image1 have equal dimensions
    if ((image0.impl().cols() != image1.impl().cols()) ||
        (image0.impl().rows() != image1.impl().rows())) {
      vw_throw( ArgumentErr() << "Primary and secondary image dimensions do not agree!" );
    }

    // Check to make sure that the images are single channel/single plane
    if (!(image0.channels() == 1 && image0.impl().planes() == 1 &&
          image1.channels() == 1 && image1.impl().planes() == 1)) {
      vw_throw( ArgumentErr() << "Both images must be single channel/single plane images!" );
    }

    // Rasterize the views and pass them along to be correlated.
    vw_out() << "Reference Correlator\n";
    typename PreProcFilterT::result_type l_image = preproc_filter(channels_to_planes(image0));
    typename PreProcFilterT::result_type r_image = preproc_filter(channels_to_planes(image1));

    int width = image0.impl().cols();
    int height = image0.impl().rows();

    //Run the correlator and record how long it takes to run.
    Stopwatch timer;
    timer.start();

    vw_out() << "\tRunning left-to-right correlation... " << std::flush;
    ImageView<PixelMask<Vector2f> > resultL2R = this->correlate(l_image, r_image,false,PreProcFilterT::use_bit_image());
    vw_out() << "done.\n\tRunning right-to-left correlation... " << std::flush;
    ImageView<PixelMask<Vector2f> > resultR2L = this->correlate(r_image, l_image,true,PreProcFilterT::use_bit_image());
    vw_out() << "done.\n";

    // Cross check the left and right disparity maps
    vw_out() << "Cross-correlation consistency check.\n";
    cross_corr_consistency_check(resultL2R, resultR2L,m_crossCorrThreshold);

    // Do subpixel correlation
    vw_out() << "Subpixel.\n";
    if (m_do_affine_subpixel)
      subpixel_correlation_affine_2d(resultL2R, l_image, r_image, m_lKernWidth,
                                     m_lKernHeight, m_useHorizSubpixel, m_useVertSubpixel);
    else
      subpixel_correlation_parabola(resultL2R, l_image, r_image, m_lKernWidth,
                                    m_lKernHeight, m_useHorizSubpixel, m_useVertSubpixel);
    vw_out() << "Done Subpixel.\n";

    int matched = 0;
    int total = 0;
    int nn = 0;
    for (int j = 0; j < resultL2R.rows(); j++) {
      for (int i = 0; i < resultL2R.cols(); i++) {
        total++;
        if ( is_valid(resultL2R(i,j)) ) {
          matched++;
        }
        nn++;
      }
    }

    timer.stop();
    double lapse__ = timer.elapsed_seconds();
    if (m_verbose) {
      vw_out() << "\tTotal correlation + subpixel took " << lapse__ << " sec";
      double nTries = 2.0 * (m_lMaxV - m_lMinV + 1) * (m_lMaxH - m_lMinH + 1);
      double rate = nTries * width * height / lapse__ / 1.0e6;
      vw_out() << "\t(" << rate << " M disparities/second)" << std::endl;
    }

    double score = (100.0 * matched) / total;
    vw_out() << "\tCorrelation rate: " << score << "\n\n";

    return resultL2R;
  }
};



}} // namespace vw::stereo

#endif // __VW_STEREO_REFERENCE_CORRELATOR__
