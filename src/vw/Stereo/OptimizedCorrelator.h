#ifndef __VW_STEREO_OPTIMIZED_CORRELATOR__
#define __VW_STEREO_OPTIMIZED_CORRELATOR__

#include <vw/Image/ImageView.h>
#include <vw/Image/Manipulation.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Stereo/Correlate.h>

// Boost
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/xtime.hpp>

namespace vw { 
namespace stereo {

  // Forward Declarations
  class CorrelationWorkThreadBase;
  void millisecond_sleep(unsigned long milliseconds);

  class OptimizedCorrelator {
    
    int m_lKernWidth, m_lKernHeight;
    int m_lMinH, m_lMaxH, m_lMinV, m_lMaxV;
    int m_verbose;

  public:
    OptimizedCorrelator();
    OptimizedCorrelator(int minH,	/* left bound disparity search window*/
                        int maxH,	/* right bound disparity search window*/
                        int minV,	/* bottom bound disparity search window */ 
                        int maxV,	/* top bound disparity search window */
                        int kernWidth,	/* size of the kernel */
                        int kernHeight,       
                        int verbose,
                        double crosscorrThreshold,
                        float corrscore_rejection_threshold,
                        int useSubpixelH,
                        int useSubpixelV);
    
    template <class ViewT, class PreProcFilterT>
    ImageView<PixelDisparity<float> > operator()(ImageViewBase<ViewT> const& image0,
                                                 ImageViewBase<ViewT> const& image1,
                                                 PreProcFilterT const& preproc_filter) {
      typedef typename ViewT::pixel_type pixel_type;
      
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
      typename PreProcFilterT::result_type left_image = preproc_filter(channels_to_planes(image0));
      typename PreProcFilterT::result_type right_image = preproc_filter(channels_to_planes(image1));
      return do_correlation(left_image,right_image,preproc_filter.use_bit_image());
    }
      
    void register_worker_thread(int id, CorrelationWorkThreadBase* child) {
      boost::mutex::scoped_lock lock(m_mutex);
      m_children[id] = child;
    }

    CorrelationWorkThreadBase* get_worker_thread(int id) {
      boost::mutex::scoped_lock lock(m_mutex);
      return m_children[id];
    }
      
  private:
      
    template <class ChannelT>
    ImageView<PixelDisparity<float> > do_correlation(ImageView<ChannelT> &left_image, ImageView<ChannelT> &right_image, bool bit_image);

    std::vector<CorrelationWorkThreadBase*> m_children;
    double m_crossCorrThreshold;
    float m_corrscore_rejection_threshold;
    int m_useHorizSubpixel;
    int m_useVertSubpixel;
    mutable boost::mutex m_mutex;
      
  };


}}   // namespace vw::stereo

#endif /* __OptimizedCorrelator_h__ */
