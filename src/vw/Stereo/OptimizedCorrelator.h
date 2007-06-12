#ifndef __VW_STEREO_OPTIMIZED_CORRELATOR__
#define __VW_STEREO_OPTIMIZED_CORRELATOR__

#include <vw/Image/ImageView.h>
#include <vw/Image/Manipulation.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Stereo/Correlate.h>

// Boost
#include <boost/thread/mutex.hpp>

namespace vw { 
namespace stereo {

  // Forward Declarations
  class CorrelationWorkThread;
  
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
    
    template <class PixelT>
    ImageView<PixelDisparity<float> > operator()(vw::ImageView<PixelT>& image0,
                                                 vw::ImageView<PixelT>& image1, 
                                                 bool bit_image = false) {
        
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
      return correlation_2D(l_img, r_img, bit_image);
    }
      
    void set_kernel(int lWidth, int lHeight) {
      m_lKernWidth = lWidth;
      m_lKernHeight = lHeight;
    }
      
    void set_horizontal_range(int lMin, int lMax) {
      m_lMinH = lMin;
      m_lMaxH = lMax;
    }
    
    void set_vertical_range(int lMin, int lMax) {
      m_lMinV = lMin;
      m_lMaxV = lMax;
    }
      
    void set_verbose(int enable) { m_verbose = enable; }
    int get_verbose() const { return m_verbose; }

    void set_horiz_subpixel(int enable) { m_useHorizSubpixel= enable; }
    int get_horiz_subpixel() const { return m_useHorizSubpixel; }
      
    void set_vert_subpixel(int enable) { m_useVertSubpixel= enable; }
    int get_vert_subpixel() const { return m_useVertSubpixel; }
      
    double get_cross_corr_threshold() const { return m_crossCorrThreshold; }
    void set_cross_corr_threshold(double thresh) { m_crossCorrThreshold = thresh; }    
      
    void register_worker_thread(int id, CorrelationWorkThread* child);
    CorrelationWorkThread* get_worker_thread(int id);
      
  private:
      
    std::vector<CorrelationWorkThread*> m_children;
    double m_crossCorrThreshold;
    float m_corrscore_rejection_threshold;
    int m_useHorizSubpixel;
    int m_useVertSubpixel;
    mutable boost::mutex m_mutex;
      
    template <class ChannelT>
    ImageView<PixelDisparity<float> > correlation_2D( ImageView<ChannelT> &left_image,
                                                       ImageView<ChannelT> &right_image, 
                                                       bool bit_image = false);
  };


}}   // namespace vw::stereo

#endif /* __OptimizedCorrelator_h__ */
