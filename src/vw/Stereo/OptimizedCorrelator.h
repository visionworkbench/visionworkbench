#ifndef __VW_STEREO_OPTIMIZED_CORRELATOR_H__
#define __VW_STEREO_OPTIMIZED_CORRELATOR_H__

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
                        int useSubpixelH,
                        int useSubpixelV);
    
    template <class PixelT>
    ImageView<PixelDisparity<float> > operator()(vw::ImageView<PixelT>& image0,
                                                 vw::ImageView<PixelT>& image1, 
                                                 bool bit_image = false) {
        
      // Check to make sure that image0 and image1 have equal dimensions 
      if ((image0.cols() != image1.cols()) ||
          (image0.rows() != image1.rows())) {
        throw ArgumentErr() << "Primary and secondary image dimensions do not agree\n";
      }
        
      // Check to make sure that the images are single channel/single plane
      if (!(image0.channels() == 1 && image0.planes() == 1 &&
            image1.channels() == 1 && image1.planes() == 1)) {
        throw ArgumentErr() << "Both images must be single channel/single plane images\n";
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
      
    typedef struct soadStruct {
      double hDisp;	/* disparity in x for the best match */
      double vDisp;	/* disparity in y for the best match */
      double best;             	   /* soad at the best match */
    } soad;
      
    void register_worker_thread(int id, CorrelationWorkThread* child);
    CorrelationWorkThread* get_worker_thread(int id);
      
  private:
      
    std::vector<CorrelationWorkThread*> m_children;
    double m_crossCorrThreshold;
    int m_useHorizSubpixel;
    int m_useVertSubpixel;
    mutable boost::mutex m_mutex;
      
    template <class ChannelT>
    ImageView<PixelDisparity<float> > correlation_2D( ImageView<ChannelT> &left_image,
                                                       ImageView<ChannelT> &right_image, 
                                                       bool bit_image = false);
  };

  /// Work thread class that does the actual correlation
  struct CorrelationWorkThread {
    CorrelationWorkThread();
    CorrelationWorkThread(const CorrelationWorkThread& other);
    const CorrelationWorkThread& operator=(const CorrelationWorkThread& other);
    CorrelationWorkThread(OptimizedCorrelator* parent,
                          int id,
                          int min_h, int max_h,
                          int min_v, int max_v,
                          int image_width, int image_height,
                          int kern_width, int kern_height,
                          float* left_image, float* right_image);
    CorrelationWorkThread(OptimizedCorrelator* parent, 
                          int id,
                          int min_h, int max_h,
                          int min_v, int max_v,
                          int image_width, int image_height,
                          int kern_width, int kern_height,
                          unsigned char* left_bit_image, unsigned char* right_bit_image);

    /// Call this operator to start the correlation operation.
    void operator() ();

    /// Returns a number between 0.0 and 100.0
    double progress_percent();
    std::string correlation_rate_string();
    std::string progress_string();
    bool is_done();
    void terminate();
    ImageView<PixelDisparity<float> > result() { return m_result; }

  private:
    OptimizedCorrelator::soad *fast2Dcorr_optimized(
                                                   int minDisp,	/* left bound disparity search */
                                                   int maxDisp,	/* right bound disparity search */
                                                   int topDisp,	/* top bound disparity search window */
                                                   int btmDisp,	/* bottom bound disparity search window */ 
                                                   int height,	/* image height */
                                                   int width,	/* image width */
                                                   int kernel_height, 
                                                   int kernel_width,
                                                   const unsigned char *Rimg,	/* reference image fixed */
                                                   const unsigned char *Simg	/* searched image sliding */
                                                   );
      
    OptimizedCorrelator::soad * fast2Dcorr(int minDisp,	/* left bound disparity search */
                                          int maxDisp,	/* right bound disparity search */
                                          int topDisp,	/* top bound disparity search window */
                                          int btmDisp,	/* bottom bound disparity search window */ 
                                          int height,	/* image height */
                                          int width,	/* image width */
                                          int vKern,  /* kernel height */
                                          int hKern,  /* kernel width */
                                          const float *Rimg,	/* reference image fixed */
                                          const float *Simg	/* searched image sliding */
                                          );

    void set_correlation_rate_string(std::string str);
    void set_progress_percent(double percent);
    void set_progress_string(std::string str);
    void set_done(bool val);
      
    OptimizedCorrelator* m_parent;
    int m_id;
    ImageView<PixelDisparity<float> > m_result;
      
    int m_min_h, m_max_h;
    int m_min_v, m_max_v;
    int m_image_height, m_image_width;
    int m_kern_height, m_kern_width;
    unsigned char* m_bit_image_left;
    unsigned char* m_bit_image_right;
    float* m_image_left;
    float* m_image_right;
      
    std::string m_correlation_rate_string;
    double m_progress_percent;
    std::string m_progress_string;
    bool m_done;
    bool m_should_terminate;
    mutable boost::mutex m_mutex;
    char temp_progress[100];
  };

}}   // namespace vw::stereo

#endif /* __OptimizedCorrelator_h__ */
