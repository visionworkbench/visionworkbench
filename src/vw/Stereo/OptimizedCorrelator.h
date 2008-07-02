#ifndef __VW_STEREO_OPTIMIZED_CORRELATOR__
#define __VW_STEREO_OPTIMIZED_CORRELATOR__

#include <vw/Core/Stopwatch.h>
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

  // These functors allow us to specialize the behavior of the image
  // differencing operation, which is part of measuring the sum of
  // absolute difference (SOAD) between two images.  For 8-bit slog
  // images, we take the xor (^) of the two images.
  struct StandardAbsDifferenceFunc {
    static inline float absdiff (const float val1, const float val2) { return fabs(val1-val2); }
    static inline double absdiff (const double val1, const double val2) { return fabs(val1-val2); }
    static inline uint8 absdiff (const uint8 val1, const uint8 val2) { return val1 > val2 ? (val1-val2) : (val2-val1); }
    static inline uint16 absdiff (const uint16 val1, const uint16 val2) { return val1 > val2 ? (val1-val2) : (val2-val1); }
    static inline uint32 absdiff (const uint32 val1, const uint32 val2) { return val1 > val2 ? (val1-val2) : (val2-val1); }
    static inline int8 absdiff (const int8 val1, const int8 val2) { return abs(val1-val2); }
    static inline int16 absdiff (const int16 val1, const int16 val2) { return abs(val1-val2); }
    static inline int32 absdiff (const int32 val1, const int32 val2) { return abs(val1-val2); }
  };

  struct BitAbsDifferenceFunc {
    static inline uint8 absdiff (const uint8 val1, const uint8 val2) { return val1^val2; }
  };

  /// Given a type, these traits classes help to determine a suitable
  /// working type for accumulation operations in the correlator
  template <class T> struct CorrelatorAccumulatorType {};
  template <> struct CorrelatorAccumulatorType<vw::uint8>   { typedef vw::uint16  type; };
  template <> struct CorrelatorAccumulatorType<vw::int8>    { typedef vw::uint16  type; };
  template <> struct CorrelatorAccumulatorType<vw::uint16>  { typedef vw::uint32  type; };
  template <> struct CorrelatorAccumulatorType<vw::int16>   { typedef vw::uint32  type; };
  template <> struct CorrelatorAccumulatorType<vw::uint32>  { typedef vw::uint64  type; };
  template <> struct CorrelatorAccumulatorType<vw::int32>   { typedef vw::uint64  type; };
  template <> struct CorrelatorAccumulatorType<vw::float32> { typedef vw::float32 type; };
  template <> struct CorrelatorAccumulatorType<vw::float64> { typedef vw::float64 type; };


  // Threading classes
  class CorrelationWorkThreadBase {
  protected:
    ImageView<PixelDisparity<float> > m_result;

    std::string m_correlation_rate_string;
    double m_progress_percent;
    std::string m_progress_string;
    bool m_done;
    mutable Mutex m_mutex;
    char temp_progress[100];
  public:

    CorrelationWorkThreadBase() : m_progress_percent(0), m_done(false) {}
    virtual ~CorrelationWorkThreadBase() {}
    virtual void operator()() = 0;

    virtual ImageView<PixelDisparity<float> > result() { return m_result; }
  
    std::string correlation_rate_string() const {
      Mutex::Lock lock(m_mutex);
      return m_correlation_rate_string;
    }
  
    std::string progress_string() const {
      Mutex::Lock lock(m_mutex);
      return m_progress_string;
    }
  
    /// Returns a number between 0.0 and 100.0
    double progress_percent() { return m_progress_percent; }
  
    bool is_done() const {
      return m_done;
    }

    void set_correlation_rate_string(std::string str) {
      Mutex::Lock lock(m_mutex);
      m_correlation_rate_string = str;
    }
  
    void set_progress_string(std::string str) {
      Mutex::Lock lock(m_mutex);
      m_progress_string = str;
    }
  
    void set_done(bool val) {
      m_done = val;
    }
  };

  /// The correlation work thread is a functor that does the actual
  /// correlation processing.  This functor can be used to easily spin
  /// off multiple threads of correlation on a multiprocessor machine.
  template <class ViewT, class AbsDifferenceT>
  class CorrelationWorkThread : public CorrelationWorkThreadBase {

    // Handy typedefs
    typedef typename ViewT::pixel_type pixel_type;
    typedef typename PixelChannelType<pixel_type>::type channel_type;
    
    // This is used as an intermediate format for the correlator
    typedef struct soadStruct {
      double hDisp; // disparity in x for the best match 
      double vDisp; // disparity in y for the best match 
      double best;  // soad at the best match 
    } soad;
  
    int m_min_h, m_max_h;
    int m_min_v, m_max_v;
    int m_kern_height, m_kern_width;
    float m_corrscore_rejection_threshold;

    ViewT m_left_image;
    ViewT m_right_image;

    template <class ChannelT>
    soad *fast2Dcorr_optimized(int minDisp,	/* left bound disparity search */
                               int maxDisp,	/* right bound disparity search */
                               int topDisp,	/* top bound disparity search window */
                               int btmDisp,	/* bottom bound disparity search window */ 
                               int height,	/* image height */
                               int width,	/* image width */
                               int vKern,  /* kernel height */
                               int hKern,  /* kernel width */
                               float corrscore_rejection_threshold, /* correlation score fitness score (1.0 disables, 1.5-2.0 is a good value) */
                               const ChannelT* Rimg,	/* reference image fixed */
                               const ChannelT* Simg	  /* searched image sliding */
                               ) {
      
      typedef typename CorrelatorAccumulatorType<ChannelT>::type accumulator_type;

      // Start a timer
      vw::Stopwatch timer;
      timer.start();

      int numCorrTries = (btmDisp - topDisp + 1) * (maxDisp - minDisp + 1);

      ChannelT *diff = new ChannelT[width*height]; // buffer containing results of substraction of L/R images
      VW_ASSERT( diff, NullPtrErr() << "cannot allocate the correlator's difference buffer!" );

      accumulator_type *cSum = new accumulator_type[width]; // sum per column (kernel hight)
      VW_ASSERT( cSum, NullPtrErr() << "cannot allocate the correlator's column sum buffer!" );

      struct local_result {
        accumulator_type best, worst;
        uint16 hvdsp;
      };
      local_result *result_buf = new local_result[width * height]; // correlation result buffer
      VW_ASSERT( result_buf, NullPtrErr() << "cannot allocate the correlator's local result buffer!" );

      // fill the result_buf with default values
      for (int nn = 0; nn < (width*height); ++nn) {
        result_buf[nn].worst = ScalarTypeLimits<accumulator_type>::lowest();
        result_buf[nn].best = ScalarTypeLimits<accumulator_type>::highest();
        result_buf[nn].hvdsp = ScalarTypeLimits<accumulator_type>::highest();
      }

      // for each allowed disparity...
      for( int dsy=topDisp; dsy<=btmDisp; ++dsy ) {
        for( int ds=minDisp; ds<=maxDisp; ++ds ) {
          set_progress_string(str(boost::format("V: [%1%,%2%] H: [%3%,%4%] processing %5% %6%")
                                  % topDisp % btmDisp % minDisp % maxDisp % dsy % ds));
      
          uint16 ds_combined = (ds-minDisp) + (dsy - topDisp) * (maxDisp - minDisp + 1);
          
          // compute the region of correlation
          int yStart = (dsy<0) ? (-dsy) : 0;
          int xStart = (ds<0) ? (-ds) : 0;
          int yEnd = yStart + height - abs(dsy) - (vKern-1);
          int xEnd = xStart + width - abs(ds) - (hKern-1);
      
          // compute the difference buffer
          ChannelT *diff_row = diff + yStart*width + xStart;
          const ChannelT *rimg_row = Rimg + yStart*width + xStart;
          const ChannelT *simg_row = Simg + (yStart+dsy)*width + (xStart+ds);
          for( int j=yEnd-yStart+(vKern-1); j; --j ) {
            ChannelT *diff_ptr = diff_row;
            const ChannelT *rimg_ptr = rimg_row;
            const ChannelT *simg_ptr = simg_row;
            for( int i=xEnd-xStart+(hKern-1); i; --i ) {
              *(diff_ptr++) = AbsDifferenceT::absdiff( (*(rimg_ptr++)), (*(simg_ptr++)) );
            }
            diff_row += width;
            rimg_row += width;
            simg_row += width;
          }

          // seed the column sum buffer
          for( int i=0; i<width; ++i ) cSum[i] = 0;
          diff_row = diff + yStart*width + xStart;
          accumulator_type *csum_row = cSum + xStart;
          for( int j=vKern; j; --j ) {
            ChannelT *diff_ptr = diff_row;
            accumulator_type *csum_ptr = csum_row;
            for( int i=xEnd-xStart+(hKern-1); i; --i ) {
              *(csum_ptr++) += *(diff_ptr++);
            }
            diff_row += width;
          }

          // perform the correlation
          local_result *result_row = result_buf+(xStart+hKern/2)+(yStart+vKern/2)*width;
          for( int j=yStart; ; ++j, result_row+=width ) {
            local_result *result_ptr = result_row;

            // seed the row sum (i.e. soad)
            typename AccumulatorType<accumulator_type>::type rsum=0;
            accumulator_type *csum_ptr = cSum + xStart;
            for( int i=hKern; i; --i )
              rsum += *(csum_ptr++);

            // correlate the row
            accumulator_type *csum_tail = cSum + xStart;
            accumulator_type *csum_head = csum_tail + hKern;
            for( int i=xEnd-xStart; i; --i, ++result_ptr ) {
              // store the result if better
              if( rsum < result_ptr->best ) {
                result_ptr->best = rsum;
                result_ptr->hvdsp = ds_combined;
              }

              if( rsum > result_ptr->worst ) 
                result_ptr->worst = rsum;

              // update the row sum
              rsum += *(csum_head++) - *(csum_tail++);
            }

            // break if this is the last row
            if( j+1 == yEnd ) break;

            // update the column sum
            csum_ptr = cSum + xStart;
            ChannelT *diff_tail = diff + xStart + j*width;
            ChannelT *diff_head = diff_tail + vKern*width;
            for( int i=xEnd-xStart+(hKern-1); i; --i ) {
              *(csum_ptr++) += *(diff_head++) - *(diff_tail++);
            }
          }
        }
      }
      delete[] diff;
      delete[] cSum;

      soad *result;      // soad buffer to be returned
      result = (soad *)malloc (width*height*sizeof(soad));
      VW_ASSERT( result, NullPtrErr() << "cannot allocate the correlator's SOAD buffer!" );
  
      // convert from the local result buffer to the return format and
      // reject any pixels that do not have a sufficiently distincive
      // peak in the correlation space.
      for( int j = 0; j < height; ++j) {
        for( int i = 0; i < width; ++i) {
          float ratio = float(result_buf[j*width+i].worst+1) / float(result_buf[j*width+i].best+1);
          int nn = width*j+i;
          if( result_buf[nn].best == USHRT_MAX || ratio <= m_corrscore_rejection_threshold ) {
            result[nn].best = VW_STEREO_MISSING_PIXEL;
            result[nn].hDisp = VW_STEREO_MISSING_PIXEL;
            result[nn].vDisp = VW_STEREO_MISSING_PIXEL;
          } else {
            result[nn].best = result_buf[nn].best;
            int hvDisp = result_buf[nn].hvdsp;
            int hDisp = hvDisp % (maxDisp-minDisp+1) + minDisp;
            int vDisp = hvDisp / (maxDisp-minDisp+1) + topDisp;
            result[nn].hDisp = hDisp;
            result[nn].vDisp = vDisp;
          }
        }
      }
      delete[] result_buf;

      set_progress_string(str(boost::format("V: [%1%,%2%] H: [%3%,%4%] processed successfully.")
                              % topDisp % btmDisp % minDisp % maxDisp));
  
      timer.stop();
      double duration = timer.elapsed_seconds();
      double mdisp_per_sec= ((double)numCorrTries * width * height)/duration/1e6;
      set_correlation_rate_string(str(boost::format("%1% seconds (%2% M disparities/second)")
                                      % (int)duration % mdisp_per_sec));

      return(result);
    }

  public:

    CorrelationWorkThread(int min_h, int max_h,
                          int min_v, int max_v,
                          int kern_width, int kern_height,
                          float corrscore_rejection_threshold,
                          ViewT const& left_image, 
                          ViewT const& right_image) {
      m_min_h = min_h;
      m_max_h = max_h;
      m_min_v = min_v;
      m_max_v = max_v;
      m_kern_height = kern_height;
      m_kern_width = kern_width;
      m_left_image = left_image;
      m_right_image = right_image;
      m_corrscore_rejection_threshold = corrscore_rejection_threshold;
    }

  public:
    virtual ~CorrelationWorkThread() {}
  
    /// Call this operator to start the correlation operation.
    virtual void operator()() {    
      int width = m_left_image.cols();
      int height = m_left_image.rows();

      soad* result = fast2Dcorr_optimized(m_min_h, m_max_h, m_min_v, m_max_v, height, width, m_kern_height, m_kern_width, m_corrscore_rejection_threshold, &(m_left_image(0,0)), &(m_right_image(0,0)));
    
      m_result.set_size(width, height);
      for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
          if (result[j*width+i].hDisp == VW_STEREO_MISSING_PIXEL) {
            m_result(i,j) = PixelDisparity<float>();  // Default constructor creates a missing pixel
          } else {
            m_result(i,j) = PixelDisparity<float>(result[j*width+i].hDisp,
                                                  result[j*width+i].vDisp);
          }
        }
      }
    
      delete [] result;  
      set_done(true);
    }

  };

  class OptimizedCorrelator {
    
    int m_lKernWidth, m_lKernHeight;
    int m_lMinH, m_lMaxH, m_lMinV, m_lMaxV;
    int m_verbose;
    double m_crossCorrThreshold;
    float m_corrscore_rejection_threshold;
    int m_useHorizSubpixel, m_useVertSubpixel;
    mutable Mutex m_mutex;

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

      typedef typename CompoundChannelType<typename ViewT::pixel_type>::type channel_type;

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

      typename PreProcFilterT::result_type left_image = preproc_filter(channels_to_planes(image0));
      typename PreProcFilterT::result_type right_image = preproc_filter(channels_to_planes(image1));  

      //Run the correlator and record how long it takes to run.
      Stopwatch timer;
      timer.start();
      
      // Configure the workers
      boost::shared_ptr<CorrelationWorkThreadBase> worker0, worker1;
      if( preproc_filter.use_bit_image() ) {
        worker0.reset( new CorrelationWorkThread<ImageView<channel_type>, BitAbsDifferenceFunc>(-m_lMaxH, -m_lMinH, -m_lMaxV, -m_lMinV, m_lKernWidth, m_lKernHeight, m_corrscore_rejection_threshold, right_image, left_image) );
        worker1.reset( new CorrelationWorkThread<ImageView<channel_type>, BitAbsDifferenceFunc>(m_lMinH,  m_lMaxH,  m_lMinV,  m_lMaxV, m_lKernWidth, m_lKernHeight, m_corrscore_rejection_threshold, left_image, right_image) );
      } else {
        worker0.reset( new CorrelationWorkThread<ImageView<channel_type>, StandardAbsDifferenceFunc>(-m_lMaxH, -m_lMinH, -m_lMaxV, -m_lMinV, m_lKernWidth, m_lKernHeight, m_corrscore_rejection_threshold, right_image, left_image) );
        worker1.reset( new CorrelationWorkThread<ImageView<channel_type>, StandardAbsDifferenceFunc>(m_lMinH,  m_lMaxH,  m_lMinV,  m_lMaxV, m_lKernWidth, m_lKernHeight, m_corrscore_rejection_threshold, left_image, right_image) );
      }

      try {
        // Launch the threads
        Thread thread0(worker0);
        Thread thread1(worker1);

        // Update the progress info while the threads are running.
        while (!worker0->is_done() || !worker1->is_done()) {
          if (m_verbose) {
            std::cout << boost::format("\t%1$-50s %2$-50s               \r") 
              % worker0->progress_string() 
              % worker1->progress_string();
            fflush(stdout);
          }
          Thread::sleep_ms(100);
        }

        // The threads are destroyed, but the worker objects remain.
        thread0.join();
        thread1.join();
      }
      // FIXME This probably belongs over in Thread.h?
      catch (boost::thread_resource_error &e) {
        vw_throw( LogicErr() << "OptimizedCorrelator: Could not create correlation threads." );
      }

      // Print out some final summary statistics
      if (m_verbose) {
        std::cout << boost::format("\t%1$-50s %2$-50s               \n") 
          % worker0->progress_string() 
          % worker1->progress_string();
        std::cout << boost::format("\t%1$-50s %2$-50s \n") 
          % worker0->correlation_rate_string()
          % worker1->correlation_rate_string();
      }

      // Grab the results
      ImageView<PixelDisparity<float> > resultR2L = worker0->result();
      ImageView<PixelDisparity<float> > resultL2R = worker1->result();

      // Cross check the left and right disparity maps
      cross_corr_consistency_check(resultL2R, resultR2L, m_crossCorrThreshold, m_verbose);

      // Do subpixel correlation
//       NullStereoPreprocessingFilter testfilt;;
//       typename NullStereoPreprocessingFilter::result_type left_subpixel_image = testfilt(channels_to_planes(image0));
//       typename NullStereoPreprocessingFilter::result_type right_subpixel_image = testfilt(channels_to_planes(image1));  
      subpixel_correlation(resultL2R, left_image, right_image, m_lKernWidth, m_lKernHeight, m_useHorizSubpixel, m_useVertSubpixel, m_verbose);

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

      timer.stop();
      double lapse__ = timer.elapsed_seconds();
      if (m_verbose) {
        vw_out(InfoMessage, "stereo") << "\tTotal correlation + subpixel took " << lapse__ << " sec";
        double nTries = 2.0 * (m_lMaxV - m_lMinV + 1) * (m_lMaxH - m_lMinH + 1);
        double rate = nTries * left_image.cols() * left_image.rows() / lapse__ / 1.0e6;
        vw_out(InfoMessage, "stereo") << "\t(" << rate << " M disparities/second)\n";

        double score = (100.0 * matched) / total;
        vw_out(InfoMessage, "stereo") << "\tCorrelation rate: " << score << "\n\n";
      }

      return resultL2R;
    }      
  };


}}   // namespace vw::stereo

#endif /* __OptimizedCorrelator_h__ */
