#include <vw/Image/ImageView.h>
#include <vw/Image/Algorithms.h>
#include <vw/Stereo/OptimizedCorrelator.h>

using namespace vw;

// Boost
#include <boost/format.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/xtime.hpp>
#include <boost/utility/enable_if.hpp>

using namespace std;
using namespace vw;
using namespace vw::stereo;

// Some useful default values and constants 
#define DEFAULT_KERN_WIDTH 29
#define DEFAULT_KERN_HEIGHT 21
#define DEFAULT_MIN_H -50
#define DEFAULT_MAX_H 50
#define DEFAULT_MIN_V -1
#define DEFAULT_MAX_V 1
#define DEFAULT_CROSSCORR_THRESHOLD 2
#define DEFAULT_USE_SUBPIXEL_H true
#define DEFAULT_USE_SUBPIXEL_V false

// Switches to turn on/off debugging 
#define VW_DEBUG_CORRELATOR 0


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
template <> struct CorrelatorAccumulatorType<vw::uint8>   { typedef vw::uint16   type; };
template <> struct CorrelatorAccumulatorType<vw::int8>    { typedef vw::uint16   type; };
template <> struct CorrelatorAccumulatorType<vw::uint16>  { typedef vw::uint32   type; };
template <> struct CorrelatorAccumulatorType<vw::int16>   { typedef vw::uint32   type; };
template <> struct CorrelatorAccumulatorType<vw::uint32>  { typedef vw::uint64   type; };
template <> struct CorrelatorAccumulatorType<vw::int32>   { typedef vw::uint64   type; };
template <> struct CorrelatorAccumulatorType<vw::float32> { typedef vw::float32 type; };
template <> struct CorrelatorAccumulatorType<vw::float64> { typedef vw::float64 type; };

namespace vw {
namespace stereo {

  // Handy, platform independent utility function for putting boost
  // threads to sleep for a specific number of milliseconds.
  void millisecond_sleep(unsigned long milliseconds) {
    
    static long const nanoseconds_per_second = 1000L*1000L*1000L;
    boost::xtime xt;
    boost::xtime_get(&xt, boost::TIME_UTC);
    xt.nsec+=1000*1000*milliseconds;
    while (xt.nsec > nanoseconds_per_second) {
      xt.nsec -= nanoseconds_per_second;
      xt.sec++;
    }
    
    boost::thread::sleep(xt);
  }
  
  class CorrelationWorkThreadBase {
  protected:
    OptimizedCorrelator* m_parent;
    int m_id;
    ImageView<PixelDisparity<float> > m_result;

    std::string m_correlation_rate_string;
    double m_progress_percent;
    std::string m_progress_string;
    bool m_done;
    bool m_should_terminate;
    mutable boost::mutex m_mutex;
    char temp_progress[100];
  public:
    virtual bool is_done() {
      boost::mutex::scoped_lock lock(m_mutex);
      return m_done;
    }

    virtual ImageView<PixelDisparity<float> > result() { return m_result; }
  
    virtual void terminate() {
      boost::mutex::scoped_lock lock(m_mutex);
      m_should_terminate = true;
    }
  
    virtual std::string correlation_rate_string() {
      boost::mutex::scoped_lock lock(m_mutex);
      return m_correlation_rate_string;
    }
  
    virtual std::string progress_string() {
      boost::mutex::scoped_lock lock(m_mutex);
      return m_progress_string;
    }
  
    /// Returns a number between 0.0 and 100.0
    virtual double progress_percent() { return m_progress_percent; }
  
  
    virtual void set_correlation_rate_string(std::string str) {
      boost::mutex::scoped_lock lock(m_mutex);
      m_correlation_rate_string = str;
    }
  
    virtual void set_progress_string(std::string str) {
      boost::mutex::scoped_lock lock(m_mutex);
      m_progress_string = str;
    }
  
    virtual void set_done(bool val) {
      boost::mutex::scoped_lock lock(m_mutex);
      m_done = val;
    }
  };


  /// The correlation work thread is a functor that does the actual
  /// correlation processing.  This functor can be used to easily spin
  /// off multiple threads of correlation on a multiprocessor machine.
  template <class ViewT, class AbsDifferenceT>
  class CorrelationWorkThread : CorrelationWorkThreadBase {

    // Handy typedefs
    typedef typename ViewT::pixel_type pixel_type;
    typedef typename PixelChannelType<pixel_type>::type channel_type;
    
    // This is used as an intermediate format for the correlator
    typedef struct soadStruct {
      double hDisp;	// disparity in x for the best match 
      double vDisp;	// disparity in y for the best match 
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

      double before = Time();
      int numCorrTries = (btmDisp - topDisp + 1) * (maxDisp - minDisp + 1);

      ChannelT *diff = new ChannelT[width*height]; // buffer containing results of substraction of L/R images
      VW_ASSERT( diff, NullPtrErr() << "cannot allocate the correlator's difference buffer!" );

      accumulator_type *cSum = new accumulator_type[width];      // sum per column (kernel hight)
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
  
      double duration= Time()-before;
      double mdisp_per_sec= ((double)numCorrTries * width * height)/duration/1e6;
      set_correlation_rate_string(str(boost::format("%1% seconds (%2% M disparities/second)")
                                      % (int)duration % mdisp_per_sec));

      return(result);
    }


  
    // Hide this constructor -- it should never be called.
    CorrelationWorkThread() {}

  public:

    // The copy constructor must be coded explicitly be cause the mutex
    // object is non-copyable.
    CorrelationWorkThread(const CorrelationWorkThread& other) {
      boost::mutex::scoped_lock lock(other.m_mutex);
      m_parent = other.m_parent;
      m_id = other.m_id;
      m_result = other.m_result;
      m_min_h = other.m_min_h;
      m_max_h = other.m_max_h;
      m_min_v = other.m_min_v;
      m_max_v = other.m_max_v;
      m_kern_height = other.m_kern_height;
      m_kern_width = other.m_kern_width;
      m_left_image = other.m_left_image;
      m_right_image = other.m_right_image;
      m_done = other.m_done;
      m_should_terminate = other.m_should_terminate;
      m_corrscore_rejection_threshold = other.m_corrscore_rejection_threshold;
    }

    // For assignment we need to synchronize both objects!
    const CorrelationWorkThread& operator=(const CorrelationWorkThread& other) {
      if (this == &other)
        return *this;
      boost::mutex::scoped_lock lock1(&m_mutex < &other.m_mutex ? m_mutex : other.m_mutex);
      boost::mutex::scoped_lock lock2(&m_mutex > &other.m_mutex ? m_mutex : other.m_mutex);
      m_parent = other.m_parent;
      m_id = other.m_id;
      m_result = other.m_result;
      m_min_h = other.m_min_h;
      m_max_h = other.m_max_h;
      m_min_v = other.m_min_v;
      m_max_v = other.m_max_v;
      m_kern_height = other.m_kern_height;
      m_kern_width = other.m_kern_width;
      m_left_image = other.m_left_image;
      m_right_image = other.m_right_image;
      m_done = other.m_done;
      m_should_terminate = other.m_should_terminate;
      m_corrscore_rejection_threshold = other.m_corrscore_rejection_threshold;
      return *this;
    }

    CorrelationWorkThread(OptimizedCorrelator* parent,
                          int id,
                          int min_h, int max_h,
                          int min_v, int max_v,
                          int kern_width, int kern_height,
                          float corrscore_rejection_threshold,
                          ViewT const& left_image, 
                          ViewT const& right_image) {
      m_parent = parent;
      m_id = id;
      m_min_h = min_h;
      m_max_h = max_h;
      m_min_v = min_v;
      m_max_v = max_v;
      m_kern_height = kern_height;
      m_kern_width = kern_width;
      m_left_image = left_image;
      m_right_image = right_image;
      m_corrscore_rejection_threshold = corrscore_rejection_threshold;

      m_done = false;
      m_should_terminate = false;
      m_progress_percent = 0.0;
    }
  
    /// Call this operator to start the correlation operation.
    void operator() () {
      set_done(false);
      m_parent->register_worker_thread(m_id, this);
    
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
    
      // After the thread has finished doing the work, it should hang
      // around so that the parent has a chance to request any final state
      // from the child.  When the parent is finished, it can dismiss the
      // child with the terminate() method.
      while (!m_should_terminate) { millisecond_sleep(100); }
    }

  };
}} // namespace vw::stereo


/********************************************************************
 *                     Constructors                                 *
 *******************************************************************/
OptimizedCorrelator::OptimizedCorrelator()
{  

  m_lKernWidth = DEFAULT_KERN_WIDTH;
  m_lKernHeight = DEFAULT_KERN_HEIGHT;

  m_lMinH = DEFAULT_MIN_H;
  m_lMaxH = DEFAULT_MAX_H;
  m_lMinV = DEFAULT_MIN_V;
  m_lMaxV = DEFAULT_MAX_V;

  m_verbose = true;

  m_crossCorrThreshold = DEFAULT_CROSSCORR_THRESHOLD;
  m_corrscore_rejection_threshold = 1.0; // 1.0 turns of this type of blunder rejection 
  m_useHorizSubpixel = DEFAULT_USE_SUBPIXEL_H;
  m_useVertSubpixel = DEFAULT_USE_SUBPIXEL_V;
}

OptimizedCorrelator::OptimizedCorrelator(int minH,	/* left bound disparity search window*/
                                         int maxH,	/* right bound disparity search window*/
                                         int minV,	/* bottom bound disparity search window */ 
                                         int maxV,	/* top bound disparity search window */
                                         int kernWidth,	/* size of the kernel */
                                         int kernHeight,       
                                         int verbose,
                                         double crosscorrThreshold,
                                         float corrscore_rejection_threshold,
                                         int useSubpixelH,
                                         int useSubpixelV)
{
  m_lKernWidth = kernWidth;
  m_lKernHeight = kernHeight;
  m_lMinH = minH;
  m_lMaxH = maxH;
  m_lMinV = minV;
  m_lMaxV = maxV;  
  m_verbose = verbose;

  m_crossCorrThreshold = crosscorrThreshold;
  m_corrscore_rejection_threshold = corrscore_rejection_threshold;
  m_useHorizSubpixel = useSubpixelH;
  m_useVertSubpixel = useSubpixelV;
  m_children = std::vector<CorrelationWorkThreadBase*>(2);
  m_children[0] = m_children[1] = NULL;
}

template <class ChannelT>
ImageView<PixelDisparity<float> > OptimizedCorrelator::do_correlation(ImageView<ChannelT> &left_image, ImageView<ChannelT> &right_image, bool bit_image) {

  //Run the correlator and record how long it takes to run.
  double begin__ = Time();
      
  try {
    boost::thread_group threads;
    if (bit_image) {
      threads.create_thread(CorrelationWorkThread<ImageView<ChannelT>, BitAbsDifferenceFunc>(this, 0, -m_lMaxH, -m_lMinH, -m_lMaxV, -m_lMinV, m_lKernWidth, m_lKernHeight, m_corrscore_rejection_threshold, right_image, left_image));
      threads.create_thread(CorrelationWorkThread<ImageView<ChannelT>, BitAbsDifferenceFunc>(this, 1,  m_lMinH,  m_lMaxH,  m_lMinV,  m_lMaxV, m_lKernWidth, m_lKernHeight, m_corrscore_rejection_threshold, left_image, right_image));
    } else {
      threads.create_thread(CorrelationWorkThread<ImageView<ChannelT>, StandardAbsDifferenceFunc>(this, 0, -m_lMaxH, -m_lMinH, -m_lMaxV, -m_lMinV, m_lKernWidth, m_lKernHeight, m_corrscore_rejection_threshold, right_image, left_image));
      threads.create_thread(CorrelationWorkThread<ImageView<ChannelT>, StandardAbsDifferenceFunc>(this, 1,  m_lMinH,  m_lMaxH,  m_lMinV,  m_lMaxV, m_lKernWidth, m_lKernHeight, m_corrscore_rejection_threshold, left_image, right_image));
    }

    // Wait until the child threads register with their parent, then begin to query them 
    // for information about their current state.
    while (!get_worker_thread(0) || !get_worker_thread(1));
    while (!get_worker_thread(0)->is_done() || !get_worker_thread(1)->is_done()) {
      if (m_verbose) {
        std::cout << boost::format("\t%1$-50s %2$-50s               \r") 
          % get_worker_thread(0)->progress_string() 
          % get_worker_thread(1)->progress_string();
        fflush(stdout);
      }
      millisecond_sleep(100);
    }

    // Once the threads have complete their task, they hang around so
    // that we can print out some summary statistics.  We then
    // terminate them and join them back into the main thread.
    if (m_verbose) {
      std::cout << boost::format("\t%1$-50s %2$-50s               \n") 
        % get_worker_thread(0)->progress_string() 
        % get_worker_thread(1)->progress_string();
      std::cout << boost::format("\t%1$-50s %2$-50s \n") 
        % get_worker_thread(0)->correlation_rate_string()
        % get_worker_thread(1)->correlation_rate_string();
    }

    // Ask the worker threads for the actual results of the disparity correlation
    ImageView<PixelDisparity<float> > resultR2L = get_worker_thread(0)->result();
    ImageView<PixelDisparity<float> > resultL2R = get_worker_thread(1)->result();

    // We now have all the information we need.  Shut down the worker
    // threads and wait for them to finish terminating.
    get_worker_thread(0)->terminate();
    get_worker_thread(1)->terminate();
    threads.join_all(); 
  
    // Cross check the left and right disparity maps
    cross_corr_consistency_check(resultL2R, resultR2L,m_crossCorrThreshold, m_verbose);

    // Do subpixel correlation
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

    double lapse__ = Time() - begin__;
    if (m_verbose) {
      vw_out(InfoMessage) << "\tTotal correlation + subpixel took " << lapse__ << " sec";
      double nTries = 2.0 * (m_lMaxV - m_lMinV + 1) * (m_lMaxH - m_lMinH + 1);
      double rate = nTries * left_image.cols() * left_image.rows() / lapse__ / 1.0e6;
      vw_out(InfoMessage) << "\t(" << rate << " M disparities/second)\n";

      double score = (100.0 * matched) / total;
      vw_out(InfoMessage) << "\tCorrelation rate: " << score << "\n\n";
    }


    return resultL2R;
    
  } catch (boost::thread_resource_error &e) {
    vw_throw( LogicErr() << "OptimizedCorrelator: Could not create correlation threads." );
  }

}

// Explicit template initialization
namespace vw { namespace stereo {
    template ImageView<PixelDisparity<float> > OptimizedCorrelator::do_correlation<uint8>(ImageView<uint8> &left_image, ImageView<uint8> &right_image, bool bit_image);
    template ImageView<PixelDisparity<float> > OptimizedCorrelator::do_correlation<float>(ImageView<float> &left_image, ImageView<float> &right_image, bool bit_image);
}}
