#include <vw/Image/ImageView.h>
#include <vw/Image/Algorithms.h>
#include <vw/Stereo/OptimizedCorrelator.h>

using namespace vw;

// Boost
#include <boost/format.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/xtime.hpp>

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

// Handy, platform independent utility function for putting boost
// threads to sleep for a specific number of milliseconds.
static void millisecond_sleep(unsigned long milliseconds) {
  
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


namespace vw {
namespace stereo {

  /// The correlation work thread is a functor that does the actual
  /// correlation processing.  This functor can be used to easily spin
  /// off multiple threads of correlation on a multiprocessor machine.
  struct CorrelationWorkThread {
    CorrelationWorkThread() {
      m_parent = NULL;
      m_bit_image_left = NULL; 
      m_bit_image_right = NULL;
      m_image_left = NULL;
      m_image_right = NULL;
      m_done = false;
      m_should_terminate = false;
    }

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
      m_image_width = other.m_image_width;
      m_image_height = other.m_image_height;
      m_kern_height = other.m_kern_height;
      m_kern_width = other.m_kern_width;
      m_bit_image_left = other.m_bit_image_left; 
      m_bit_image_right = other.m_bit_image_right;
      m_image_left = other.m_image_left;
      m_image_right = other.m_image_right;
      m_done = other.m_done;
      m_should_terminate = other.m_should_terminate;
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
      m_image_width = other.m_image_width;
      m_image_height = other.m_image_height;
      m_kern_height = other.m_kern_height;
      m_kern_width = other.m_kern_width;
      m_bit_image_left = other.m_bit_image_left; 
      m_bit_image_right = other.m_bit_image_right;
      m_image_left = other.m_image_left;
      m_image_right = other.m_image_right;
      m_done = other.m_done;
      m_should_terminate = other.m_should_terminate;
      return *this;
    }

    CorrelationWorkThread(OptimizedCorrelator* parent,
                          int id,
                          int min_h, int max_h,
                          int min_v, int max_v,
                          int image_width, int image_height,
                          int kern_width, int kern_height,
                          float* left_image, float* right_image) {
      m_parent = parent;
      m_id = id;
      m_min_h = min_h;
      m_max_h = max_h;
      m_min_v = min_v;
      m_max_v = max_v;
      m_image_width = image_width;
      m_image_height = image_height;
      m_kern_height = kern_height;
      m_kern_width = kern_width;
      m_bit_image_left = NULL; 
      m_bit_image_right = NULL;
      m_image_left = left_image;
      m_image_right = right_image;

      m_done = false;
      m_should_terminate = false;
      m_progress_percent = 0.0;
    }

    CorrelationWorkThread(OptimizedCorrelator* parent, 
                          int id,
                          int min_h, int max_h,
                          int min_v, int max_v,
                          int image_width, int image_height,
                          int kern_width, int kern_height,
                          unsigned char* left_bit_image, unsigned char* right_bit_image) {
      m_parent = parent;
      m_id = id;
      m_min_h = min_h;
      m_max_h = max_h;
      m_min_v = min_v;
      m_max_v = max_v;
      m_image_width = image_width;
      m_image_height = image_height;
      m_kern_height = kern_height;
      m_kern_width = kern_width;
      m_bit_image_left = left_bit_image; 
      m_bit_image_right = right_bit_image;
      m_image_left = NULL;
      m_image_right = NULL;

      m_done = false;
      m_should_terminate = false;
      m_progress_percent = 0.0;
    }
  
    /// Call this operator to start the correlation operation.
    void operator() () {
      set_done(false);
      m_parent->register_worker_thread(m_id, this);
    
      soad* result;
    
      if (m_bit_image_right && m_bit_image_left) {
        result = fast2Dcorr_optimized(m_min_h, m_max_h, m_min_v, m_max_v, m_image_height, m_image_width, m_kern_height, m_kern_width, m_bit_image_left, m_bit_image_right);    
      } else if (m_image_left && m_image_right) {
        result = fast2Dcorr(m_min_h, m_max_h, m_min_v, m_max_v, m_image_height, m_image_width, m_kern_height, m_kern_width, m_image_left, m_image_right);    
      } else {
        vw_throw( vw::LogicErr() << "CorrelationWorkThread: Has not been properly initialized." );
      }
    
      m_result.set_size(m_image_width, m_image_height);
      for (int j = 0; j < m_image_height; j++) {
        for (int i = 0; i < m_image_width; i++) {
          if (result[j*m_image_width+i].hDisp == VW_STEREO_MISSING_PIXEL) {
            m_result(i,j) = PixelDisparity<float>();  // Default constructor creates a missing pixel
          } else {
            m_result(i,j) = PixelDisparity<float>(result[j*m_image_width+i].hDisp,
                                                  result[j*m_image_width+i].vDisp);
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

    bool is_done() {
      boost::mutex::scoped_lock lock(m_mutex);
      return m_done;
    }

    ImageView<PixelDisparity<float> > result() { return m_result; }
  
    void terminate() {
      boost::mutex::scoped_lock lock(m_mutex);
      m_should_terminate = true;
    }
  
    std::string correlation_rate_string() {
      boost::mutex::scoped_lock lock(m_mutex);
      return m_correlation_rate_string;
    }
  
    std::string progress_string() {
      boost::mutex::scoped_lock lock(m_mutex);
      return m_progress_string;
    }
  
    /// Returns a number between 0.0 and 100.0
    double progress_percent() { return m_progress_percent; }
  
  
    void set_correlation_rate_string(std::string str) {
      boost::mutex::scoped_lock lock(m_mutex);
      m_correlation_rate_string = str;
    }
  
    void set_progress_string(std::string str) {
      boost::mutex::scoped_lock lock(m_mutex);
      m_progress_string = str;
    }
  
    void set_done(bool val) {
      boost::mutex::scoped_lock lock(m_mutex);
      m_done = val;
    }
  
  private:
    typedef struct soadStruct {
      double hDisp;	// disparity in x for the best match 
      double vDisp;	// disparity in y for the best match 
      double best;  // soad at the best match 
    } soad;

    soad *fast2Dcorr_optimized(int minDisp,	/* left bound disparity search */
                               int maxDisp,	/* right bound disparity search */
                               int topDisp,	/* top bound disparity search window */
                               int btmDisp,	/* bottom bound disparity search window */ 
                               int height,	/* image height */
                               int width,	/* image width */
                               int vKern,  /* kernel height */
                               int hKern,  /* kernel width */
                               const unsigned char *Rimg,	/* reference image fixed */
                               const unsigned char *Simg	/* searched image sliding */
                               ) {
      
      double before = Time();
      int numCorrTries = (btmDisp - topDisp + 1) * (maxDisp - minDisp + 1);

      uint8 *diff = new uint8[width*height]; // buffer containing results of substraction of L/R images
      VW_ASSERT( diff, NullPtrErr() << "cannot allocate the correlator's difference buffer!" );

      uint16 *cSum = new uint16[width];      // sum per column (kernel hight)
      VW_ASSERT( cSum, NullPtrErr() << "cannot allocate the correlator's column sum buffer!" );

      struct local_result {
        uint16 best, hvdsp;
      };
      local_result *result_buf = new local_result[width * height]; // correlation result buffer
      VW_ASSERT( result_buf, NullPtrErr() << "cannot allocate the correlator's local result buffer!" );
      std::fill((uint16*)result_buf, (uint16*)(result_buf + height * width), uint16(USHRT_MAX));

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
          uint8 *diff_row = diff + yStart*width + xStart;
          const uint8 *rimg_row = Rimg + yStart*width + xStart;
          const uint8 *simg_row = Simg + (yStart+dsy)*width + (xStart+ds);
          for( int j=yEnd-yStart+(vKern-1); j; --j ) {
            uint8 *diff_ptr = diff_row;
            const uint8 *rimg_ptr = rimg_row;
            const uint8 *simg_ptr = simg_row;
            for( int i=xEnd-xStart+(hKern-1); i; --i ) {
              *(diff_ptr++) = (*(rimg_ptr++))^(*(simg_ptr++));
            }
            diff_row += width;
            rimg_row += width;
            simg_row += width;
          }

          // seed the column sum buffer
          for( int i=0; i<width; ++i ) cSum[i] = 0;
          diff_row = diff + yStart*width + xStart;
          uint16 *csum_row = cSum + xStart;
          for( int j=vKern; j; --j ) {
            uint8 *diff_ptr = diff_row;
            uint16 *csum_ptr = csum_row;
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
            uint32 rsum=0;
            uint16 *csum_ptr = cSum + xStart;
            for( int i=hKern; i; --i )
              rsum += *(csum_ptr++);

            // correlate the row
            uint16 *csum_tail = cSum + xStart;
            uint16 *csum_head = csum_tail + hKern;
            for( int i=xEnd-xStart; i; --i, ++result_ptr ) {
              // store the result if better
              if( rsum < result_ptr->best ) {
                result_ptr->best = rsum;
                result_ptr->hvdsp = ds_combined;
              }
              // update the row sum
              rsum += *(csum_head++) - *(csum_tail++);
            }

            // break if this is the last row
            if( j+1 == yEnd ) break;

            // update the column sum
            csum_ptr = cSum + xStart;
            uint8 *diff_tail = diff + xStart + j*width;
            uint8 *diff_head = diff_tail + vKern*width;
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
  
      // convert from the local result buffer to the return format
      for( int nn = 0; nn < height * width; nn++) {
        if( result_buf[nn].best == USHRT_MAX ) {
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
      delete[] result_buf;

      set_progress_string(str(boost::format("V: [%1%,%2%] H: [%3%,%4%] processed successfully.")
                              % topDisp % btmDisp % minDisp % maxDisp));
  
      double duration= Time()-before;
      double mdisp_per_sec= ((double)numCorrTries * width * height)/duration/1e6;
      set_correlation_rate_string(str(boost::format("%1% seconds (%2% M disparities/second)")
                                      % (int)duration % mdisp_per_sec));

      return(result);
    }
  
  
    soad * fast2Dcorr(int minDisp,	/* left bound disparity search */
                      int maxDisp,	/* right bound disparity search */
                      int topDisp,	/* top bound disparity search window */
                      int btmDisp,	/* bottom bound disparity search window */ 
                      int height,	/* image height */
                      int width,	/* image width */
                      int vKern,  /* kernel height */
                      int hKern,  /* kernel width */
                      const float *Rimg,	/* reference image fixed */
                      const float *Simg	/* searched image sliding */
                      ) {
      double before = Time();
      int numCorrTries = (btmDisp - topDisp + 1) * (maxDisp - minDisp + 1);

      float *diff = new float[width*height]; // buffer containing results of substraction of L/R images
      VW_ASSERT( diff, NullPtrErr() << "cannot allocate the correlator's difference buffer!" );

      float *cSum = new float[width];      // sum per column (kernel hight)
      VW_ASSERT( cSum, NullPtrErr() << "cannot allocate the correlator's column sum buffer!" );
  
      soad *result;      // soad buffer to be returned
      result = (soad *)malloc (width*height*sizeof(soad));
      VW_ASSERT( result, NullPtrErr() << "cannot allocate the correlator's SOAD buffer!" );
  
      for( int nn=0; nn<height*width; nn++){
        result[nn].best = VW_STEREO_MISSING_PIXEL;
        result[nn].hDisp = VW_STEREO_MISSING_PIXEL;
        result[nn].vDisp = VW_STEREO_MISSING_PIXEL;
      }

      // for each allowed disparity...
      for( int dsy=topDisp; dsy<=btmDisp; ++dsy ) {
        for( int ds=minDisp; ds<=maxDisp; ++ds ) {
          set_progress_string(str(boost::format("V: [%1%,%2%] H: [%3%,%4%] processing %5% %6%")
                                  % topDisp % btmDisp % minDisp % maxDisp % dsy % ds));
      
          // compute the region of correlation
          int yStart = (dsy<0) ? (-dsy) : 0;
          int xStart = (ds<0) ? (-ds) : 0;
          int yEnd = yStart + height - abs(dsy) - (vKern-1);
          int xEnd = xStart + width - abs(ds) - (hKern-1);

          // compute the difference buffer
          float *diff_row = diff + yStart*width + xStart;
          const float *rimg_row = Rimg + yStart*width + xStart;
          const float *simg_row = Simg + (yStart+dsy)*width + (xStart+ds);
          for( int j=yEnd-yStart+(vKern-1); j; --j ) {
            float *diff_ptr = diff_row;
            const float *rimg_ptr = rimg_row;
            const float *simg_ptr = simg_row;
            for( int i=xEnd-xStart+(hKern-1); i; --i ) {
              *(diff_ptr++) = fabs((float)*(rimg_ptr++)-(float)*(simg_ptr++));
            }
            diff_row += width;
            rimg_row += width;
            simg_row += width;
          }

          // seed the column sum buffer
          for( int i=0; i<width; ++i ) cSum[i] = 0;
          diff_row = diff + yStart*width + xStart;
          float *csum_row = cSum + xStart;
          for( int j=vKern; j; --j ) {
            float *diff_ptr = diff_row;
            float *csum_ptr = csum_row;
            for( int i=xEnd-xStart+(hKern-1); i; --i ) {
              *(csum_ptr++) += *(diff_ptr++);
            }
            diff_row += width;
          }

          // perform the correlation
          soad *result_row = result+(xStart+hKern/2)+(yStart+vKern/2)*width;
          for( int j=yStart; ; ++j, result_row+=width ) {
            soad *result_ptr = result_row;

            // seed the row sum (i.e. soad)
            float rsum=0;
            float *csum_ptr = cSum + xStart;
            for( int i=hKern; i; --i )
              rsum += *(csum_ptr++);

            // correlate the row
            float *csum_tail = cSum + xStart;
            float *csum_head = csum_tail + hKern;
            for( int i=xEnd-xStart; i; --i, ++result_ptr ) {
              // store the result if better
              if( result_ptr->best==VW_STEREO_MISSING_PIXEL || rsum < result_ptr->best ) {
                result_ptr->best = rsum;
                result_ptr->hDisp = ds;
                result_ptr->vDisp = dsy;
              }
              // update the row sum
              rsum += *(csum_head++) - *(csum_tail++);
            }

            // break if this is the last row
            if( j+1 == yEnd ) break;

            // update the column sum
            csum_ptr = cSum + xStart;
            float *diff_tail = diff + xStart + j*width;
            float *diff_head = diff_tail + vKern*width;
            for( int i=xEnd-xStart+(hKern-1); i; --i ) {
              *(csum_ptr++) += *(diff_head++) - *(diff_tail++);
            }
          }
        }
      }

      delete[] diff;
      delete[] cSum;

      set_progress_string(str(boost::format("V: [%1%,%2%] H: [%3%,%4%] processed successfully.")
                              % topDisp % btmDisp % minDisp % maxDisp));
	
      double duration= Time()-before;
      double mdisp_per_sec= ((double)numCorrTries * width * height)/duration/1e6;
      set_correlation_rate_string(str(boost::format("%1% seconds (%2% M disparities/second)")
                                      % (int)duration % mdisp_per_sec));

      return(result);
    }


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
  m_useHorizSubpixel = useSubpixelH;
  m_useVertSubpixel = useSubpixelV;
  m_children = std::vector<CorrelationWorkThread*>(2);
  m_children[0] = m_children[1] = NULL;
}

/********************************************************************
 *                     Private Methods                              *
 *******************************************************************/

void OptimizedCorrelator::register_worker_thread(int id, CorrelationWorkThread* child) {
  boost::mutex::scoped_lock lock(m_mutex);
  m_children[id] = child;
}

CorrelationWorkThread* OptimizedCorrelator::get_worker_thread(int id) {
  boost::mutex::scoped_lock lock(m_mutex);
  return m_children[id];
}


template <class ChannelT>
ImageView<PixelDisparity<float> > OptimizedCorrelator::correlation_2D(ImageView<ChannelT> &left_image,
                                                                      ImageView<ChannelT> &right_image, 
                                                                      bool bit_image) {

  VW_ASSERT(left_image.cols() == right_image.cols() &&
            left_image.rows() == right_image.rows(),
            ArgumentErr() << "subpixel_correlation: input image dimensions do not agree.\n");


  int width = left_image.cols();
  int height = left_image.rows();
      
  //Run the correlator and record how long it takes to run.
  double begin__ = Time();
  
  try {
    boost::thread_group threads;
    if (bit_image == true) {
      if (m_verbose) 
        std::cout << "\tUsing bit image optimized correlator\n";
      uint8* bit_img0 = (uint8*)&(left_image(0,0));
      uint8* bit_img1 = (uint8*)&(right_image(0,0));
      threads.create_thread(CorrelationWorkThread(this, 0, -m_lMaxH, -m_lMinH, -m_lMaxV, -m_lMinV, width, height, m_lKernWidth, m_lKernHeight, bit_img1, bit_img0));
      threads.create_thread(CorrelationWorkThread(this, 1,  m_lMinH,  m_lMaxH,  m_lMinV,  m_lMaxV, width, height, m_lKernWidth, m_lKernHeight, bit_img0, bit_img1));
    } else {
      if (m_verbose)
        std::cout << "\tUsing standard correlator\n";
      float* pBuffer0 = (float*)&(left_image(0,0));
      float* pBuffer1 = (float*)&(right_image(0,0));
      threads.create_thread(CorrelationWorkThread(this, 0, -m_lMaxH, -m_lMinH, -m_lMaxV, -m_lMinV, width, height, m_lKernWidth, m_lKernHeight, pBuffer1, pBuffer0));
      threads.create_thread(CorrelationWorkThread(this, 1,  m_lMinH,  m_lMaxH,  m_lMinV,  m_lMaxV, width, height, m_lKernWidth, m_lKernHeight, pBuffer0, pBuffer1));
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
    subpixel_correlation(resultL2R, left_image, right_image, m_lKernWidth, m_lKernHeight, m_useHorizSubpixel, m_useVertSubpixel);

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
      cout << "\tTotal correlation + subpixel took " << lapse__ << " sec";
      double nTries = 2.0 * (m_lMaxV - m_lMinV + 1) * (m_lMaxH - m_lMinH + 1);
      double rate = nTries * width * height / lapse__ / 1.0e6;
      cout << "\t(" << rate << " M disparities/second)" << endl;

      double score = (100.0 * matched) / total;
      printf("\tCorrelation rate: %6.4f\n\n", score);
    }


    return resultL2R;
    
  } catch (boost::thread_resource_error &e) {
    vw_throw( LogicErr() << "OptimizedCorrelator: Could not create correlation threads." );
  }
  
}

// Explicit Intantiation
namespace vw {
namespace stereo {

template
ImageView<PixelDisparity<float> > OptimizedCorrelator::correlation_2D<float>( ImageView<float> &left_image,
                                                                             ImageView<float> &right_image, 
                                                                             bool bit_image);

template
ImageView<PixelDisparity<float> > OptimizedCorrelator::correlation_2D<uint8>( ImageView<uint8> &left_image,
                                                                             ImageView<uint8> &right_image, 
                                                                             bool bit_image);

}}

