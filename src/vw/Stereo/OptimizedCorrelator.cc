
#include <vw/Image/ImageView.h>
#include <vw/Image/Algorithms.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>
#include <vw/Stereo/OptimizedCorrelator.h>

using namespace vw;

#include <sys/time.h>

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

/* For Eric's correlator */
#define DEFAULT_DIFF 30000

/* Switches to turn on/off debugging */
#define VW_DEBUG_CORRELATOR 0

// Handy, platform independent utility function for putting boost
// threads to sleep for a specific number of milliseconds.
static void 
millisecond_sleep(unsigned long milliseconds) {
  
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

// The following has been retained for backwards compatability, 
// but is deprecated and should be eliminated.  In the modern 
// era, a class like cl_bit_xor has no reason to be.  It 
// should surely be replaced with an unnamed boost functor 
// of the form (_1^_2).

/* Adapted somewhat from Claraty's share/claraty_functors.h */

template <class arg1, class arg2, class res>
struct cl_binary_function {
  typedef arg1 first_argument_type;
  typedef arg2 second_argument_type;
  typedef res result_type;
};

#define __cl_bin(A1, A2, R) cl_binary_function<A1, A2, R>

template <class T>
struct cl_bit_xor : public __cl_bin(T, T, T) {
  T operator()(const T& x, const T& y) { return x ^ y; }
};

/// The correlation work thread is a functor that does the actual
/// correlation processing.  This functor can be used to easily spin
/// off multiple threads of correlation on a multiprocessor machine.

// The copy constructor must be coded explicitly be cause the mutex
// object is non-copyable.
CorrelationWorkThread::CorrelationWorkThread(const CorrelationWorkThread& other) {
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
const CorrelationWorkThread& CorrelationWorkThread::operator=(const CorrelationWorkThread& other) {
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

CorrelationWorkThread::CorrelationWorkThread() {
    m_parent = NULL;
    m_bit_image_left = NULL; 
    m_bit_image_right = NULL;
    m_image_left = NULL;
    m_image_right = NULL;
    m_done = false;
    m_should_terminate = false;
  }

CorrelationWorkThread::CorrelationWorkThread(SubpixelCorrelator* parent,
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

CorrelationWorkThread::CorrelationWorkThread(SubpixelCorrelator* parent, 
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

void CorrelationWorkThread::operator() () {
  set_done(false);
  m_parent->register_worker_thread(m_id, this);
  
  SubpixelCorrelator::soad* result;

  if (m_bit_image_right && m_bit_image_left) {
    result = fast2Dcorr_optimized(m_min_h, m_max_h, m_min_v, m_max_v, m_image_height, m_image_width, m_kern_height, m_kern_width, m_bit_image_left, m_bit_image_right);    
  } else if (m_image_left && m_image_right) {
    result = fast2Dcorr(m_min_h, m_max_h, m_min_v, m_max_v, m_image_height, m_image_width, m_kern_height, m_kern_width, m_image_left, m_image_right);    
  } else {
    throw vw::LogicErr() << "CorrelationWorkThread: Has not been properly initialized.";
  }

  m_result.set_size(m_image_width, m_image_height);
  for (int i = 0; i < m_image_width; i++) {
    for (int j = 0; j < m_image_height; j++) {
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
    
/// Randy's 2D correlator capitolizes on bit images (i.e. slog) that can be 
/// efficiently differenced using an xor operation.
SubpixelCorrelator::soad *CorrelationWorkThread::fast2Dcorr_optimized(
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
  ) {

    int nn,xx,yy, index;	/* loops counters */
    int yTop, yBtm;
    int xStart, xEnd;	    /* end of loop value: L/R pixel difference */
    int sumStart, sumEnd;	/* end of summing loop value */
    int ds, dsy;		      /* disparity shift between the two images */
    int *cSum;		        /* sum per column (kernel hight) */
    SubpixelCorrelator::soad *result;		/* soad buffer to be returned */
    int ii;		            /* loop counter */
    unsigned short sum;	/* sum over the kernel */
    unsigned short *diff;	/* buffer containing results of subtraction of L/R images */

    unsigned short *result_best;
    unsigned short *result_hvdisp;
    unsigned short ds_combined;
    const unsigned short MAX_BEST = USHRT_MAX;
    const unsigned short MISSING_HVDISP = USHRT_MAX;

    double before = Time();
    int numCorrTries = (btmDisp - topDisp + 1) * (maxDisp - minDisp + 1);

    /* initialization */
    cSum          = new int[width];
    result        = new SubpixelCorrelator::soad[width * height];
    diff          = new unsigned short[width * height];
    result_best   = new unsigned short[width * height];
    result_hvdisp = new unsigned short[width * height];
    
    std::fill(result_best, result_best + height * width, MAX_BEST);
    std::fill(result_hvdisp, result_hvdisp + height * width, MISSING_HVDISP);
    std::fill(cSum, cSum + width, 0);

    /* subtract the two images with the desired disparity shift */
    for(dsy = topDisp; dsy <= btmDisp; dsy++) {
      for(ds = minDisp; ds <= maxDisp; ds++) {
        set_progress_string(str(boost::format("H: [%1%,%2%] V: [%3%,%4%] processing %5% %6%")
                                % topDisp % btmDisp % minDisp % maxDisp % dsy % ds));
        
        ds_combined= (ds-minDisp) + (dsy - topDisp) * (maxDisp - minDisp + 1);
        sumStart = max(-(dsy * width + ds), 0);
        sumEnd = min(width - (dsy * width + ds), width);
        xStart = max(-ds, 0);
        xEnd = min(width - ds, width);
                        
        // initialize the cSum and diff buffers
        std::fill(cSum, cSum + width, 0);
        std::fill(diff, diff + sumStart, 1);

        // diff the reference and sliding buffers. This is rather ugly, but
        // allows the optimizer to go crazy, particularly with loop
        // unrolling. What we're basically doing is xoring the shifted pixels,
        // and filling the ones that don't overlap with 1.
        std::transform(Rimg + sumStart, Rimg + (height - 1) * width + sumEnd,
                       Simg + sumStart + dsy * width + ds, diff + sumStart,
                       cl_bit_xor<unsigned char>());
        
        std::fill(diff + (height - 1) * width + sumEnd, diff + width * height, 1);
      
        /* ramp up - fill cSum buffer at top of image & update first row */
        for(xx=xStart; xx<xEnd; xx++)
          /* update the cSum buffer */
          for(yy=0, index = xx; yy<kernel_height; yy++, index+=width)
            cSum[xx] += diff[index];
      
        /* 
         * Compute soads for the current row (kernel_width/2) 
         * getSoad2D(next, cSum, xStart, xEnd, kernel_width); inlined 
         * do the first point.
         */
        sum=0;
        for (ii=xStart; ii<xStart+kernel_width; ii++)
          sum += cSum[ii];
      
        if(result_best[(kernel_width/2)*width+(xStart+kernel_width/2)] > sum) {
          result_best[(kernel_width/2)*width+(xStart+kernel_width/2)] = sum;
          result_hvdisp[(kernel_width/2)*width+(xStart+kernel_width/2)] = ds_combined;
        }
      
        /* do the rest of the line */
        for(ii=xStart+kernel_width/2+1; ii<xEnd-kernel_width/2; ii++) {
          sum += (cSum[ii+kernel_width/2]-cSum[ii-kernel_width/2-1]);
          if(result_best[(kernel_width/2)*width+ii] > sum) {
            result_best[(kernel_width/2)*width+ii] = sum;
            result_hvdisp[(kernel_width/2)*width+ii] = ds_combined;
          }
        }
      
        /* do the rest of the image */
        for(yy=kernel_height/2+1, yTop = 0, yBtm = kernel_height*width; yy<height-kernel_height/2; yy++, yBtm+=width, yTop+=width) {
          /* update cSum buffer - remove top kernel row and add bottom one */
          for(xx=xStart; xx<xEnd; xx++) {
            cSum[xx] += (diff[yBtm+xx]-diff[yTop+xx]);
          }
	
          /* compute soads for the current row (yy) */
          /*getSoad2D(next, cSum, xStart, xEnd, kernel_width); inlined */
          /* do the first point */
          sum=0;
        
          for (ii=xStart; ii<xStart+kernel_width; ii++) {
            sum += cSum[ii];
          }
        
          if(result_best[yy*width+(xStart+kernel_width/2)] > sum) {
            result_best[yy*width+(xStart+kernel_width/2)] = sum;
            result_hvdisp[yy*width+(xStart+kernel_width/2)] = ds_combined;
          }
        
          /* do the rest of the line */
          for(ii=xStart+kernel_width/2+1; ii<xEnd-kernel_width/2; ii++) {
            sum += (cSum[ii+kernel_width/2]-cSum[ii-kernel_width/2-1]);
            if(result_best[yy*width+ii] > sum) {
              result_best[yy*width+ii] = sum;
              result_hvdisp[yy*width+ii] = ds_combined;
            }
          }
        }
      }
    }
    delete [] diff;
    delete [] cSum;
	
    set_progress_string(str(boost::format("H: [%1%,%2%] V: [%3%,%4%] processed succefully.")
                            % topDisp % btmDisp % minDisp % maxDisp));
	
    double duration = Time() - before;
    double mdisp_per_sec = (double)numCorrTries * width * height / duration / 1.0e6;
    set_correlation_rate_string(str(boost::format("%1% seconds (%2% M disparities/second)")
                                    % (int)duration % mdisp_per_sec));
    
    for(nn = 0; nn < height * width; nn++) {
      if (result_best[nn] == MAX_BEST) {
        result[nn].best = VW_STEREO_MISSING_PIXEL;
        result[nn].hDisp = VW_STEREO_MISSING_PIXEL;
        result[nn].vDisp = VW_STEREO_MISSING_PIXEL;
      } else {
        result[nn].best = result_best[nn];
      
        int hvDisp= result_hvdisp[nn];
        int hDisp= hvDisp % (maxDisp-minDisp+1) + minDisp;
        int vDisp= hvDisp / (maxDisp-minDisp+1) + topDisp;
        result[nn].hDisp = hDisp;
        result[nn].vDisp = vDisp;
      }
    }
  
    delete [] result_hvdisp;
    delete [] result_best;
  
    return result;
}

/// Eric's 2D correlator does not assume anything about the type of
/// the image, but it is also not optimized, and is about 2 or 3
/// times slower than Randy's correlator.
SubpixelCorrelator::soad *CorrelationWorkThread::fast2Dcorr(int minDisp,	/* left bound disparity search */
                                                            int maxDisp,	/* right bound disparity search */
                                                            int topDisp,	/* top bound disparity search window */
                                                            int btmDisp,	/* bottom bound disparity search window */ 
                                                            int height,	/* image height */
                                                            int width,	/* image width */
                                                            int vKern,  /* kernel height */
                                                            int hKern,  /* kernel width */
                                                            float *Rimg,	/* reference image fixed */
                                                            float *Simg	/* searched image sliding */
                                                            ) {

    int verbose = 1;
    int nn,xx,yy, index;	/* loops counters */
    int yTop, yBtm;
    int xStart, xEnd;	/* end of loop value: L/R pixel difference */
    int sumStart, sumEnd;	/* end of summing loop value */
    int ds, dsy;	/* disparity shift between the two images */
    float *diff;	/* buffer containing results of substraction of L/R images */
    float *cSum;	/* sum per column (kernel hight) */
    float *next;	/* work in buffer (soad results of current pass (ds) ) */
    SubpixelCorrelator::soad *result;	/* soad buffer to be returned */
    int ii;	/* loop counter */
    float sum;	/* sum over the kernel */

    double before = Time();
    int numCorrTries = (btmDisp - topDisp + 1) * (maxDisp - minDisp + 1);

    /* initialization */
    if ( !( diff = (float *)malloc (width*height*sizeof(float)))) {
      printf("cannot allocate fast1Dcorr's difference buffer\n");
      exit(1);
    }
    if ( !( cSum = (float *)malloc (width*sizeof(float)))) {
      printf("cannot allocate fast1Dcorr's cSum buffer\n");
      exit(1);
    }
    if ( !( next = (float *)malloc (width*sizeof(float)))) {
      printf("cannot allocate fast1Dcorr's next buffer\n");
      exit(1);
    }
    if ( !( result = (SubpixelCorrelator::soad *)malloc (width*height*sizeof(SubpixelCorrelator::soad)))) {
      printf("cannot allocate fast1Dcorr's SOAD buffer\n");
      exit(1);
    }
    
    for(nn=0; nn<height*width; nn++){
      result[nn].best = VW_STEREO_MISSING_PIXEL;
      result[nn].hDisp = VW_STEREO_MISSING_PIXEL;
      result[nn].vDisp = VW_STEREO_MISSING_PIXEL;
    }
    for(nn=0; nn<width; nn++){
      next[nn] = VW_STEREO_MISSING_PIXEL;
      cSum[nn] = 0;
    }
    
    // subtract the two images with the desired disparity shift 
    for(dsy=topDisp; dsy<= btmDisp; dsy++) {
      for(ds=minDisp; ds<=maxDisp; ds++) {
        set_progress_string(str(boost::format("H: [%1%,%2%] V: [%3%,%4%] processing %5% %6%")
                                % topDisp % btmDisp % minDisp % maxDisp % dsy % ds));
        progress_string();

        sumStart = (dsy*width+ds<0 ? -(dsy*width+ds) : 0);
        sumEnd = (dsy*width+ds<0 ? width : width-(dsy*width+ds));
        xStart = (ds<0 ? -ds : 0);
        xEnd = (ds<0 ? width : width-ds);
#if VW_DEBUG_CORRELATOR
        assert(xStart>=0);
        assert(xStart<=width);
        assert(xEnd>=0);
        assert(xEnd<=width);
#endif
        /* initialize the cSum buffer */
        for(nn=0; nn<width; nn++)
          cSum[nn] = 0;
        
        /* initialize the diff buffer */
        for(nn=0; nn<width*height; nn++)
          diff[nn] = DEFAULT_DIFF;
        
        /* diff the reference and sliding buffers */
        for(nn=sumStart; nn<(height-1)*width+sumEnd; nn++){
#if VW_DEBUG_CORRELATOR
          assert(nn>=0);
          assert(nn<width*height);
          assert(nn+dsy*width+ds>=0);
          assert(nn+dsy*width+ds<width*height);
#endif
          diff[nn] = fabs((float)Rimg[nn]-(float)Simg[nn+dsy*width+ds]);
        }
        
        /* ramp up - fill cSum buffer at top of image & update first row */
        for(xx=xStart; xx<xEnd; xx++)
          /* update the cSum buffer */
          for(yy=0, index = xx; yy<vKern; yy++, index+=width)
            cSum[xx] += diff[index];
        /* compute soads for the current row (hKern/2) */
        /*getSoad2D(next, cSum, xStart, xEnd, hKern); inlined */
        /* do the first point */
        sum=0;
        for (ii=xStart; ii<xStart+hKern; ii++)
          sum += cSum[ii];
        next[xStart+hKern/2] = sum;
        /* do the rest of the line */
        for(ii=xStart+hKern/2+1; ii<xEnd-hKern/2; ii++) {
          sum += (cSum[ii+hKern/2]-cSum[ii-hKern/2-1]);
          next[ii] = sum;
        }
        /* update disparity results buffer */
        for(xx=xStart+hKern/2; xx<xEnd-hKern/2; xx++)
          if(result[(hKern/2)*width+xx].best > next[xx] ||
             result[(hKern/2)*width+xx].best == VW_STEREO_MISSING_PIXEL) {
            result[(hKern/2)*width+xx].best = next[xx];
            result[(hKern/2)*width+xx].hDisp = ds;
            result[(hKern/2)*width+xx].vDisp = dsy;
          }
        /* do the rest of the image */
        for(yy=vKern/2+1, yTop = 0, yBtm = vKern*width; 
            yy<height-vKern/2; yy++, yBtm+=width, yTop+=width){
          /* update cSum buffer - remove top kernel row and add bottom one */
          for(xx=xStart; xx<xEnd; xx++)
            cSum[xx] += (diff[yBtm+xx]-diff[yTop+xx]);
          /* compute soads for the current row (yy) */
          /*getSoad2D(next, cSum, xStart, xEnd, hKern); inlined */
          /* do the first point */
          sum=0;
          for (ii=xStart; ii<xStart+hKern; ii++)
            sum += cSum[ii];
          next[xStart+hKern/2] = sum;
          /* do the rest of the line */
          for(ii=xStart+hKern/2+1; ii<xEnd-hKern/2; ii++) {
            sum += (cSum[ii+hKern/2]-cSum[ii-hKern/2-1]);
            next[ii] = sum;
          }
          /* update disparity results buffer */
          for(xx=xStart+hKern/2; xx<xEnd-hKern/2; xx++) {
            if(result[yy*width+xx].best > next[xx] ||
               result[yy*width+xx].best == VW_STEREO_MISSING_PIXEL) {
              result[yy*width+xx].best = next[xx];
              assert(ds > -1000 && ds < 1000);
              result[yy*width+xx].hDisp = ds;	  
              result[yy*width+xx].vDisp = dsy;	  
            }
          }
        }
      }
    }
    free(diff);
    free(next);
    free(cSum);

    set_progress_string(str(boost::format("H: [%1%,%2%] V: [%3%,%4%] processed succefully.")
                            % topDisp % btmDisp % minDisp % maxDisp));

    double duration= Time()-before;
    double mdisp_per_sec= ((double)numCorrTries * width * height)/duration/1e6;
    set_correlation_rate_string(str(boost::format("%1% seconds (%2% M disparities/second)")
                                    % (int)duration % mdisp_per_sec));
    return(result);
}

void CorrelationWorkThread::terminate() {
  boost::mutex::scoped_lock lock(m_mutex);
  m_should_terminate = true;
}

std::string CorrelationWorkThread::correlation_rate_string() {
  boost::mutex::scoped_lock lock(m_mutex);
  return m_correlation_rate_string;
}

std::string CorrelationWorkThread::progress_string() {
  boost::mutex::scoped_lock lock(m_mutex);
  return m_progress_string;
}

bool CorrelationWorkThread::is_done() {
  boost::mutex::scoped_lock lock(m_mutex);
  return m_done;
}
  
void CorrelationWorkThread::set_correlation_rate_string(std::string str) {
  boost::mutex::scoped_lock lock(m_mutex);
  m_correlation_rate_string = str;
}
  
void CorrelationWorkThread::set_progress_string(std::string str) {
  boost::mutex::scoped_lock lock(m_mutex);
  m_progress_string = str;
}

void CorrelationWorkThread::set_done(bool val) {
  boost::mutex::scoped_lock lock(m_mutex);
  m_done = val;
}


/********************************************************************
 *                     Constructors                                 *
 *******************************************************************/
SubpixelCorrelator::SubpixelCorrelator()
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

SubpixelCorrelator::SubpixelCorrelator(int minH,	/* left bound disparity search window*/
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

void SubpixelCorrelator::register_worker_thread(int id, CorrelationWorkThread* child) {
  boost::mutex::scoped_lock lock(m_mutex);
  m_children[id] = child;
}

CorrelationWorkThread* SubpixelCorrelator::get_worker_thread(int id) {
  boost::mutex::scoped_lock lock(m_mutex);
  return m_children[id];
}


template <class ChannelT>
ImageView<PixelDisparity<float> > SubpixelCorrelator::correlation_2D( ImageView<ChannelT> &left_image,
                                                                      ImageView<ChannelT> &right_image, 
                                                                      bool bit_image) {

  VW_ASSERT(left_image.cols() == right_image.cols() &&
            left_image.rows() == right_image.rows(),
            ArgumentErr() << "subpixel_correlation: input image dimensions do not agree.\n");


  int width = left_image.cols();
  int height = left_image.rows();
      
  //Run the correlator and record how long it takes to run.
  double begin__ = Time();
  
  uint8* bit_img0 = (uint8*)&(left_image(0,0));
  uint8* bit_img1 = (uint8*)&(right_image(0,0));
  float* pBuffer0 = (float*)&(left_image(0,0));
  float* pBuffer1 = (float*)&(right_image(0,0));

  try {
    boost::thread_group threads;
    if (bit_image == true) {
      std::cout << "\tUsing bit image optimized correlator\n";
      threads.create_thread(CorrelationWorkThread(this, 0, -m_lMaxH, -m_lMinH, -m_lMaxV, -m_lMinV, width, height, m_lKernWidth, m_lKernHeight, bit_img1, bit_img0));
      threads.create_thread(CorrelationWorkThread(this, 1, m_lMinH,  m_lMaxH,  m_lMinV,  m_lMaxV, width, height, m_lKernWidth, m_lKernHeight, bit_img0, bit_img1));
    } else {
      std::cout << "\tUsing standard correlator\n";
      threads.create_thread(CorrelationWorkThread(this, 0, -m_lMaxH, -m_lMinH, -m_lMaxV, -m_lMinV, width, height, m_lKernWidth, m_lKernHeight, pBuffer1, pBuffer0));
      threads.create_thread(CorrelationWorkThread(this, 1,  m_lMinH,  m_lMaxH,  m_lMinV,  m_lMaxV, width, height, m_lKernWidth, m_lKernHeight, pBuffer0, pBuffer1));
    }
    
    // Wait until the child threads register with their parent, then begin to query them 
    // for information about their current state.
    while (!get_worker_thread(0) || !get_worker_thread(1));
    while (!get_worker_thread(0)->is_done() || !get_worker_thread(1)->is_done()) {
      std::cout << boost::format("\t%1$-50s %2$-50s               \r") 
        % get_worker_thread(0)->progress_string() 
        % get_worker_thread(1)->progress_string();
      fflush(stdout);
      millisecond_sleep(100);
    }

    // Once the threads have complete their task, they hang around so
    // that we can print out some summary statistics.  We then
    // terminate them and join them back into the main thread.
    std::cout << boost::format("\t%1$-50s %2$-50s               \n") 
      % get_worker_thread(0)->progress_string() 
      % get_worker_thread(1)->progress_string();
    std::cout << boost::format("\t%1$-50s %2$-50s \n") 
      % get_worker_thread(0)->correlation_rate_string()
      % get_worker_thread(1)->correlation_rate_string();

    // Ask the worker threads for the actual results of the disparity correlation
    ImageView<PixelDisparity<float> > resultR2L = get_worker_thread(0)->result();
    ImageView<PixelDisparity<float> > resultL2R = get_worker_thread(1)->result();

    // We now have all the information we need.  Shut down the worker
    // threads and wait for them to finish terminating.
    get_worker_thread(0)->terminate();
    get_worker_thread(1)->terminate();
    threads.join_all(); 
  
    // Cross check the left and right disparity maps
    cross_corr_consistency_check(resultL2R, resultR2L,m_crossCorrThreshold);

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
    }

    double score = (100.0 * matched) / total;
    printf("\tCorrelation rate: %6.4f\n\n", score);

    return resultL2R;
    
  } catch (boost::thread_resource_error &e) {
    throw LogicErr() << "SubpixelCorrelator: Could not create correlation threads.\n";
  }
  
}

// Explicit Intantiation
namespace vw {
namespace stereo {

template
ImageView<PixelDisparity<float> > SubpixelCorrelator::correlation_2D<float>( ImageView<float> &left_image,
                                                                             ImageView<float> &right_image, 
                                                                             bool bit_image);
    
template
ImageView<PixelDisparity<float> > SubpixelCorrelator::correlation_2D<uint8>( ImageView<uint8> &left_image,
                                                                             ImageView<uint8> &right_image, 
                                                                             bool bit_image);

}}

