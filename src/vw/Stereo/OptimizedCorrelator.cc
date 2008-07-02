#include <vw/Stereo/OptimizedCorrelator.h>

#include <boost/format.hpp>

#include <vw/Core/Thread.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/Algorithms.h>

using namespace std;
using namespace vw;
using namespace vw::stereo;

// Some useful default values and constants 
#define DEFAULT_KERN_WIDTH 24
#define DEFAULT_KERN_HEIGHT 24
#define DEFAULT_MIN_H -50
#define DEFAULT_MAX_H 50
#define DEFAULT_MIN_V -1
#define DEFAULT_MAX_V 1
#define DEFAULT_CROSSCORR_THRESHOLD 2
#define DEFAULT_USE_SUBPIXEL_H true
#define DEFAULT_USE_SUBPIXEL_V true

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
}
