#include <vw/Image/Interpolation.h>

#if defined(VW_ENABLE_SSE) && (VW_ENABLE_SSE==1)

// This coefficients allow us to compute the bicubic weights according to (((A*norm+B)*norm+C)+D)
const float vw::bicubic_coeffs[16] __attribute__ ((aligned (16))) = 
  { -0.5, 1.5, -1.5, 0.5,  1.0, -2.5, 2.0, -0.5,  -0.5, 0.0, 0.5, 0.0,  0.0, 1.0, 0.0, 0.0 };

#endif // VW_ENABLE_SSE
