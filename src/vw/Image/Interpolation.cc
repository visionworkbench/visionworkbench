// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Image/Interpolation.h>

#if defined(VW_ENABLE_SSE) && (VW_ENABLE_SSE==1)

// This coefficients allow us to compute the bicubic weights according to (((A*norm+B)*norm+C)+D)
const float vw::bicubic_coeffs[16] __attribute__ ((aligned (16))) =
  { -0.5, 1.5, -1.5, 0.5,  1.0, -2.5, 2.0, -0.5,  -0.5, 0.0, 0.5, 0.0,  0.0, 1.0, 0.0, 0.0 };

#endif // VW_ENABLE_SSE
