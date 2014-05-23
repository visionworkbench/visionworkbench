// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#include <vw/Image/Interpolation.h>
/*

// This coefficients allow us to compute the bicubic weights according to (((A*norm+B)*norm+C)+D)
const float vw::bicubic_coeffs[16] __attribute__ ((aligned (16))) =
  { -0.5, 1.5, -1.5, 0.5,  1.0, -2.5, 2.0, -0.5,  -0.5, 0.0, 0.5, 0.0,  0.0, 1.0, 0.0, 0.0 };

#endif // VW_ENABLE_SSE
*/
