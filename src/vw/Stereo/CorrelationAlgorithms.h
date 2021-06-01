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


#ifndef __VW_STEREO_CORRELATION_ALGORITHMS_H__
#define __VW_STEREO_CORRELATION_ALGORITHMS_H__

// A very small header file having a single definition. Put it here
// rather than in CorrelationView.h to avoid including that large
// header file everywhere.

namespace vw {
namespace stereo {

  enum CorrelationAlgorithm {
    VW_CORRELATION_BM        = 0, ///< Use a local sliding search window.
    VW_CORRELATION_SGM       = 1, ///< Use a Semi-Global Matching algorithm.
    VW_CORRELATION_MGM       = 2, ///< Use the More-Global Matching algorithm.
    VW_CORRELATION_FINAL_MGM = 3, ///< SGM until resolution level 0, then MGM.
    VW_CORRELATION_OTHER     = 4  ///< Some external 1D stereo algorithm.
  };

}} // namespace vw::stereo

#endif//__VW_STEREO_CORRELATION_ALGORITHMS_H__
