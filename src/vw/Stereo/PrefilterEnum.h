// __BEGIN_LICENSE__
//  Copyright (c) 2006-2026, United States Government as represented by the
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

#ifndef __VW_STEREO_PREFILTER_ENUM_H__
#define __VW_STEREO_PREFILTER_ENUM_H__

namespace vw {
namespace stereo {

  enum PrefilterModeType {
    PREFILTER_NONE    = 0,
    PREFILTER_MEANSUB = 1,
    PREFILTER_LOG     = 2
  };

}} // end namespace vw::stereo

#endif // __VW_STEREO_PREFILTER_ENUM_H__
