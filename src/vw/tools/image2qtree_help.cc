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

// The helper source file is used to instantiate different pixel
// versions of image2qtree into different object files. It works by using
// the PIXEL_TYPE macro.

#include <vw/tools/image2qtree.h>

namespace vw {
  
#define INSTANTIATE_CUSTOM_MOSAIC(PIXELTYPEMACRO, CHANNELTYPE) INSTANTIATE_CUSTOM_MOSAIC_(PIXELTYPEMACRO, CHANNELTYPE)

#define INSTANTIATE_CUSTOM_MOSAIC_(PIXELTYPEMACRO, CHANNELTYPE) INSTANTIATE_CUSTOM_MOSAIC__(PIXELTYPEMACRO, CHANNELTYPE, PIXELTYPEMACRO ## _ ## CHANNELTYPE)

#define INSTANTIATE_CUSTOM_MOSAIC__(PIXELTYPEMACRO, CHANNELTYPE, FUNCSUFFIX) \
  void do_mosaic_##FUNCSUFFIX(const Options& opt,                       \
                              const vw::ProgressCallback *progress) {   \
    do_mosaic<vw::PIXELTYPEMACRO<vw::CHANNELTYPE > >(opt, progress);    \
  }

INSTANTIATE_CUSTOM_MOSAIC( PIXEL_TYPE, CHANNEL_TYPE )

} // end namespace vw
