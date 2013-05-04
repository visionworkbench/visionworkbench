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


// The helper source file is used to instantian different pixel
// versions of geoblend into different object files. It works by using
// the PIXEL_TYPE macro.

#include <vw/tools/geoblend.h>

using namespace vw;

#define INSTANTIATE_CUSTOM_BLEND(PIXELTYPEMACRO, CHANNELTYPE) INSTANTIATE_CUSTOM_BLEND_(PIXELTYPEMACRO, CHANNELTYPE)

#define INSTANTIATE_CUSTOM_BLEND_(PIXELTYPEMACRO, CHANNELTYPE) INSTANTIATE_CUSTOM_BLEND__(PIXELTYPEMACRO, CHANNELTYPE, PIXELTYPEMACRO ## _ ## CHANNELTYPE)

#define INSTANTIATE_CUSTOM_BLEND__(PIXELTYPEMACRO, CHANNELTYPE, FUNCSUFFIX) \
  void do_blend_##FUNCSUFFIX(void) {                                    \
    do_blend<PIXELTYPEMACRO<CHANNELTYPE > >();                          \
  }

INSTANTIATE_CUSTOM_BLEND( PIXEL_TYPE, uint8 )
INSTANTIATE_CUSTOM_BLEND( PIXEL_TYPE, int16 )
INSTANTIATE_CUSTOM_BLEND( PIXEL_TYPE, uint16 )
INSTANTIATE_CUSTOM_BLEND( PIXEL_TYPE, float32 )
