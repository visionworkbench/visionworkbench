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


/// \file LocalToneMap.h
///
/// Tonemapping operators that compute the tonemapping in one small
/// image region at a time.
///
/// This file implements the following tone mapping operators.
///
///
#ifndef __VW_HDR_LOCALTONEMAP_H__
#define __VW_HDR_LOCALTONEMAP_H__

#include <vw/Image/ImageView.h>

namespace vw {
namespace hdr {

  ImageView<PixelRGB<double> > ashikhmin_tone_map(ImageView<PixelRGB<double> > hdr_image,
                                                  double threshold = 0.5);

}} // namespace vw::HDR

#endif  // __VW_HDR_LOCALTONEMAP_H__
