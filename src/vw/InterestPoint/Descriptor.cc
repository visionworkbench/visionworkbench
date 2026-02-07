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

/// \file Descriptor.cc
///
/// Basic classes and structures for storing image interest points.
///
#include <vw/InterestPoint/Descriptor.h>

namespace vw {
namespace ip {

  const uint32 SGradDescriptorGenerator::box_strt[5] = {17,15,9,6,0};
  const uint32 SGradDescriptorGenerator::box_size[5] = {2,4,8,10,14};
  const uint32 SGradDescriptorGenerator::box_half[5] = {1,2,4,5,7};

}} // namespace vw::ip
