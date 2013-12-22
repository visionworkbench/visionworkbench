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


#include <vw/Stereo/DisparityMap.h>

namespace vw {
namespace stereo {

  StdDevImageFunc::StdDevImageFunc(int32 kernel_width, int32 kernel_height) :
    m_kernel_width(kernel_width), m_kernel_height(kernel_height) {
    VW_ASSERT(m_kernel_width > 0 && m_kernel_height > 0,
              ArgumentErr() << "StdDevImageFunc: kernel sizes must be non-zero.");
  }


  BBox2i StdDevImageFunc::work_area() const {
    return BBox2i(Vector2i(-m_kernel_width/2, -m_kernel_height/2),
                  Vector2i(m_kernel_width, m_kernel_height));
  }

}}    // namespace vw::stereo
