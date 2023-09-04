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


/// \file BlockRasterize.cc
///
/// A little shared utility used by both BlockRasterizeView.h and  
/// BlockImageOperator.h

#include <vw/Image/BlockRasterize.h>

namespace vw {
  // GDAL images are sometimes stored with very tall blocks, on the order
  // of 200,000 pixels. That makes VW fail and/or slow. Adjust the 
  // size of these blocks when they are used for internal tiling
  // (which need not be the same as what is on disk).
  void adjust_block_size(vw::Vector2i & block_size) {
    block_size[0] = std::min(block_size[0], 5120);
    block_size[1] = std::min(block_size[1], 5120);
  }
} // end namespace vw
