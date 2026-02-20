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

/// \file ImageChannelRead.cc
///
/// Pre-compiled wrapper for read_channel<float> to avoid instantiating
/// the 12-type template ladder in every translation unit. The template
/// read_channel<float> in DiskImageUtils.h generates DiskImageView for
/// Vector<float, N> with N=1..12, costing ~30s of compile time per
/// caller. Keeping the instantiation here confines that cost to one
/// compilation unit in the VW library.

#include <vw/FileIO/ImageChannelRead.h>
#include <vw/FileIO/DiskImageUtils.h>

namespace vw {

ImageViewRef<float> read_float_channel(std::string const& filename, int ch) {
  return read_channel<float>(filename, ch);
}

} // namespace vw
