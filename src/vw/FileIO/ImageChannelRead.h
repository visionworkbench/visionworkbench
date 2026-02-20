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

/// \file ImageChannelRead.h

#ifndef __VW_FILEIO_IMAGECHANNELREAD_H__
#define __VW_FILEIO_IMAGECHANNELREAD_H__

#include <vw/Image/ImageViewRef.h>
#include <vw/Math/Vector.h>
#include <string>

namespace vw {

// Pre-compiled wrappers to avoid instantiating the 12-type template
// ladder (read_channels in DiskImageUtils.h) in every translation unit.
// Each read_channels call generates DiskImageView<Vector<T, N>> for
// N=1..12, costing ~30s of compile time per caller.

// Read a single channel as float. Channel starts at 0.
ImageViewRef<float> read_float_channel(std::string const& filename, int ch);

// Read 3/4/6 channels as double. Channel starts at 0.
ImageViewRef<Vector<double, 3>>
read_double_channels_3(std::string const& filename, int ch);
ImageViewRef<Vector<double, 4>>
read_double_channels_4(std::string const& filename, int ch);
ImageViewRef<Vector<double, 6>>
read_double_channels_6(std::string const& filename, int ch);

} // namespace vw

#endif
