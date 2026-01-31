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

#ifndef __VW_CARTOGRAPHY_HILLSHADE_H__
#define __VW_CARTOGRAPHY_HILLSHADE_H__

#include <vw/FileIO/GdalWriteOptions.h>
#include <string>

/// \file Hillshade.h. Hillshade images with georeference.

namespace vw {
namespace cartography {

  /// Redirect to the function with the required data type.
  void do_multitype_hillshade(std::string const& input_file,
                              std::string const& output_file,
                              double azimuth, double elevation, double scale,
                              double nodata_value, double blur_sigma,
                              bool align_to_georef,
                              vw::GdalWriteOptions const& opt);

}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_HILLSHADE_H__

