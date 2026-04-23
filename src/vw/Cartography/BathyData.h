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

// Data types that describe the bathymetry inputs to the stereo pipeline:
// a per-image water-surface plane (four coefficients in a local stereographic
// projection) and a bundle of masks + planes + refraction index consumed by
// BathyStereoModel during triangulation. The types live in their own header
// so consumers that only hold or pass bathy data (e.g. option structs,
// function signatures) do not need to pull in the full BathyStereoModel.

#ifndef __VW_CARTOGRAPHY_BATHYDATA_H__
#define __VW_CARTOGRAPHY_BATHYDATA_H__

#include <vw/Cartography/GeoReference.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/PixelMask.h>

#include <vector>

namespace vw {

struct BathyPlane {
  std::vector<double> bathy_plane;
  vw::cartography::GeoReference plane_proj;
};

// A struct to hold the bathymetry settings and data
struct BathyData {
  std::vector<vw::ImageViewRef<vw::PixelMask<float>>> bathy_masks;
  std::vector<BathyPlane> bathy_planes;
  float refraction_index;
  BathyData(): refraction_index(1.0) {}
};

} // namespace vw

#endif // __VW_CARTOGRAPHY_BATHYDATA_H__
