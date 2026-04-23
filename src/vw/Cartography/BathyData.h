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

#include <vw/Math/Vector.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/PixelMask.h>

#include <string>
#include <vector>

namespace vw {

// Per-image water surface used during bathymetry triangulation.
//
// Can hold either a best-fit plane in a local stereographic projection
// (four coefficients, with the projection stored in plane_proj), or a
// georeferenced image of water-surface heights. When the image is
// provided, the best-fit plane approximation in plane_proj's coordinates
// is also stored. If the image is empty, use the plane.
struct BathyPlane {
  // Plane coefficients a, b, c, d such that a*x + b*y + c*z + d = 0 in
  // plane_proj's local stereographic frame. Normal (a, b, c) points up.
  std::vector<double> bathy_plane;
  vw::cartography::GeoReference plane_proj;

  // Optional image of water-surface heights.
  vw::ImageView<vw::PixelMask<float>> water_surface;
};

// A struct to hold the bathymetry settings and data
struct BathyData {
  std::vector<vw::ImageViewRef<vw::PixelMask<float>>> bathy_masks;
  std::vector<BathyPlane> bathy_planes;
  float refraction_index;
  BathyData(): refraction_index(1.0) {}
};

// Check if the given left and right pixels are in the masked region (invalid in
// the mask). That will mean bathymetry correction should be applied.
bool areMasked(ImageViewRef<PixelMask<float>> const& left_mask,
               ImageViewRef<PixelMask<float>> const& right_mask,
               Vector2 const& lpix, Vector2 const& rpix);

// Read a bathy mask. Water pixels are those with non-positive values or
// matching the file's nodata value. Both are invalidated in the returned
// masked image. The returned nodata_val is suitable for writing the mask back.
vw::ImageViewRef<vw::PixelMask<float>> read_bathy_mask(std::string const& filename,
                                                       float & nodata_val);

// Read a set of bathy masks.
void read_bathy_masks(std::vector<std::string> const& mask_filenames,
                      std::vector<vw::ImageViewRef<vw::PixelMask<float>>> & bathy_masks);

// Read the bathy planes and associated data. More often than not they will be
// identical. If there is more than one bathy plane file, they are all kept in
// the same string, separated by space.
void readBathyPlanes(std::string const& bathy_plane_files,
                     int num_images,
                     std::vector<BathyPlane> & bathy_plane_vec);

// Best-fit plane through a set of 3D points via SVD (no outlier rejection).
// Returns the four coefficients (a, b, c, d) such that a*x + b*y + c*z + d = 0,
// with the normal (a, b, c) of unit length and oriented so c > 0.
void fitPlaneToPoints(std::vector<vw::Vector3> const& points,
                      std::vector<double> & bathy_plane);

} // namespace vw

#endif // __VW_CARTOGRAPHY_BATHYDATA_H__
