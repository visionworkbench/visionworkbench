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

// Bathymetry data types and a few low-level primitives. Holds the
// vw::BathyPlane struct (water-surface plane coefficients + optional
// raster + companion stereographic georef + mean height), the
// vw::BathyData bundle (masks + planes + refraction index), the mask /
// plane file readers, and the ECEF <-> projected-coordinates helpers
// (bathyProjPoint, bathyUnprojPoint). The actual ray-bending logic
// (curvedSnellLaw, etc.) lives in vw/Cartography/BathyRay.h.

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
// (four coefficients, with the projection stored in stereographic_proj),
// or a georeferenced image of water-surface heights. When the image is
// provided, the best-fit plane approximation in stereographic_proj's
// coordinates is also stored. If the image is empty, use the plane.
struct BathyPlane {
  // Plane coefficients a, b, c, d such that a*x + b*y + c*z + d = 0 in
  // stereographic_proj's meter-scale frame. Normal (a, b, c) points up.
  std::vector<double> bathy_plane;

  // Coordinate frame of the stored raster. For text input this is a
  // local stereographic projection. For raster input this is whatever
  // georef the input file carried (may be geographic lon/lat). Use it
  // for raster pixel lookups (pixel_to_point, point_to_pixel) only.
  vw::cartography::GeoReference plane_proj;

  // Always a local stereographic projection with meter-scale x, y axes.
  // Identical to plane_proj when plane_proj is stereographic (text
  // input, or a stereographic raster). For a geographic raster input,
  // derived at load time with the origin at the raster centroid, so
  // all downstream plane math (signedDistToPlane, rayPlaneIntersect,
  // Snell's law in proj coords, local-tangent fallback) can assume
  // meter-scale axes by construction.
  vw::cartography::GeoReference stereographic_proj;

  // Optional image of water-surface heights.
  vw::ImageView<vw::PixelMask<float>> water_surface;

  // Mean water-surface height in meters above the datum. Computed at load
  // time from whichever input was supplied: the plane coefficients (text)
  // or the raster pixel values (image).
  double mean_height = 0.0;
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

// Read a single bathy plane. Accepts either a text plane file or a GeoTIFF
// raster of water surface heights.
void readBathyPlane(std::string const& bathy_plane_file, BathyPlane & bp);

// Read the bathy planes and associated data. More often than not they will be
// identical. If there is more than one bathy plane file, they are all kept in
// the same string, separated by space.
void readBathyPlanes(std::string const& bathy_plane_files,
                     int num_images,
                     std::vector<BathyPlane> & bathy_plane_vec);

// Project an ECEF point to local projection coordinates using the
// projection's datum. Bathy code is the primary user; nothing here is
// bathy-specific in implementation.
vw::Vector3 bathyProjPoint(vw::cartography::GeoReference const& projection,
                           vw::Vector3 const& xyz);

// Inverse of bathyProjPoint: from local projection coordinates back to ECEF.
vw::Vector3 bathyUnprojPoint(vw::cartography::GeoReference const& projection,
                             vw::Vector3 const& proj_pt);

// Given an ECEF point xyz and two bathy planes, find if xyz is above or
// below each plane. Outputs distances[0] and distances[1] in the same
// stereographic frame as the corresponding bathy_plane coefs.
void signedDistToPlanes(std::vector<BathyPlane> const& bathy_plane_vec,
                        vw::Vector3 const& xyz,
                        std::vector<double>& distances);

} // namespace vw

#endif // __VW_CARTOGRAPHY_BATHYDATA_H__
