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
  // all downstream plane math (signed_dist_to_plane, rayPlaneIntersect,
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

// Sample 3 points on a proj-space plane, unproject to ECEF, and fit a local
// ECEF tangent plane through them. Accurate within ~20 m (flat-Earth limit);
// lets the Newton-Raphson refraction solver stay in ECEF, avoiding expensive
// per-iteration proj/unproj round-trips. offset_meters is the metric spacing
// between the three sample points in the plane. Callers must pass a
// meter-scale georef (typically BathyPlane::stereographic_proj) so this
// interpretation is physically meaningful.
std::vector<double> fitLocalEcefPlaneToProjPlane(
    std::vector<double> const& plane,
    vw::cartography::GeoReference const& plane_proj,
    vw::Vector3 const& proj_pt,
    double offset_meters);

// Fit a local ECEF tangent plane at proj_pt by bilinear-sampling three
// raster heights one pixel apart (center, x-neighbor, y-neighbor, with
// +/-1 fallback near edges). Caller must supply a bp with a non-empty
// water_surface. The raster-pixel spacing is the natural sampling scale,
// so offset_meters is only forwarded to the plane-based fallback path.
// Falls back if the center pixel is out of bounds or any sample is invalid.
std::vector<double> fitLocalEcefPlaneToProjSurface(BathyPlane const& bp,
                                                   vw::Vector3 const& proj_pt,
                                                   double offset_meters);

// Dispatcher: if bp carries a raster water surface, fit the tangent from it;
// otherwise fit from the global plane coefficients in bp.bathy_plane.
std::vector<double> fitLocalEcefPlane(BathyPlane const& bp,
                                      vw::Vector3 const& proj_pt,
                                      double offset_meters);

// Fit a plane in stereographic_proj's meter-scale frame from three
// raster neighbors near proj_pt (which is itself in stereographic_proj's
// coords). Used by the camera-to-ground refraction path to replace the
// global best-fit plane with a local raster-derived plane right around
// the ray-surface hit. Returns false (and leaves plane empty) if the
// center pixel is out of bounds, if both +1 and -1 neighbors are out of
// bounds on any axis, if any sample is invalid, or if the three samples
// are degenerate. Callers fall back to the global plane in that case.
// Requires bp.water_surface to be non-empty.
bool refineLocalPlaneFromRaster(BathyPlane const& bp,
                                vw::Vector3 const& proj_pt,
                                std::vector<double>& plane);

// Project an ECEF point to local projection coordinates.
vw::Vector3 proj_point(vw::cartography::GeoReference const& projection,
                       vw::Vector3 const& xyz);

// Unproject from local projection coordinates back to ECEF.
vw::Vector3 unproj_point(vw::cartography::GeoReference const& projection,
                         vw::Vector3 const& proj_pt);

// Given an ECEF point xyz and two bathy planes, find if xyz is above or
// below each plane. Outputs distances[0] and distances[1] in the same
// stereographic frame as the corresponding bathy_plane coefs.
void signed_distances_to_planes(std::vector<BathyPlane> const& bathy_plane_vec,
                                vw::Vector3 const& xyz,
                                std::vector<double>& distances);

// Intersect a ray (in ECEF) with a curved bathy plane. The plane is flat
// in the given local stereographic projection, so the intersection is
// iterated to stay both on the ray in ECEF and on the plane in proj coords.
// Outputs the intersection in ECEF, the same point in proj coords, and the
// ray direction in proj coords. mean_height is the physical water-surface
// height in meters above the datum (callers should pass
// vw::BathyPlane::mean_height).
bool rayBathyPlaneIntersect(vw::Vector3 const& in_ecef,
                            vw::Vector3 const& in_dir,
                            std::vector<double> const& plane,
                            vw::cartography::GeoReference const& plane_proj,
                            double mean_height,
                            vw::Vector3& intersect_ecef,
                            vw::Vector3& intersect_proj_pt,
                            vw::Vector3& intersect_proj_dir);

// Bend an ECEF ray at the bathy water surface and return the bent ECEF
// ray. The water surface is whatever bp describes: the global best-fit
// plane (text input or text-equivalent raster) or a per-pixel raster.
// When bp carries a raster, this routine first locates the approximate
// hit using the global plane, samples three raster neighbors there to
// fit a local plane, and bends with that refined local plane (falling
// back to the global plane on any sampling failure). For near-planar
// rasters the refined plane is numerically indistinguishable from the
// global plane, so output matches the no-raster behavior.
bool curvedSnellLaw(vw::Vector3 const& in_ecef,
                    vw::Vector3 const& in_dir,
                    BathyPlane const& bp,
                    double refraction_index,
                    vw::Vector3& out_ecef,
                    vw::Vector3& out_dir);

} // namespace vw

#endif // __VW_CARTOGRAPHY_BATHYDATA_H__
