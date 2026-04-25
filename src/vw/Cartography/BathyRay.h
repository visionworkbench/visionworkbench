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

// Bathy ray-water-surface logic. Curved Snell's law in ECEF against a
// vw::BathyPlane water surface (formula or per-pixel raster), plus the
// few helpers that other modules need: per-plane signed distances and
// the local-ECEF-tangent fit used by the camera-to-pixel solver.
// Everything else - the curved-plane intersection iterator, the local
// tangent variants, raster neighbor sampling - lives file-local in
// BathyRay.cc; consumers do not need it. Low-level ECEF / projected
// helpers (bathyProjPoint, bathyUnprojPoint) live in BathyData.h.
// Straight-line Snell's law is in SnellLaw.h. Camera-aware helpers
// (point_to_pixel, datumBathyIntersection) are in BathyCamera.h.
// BathyStereoModel is the triangulator.

#ifndef __VW_CARTOGRAPHY_BATHYRAY_H__
#define __VW_CARTOGRAPHY_BATHYRAY_H__

#include <vw/Math/Vector.h>
#include <vw/Cartography/BathyData.h>

#include <vector>

namespace vw {

// Given an ECEF point xyz and two bathy planes, find if xyz is above or
// below each plane. Outputs distances[0] and distances[1] in the same
// stereographic frame as the corresponding bathy_plane coefs.
void signed_distances_to_planes(std::vector<BathyPlane> const& bathy_plane_vec,
                                vw::Vector3 const& xyz,
                                std::vector<double>& distances);

// Fit a local ECEF tangent plane near proj_pt. Raster-aware: when bp
// carries a raster water surface, the tangent reflects the raster near
// proj_pt (3 pixel-aligned samples through ECEF); otherwise it reflects
// the global fitted plane. proj_pt is in bp.stereographic_proj's
// meter-scale frame. offset_meters is forwarded to the plane-based path.
std::vector<double> fitLocalEcefPlane(BathyPlane const& bp,
                                      vw::Vector3 const& proj_pt,
                                      double offset_meters);

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

#endif // __VW_CARTOGRAPHY_BATHYRAY_H__
