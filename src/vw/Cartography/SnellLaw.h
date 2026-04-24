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

// Snell's law refraction math for an ocean surface modeled as a plane
// (flat in a local stereographic projection; curved in ECEF). This file
// holds only the geometric primitives. Nothing here depends on the
// vw::BathyPlane struct; callers pass the four plane coefficients and a
// GeoReference directly. See vw/Cartography/BathyStereoModel.h for the
// higher-level types (BathyPlane, BathyData) and the BathyStereoModel
// class that consumes these primitives.

#ifndef __VW_CARTOGRAPHY_SNELLLAW_H__
#define __VW_CARTOGRAPHY_SNELLLAW_H__

#include <vw/Math/Vector.h>
#include <vw/Cartography/GeoReference.h>

#include <vector>

namespace vw {

// Signed distance from a 3D point to a plane defined by coefficients
// [a, b, c, d] with the plane equation a*x + b*y + c*z + d = 0.
double signed_dist_to_plane(std::vector<double> const& plane,
                            vw::Vector3 const& point);

// Project an ECEF point to local projection coordinates.
vw::Vector3 proj_point(vw::cartography::GeoReference const& projection,
                       vw::Vector3 const& xyz);

// Unproject from local projection coordinates back to ECEF.
vw::Vector3 unproj_point(vw::cartography::GeoReference const& projection,
                         vw::Vector3 const& proj_pt);

// Ray-plane intersection in a single coordinate system (whatever the plane
// coefficients and the ray are expressed in). Returns false if the ray does
// not descend towards the plane.
bool rayPlaneIntersect(vw::Vector3 const& in_xyz, vw::Vector3 const& in_dir,
                       std::vector<double> const& plane,
                       vw::Vector3 & out_xyz);

// Intersect a ray (in ECEF) with a curved bathy plane. The plane is flat
// in the given local stereographic projection, so the intersection is
// iterated to stay both on the ray in ECEF and on the plane in proj coords.
// Outputs the intersection in ECEF, the same point in proj coords, and the
// ray direction in proj coords. mean_height is the physical water-surface
// height in meters above the datum, used to seed the initial ray-datum
// intersection (callers should pass vw::BathyPlane::mean_height).
bool rayBathyPlaneIntersect(vw::Vector3 const& in_ecef,
                            vw::Vector3 const& in_dir,
                            std::vector<double> const& plane,
                            vw::cartography::GeoReference const& plane_proj,
                            double mean_height,
                            vw::Vector3 & intersect_ecef,
                            vw::Vector3 & intersect_proj_pt,
                            vw::Vector3 & intersect_proj_dir);

// Given a ray going down towards Earth, starting at point in_xyz and
// with unit direction in_dir, a plane 'p' to the water surface with four
// coefficients such that the plane equation is p[0] * x + p[1] * y
// + p[2] * z + p[3] = 0, the normal (p[0], p[1], p[2]) pointing
// upwards away from Earth, and water refraction index, find where
// this ray meets the water plane named out_xyz, and the ray direction out_dir
// after it bends according to Snell's law. Return true on success.
// This also works in projected coordinates.
bool snellLaw(vw::Vector3 const& in_xyz, vw::Vector3 const& in_dir,
              std::vector<double> const& plane, double refraction_index,
              vw::Vector3 & out_xyz, vw::Vector3 & out_dir);

// Like snellLaw, but for a curved water surface. The water surface is
// modeled as a plane in local stereographic projection coordinates. The
// ray is bent in that coordinate system, then transformed back to ECEF.
// mean_height is the physical water-surface height in meters above the
// datum; callers should pass vw::BathyPlane::mean_height.
bool curvedSnellLaw(vw::Vector3 const& in_ecef, vw::Vector3 const& in_dir,
                    std::vector<double> const& plane,
                    vw::cartography::GeoReference const& plane_proj,
                    double refraction_index, double mean_height,
                    vw::Vector3 & out_ecef, vw::Vector3 & out_dir);

} // namespace vw

#endif // __VW_CARTOGRAPHY_SNELLLAW_H__
