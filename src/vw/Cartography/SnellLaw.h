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

// Straight-line Snell's law logic. Curved surface algorithms are in 
// BathyData.h.

#ifndef __VW_CARTOGRAPHY_SNELLLAW_H__
#define __VW_CARTOGRAPHY_SNELLLAW_H__

#include <vw/Math/Vector.h>

#include <vector>

namespace vw {

// Signed distance from a 3D point to a plane defined by coefficients
// [a, b, c, d] with the plane equation a*x + b*y + c*z + d = 0.
double signed_dist_to_plane(std::vector<double> const& plane,
                            vw::Vector3 const& point);

// Ray-plane intersection in a single coordinate system (whatever the plane
// coefficients and the ray are expressed in). Returns false if the ray does
// not descend towards the plane.
bool rayPlaneIntersect(vw::Vector3 const& in_xyz, vw::Vector3 const& in_dir,
                       std::vector<double> const& plane,
                       vw::Vector3 & out_xyz);

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

} // namespace vw

#endif // __VW_CARTOGRAPHY_SNELLLAW_H__
