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

// Straight-line Snell's law primitives. See SnellLaw.h for the public
// interface. Anything bathy-specific (water surface as a plane in a
// local stereographic projection, raster water surface, ECEF / proj
// round-trips) lives in vw/Cartography/BathyData.cc.

#include <vw/Cartography/SnellLaw.h>

#include <vw/Math/Vector.h>

#include <cmath>

namespace vw {

// Given a plane as four values a, b, c, d, with the plane being
// a * x + b * y + c * z + d = 0, find how far off a point (x, y, z) is from the plane
// by evaluating the above expression.
double signed_dist_to_plane(std::vector<double> const& plane, vw::Vector3 const& point) {

  double ans = 0.0;
  for (unsigned coord_it = 0; coord_it < 3; coord_it++) {
    ans += plane[coord_it] * point[coord_it];
  }
  ans += plane[3];

  return ans;
}

// Find where a ray given by in_xyz + alpha * in_dir intersects a plane.
// Returns false if the ray does not descend towards the plane.
bool rayPlaneIntersect(vw::Vector3 const& in_xyz, vw::Vector3 const& in_dir,
                       std::vector<double> const& plane,
                       vw::Vector3 & out_xyz) {

  // See where the ray intersects the plane
  double cn = 0.0, dn = 0.0; // Dot product of in_xyz and in_dir with plane normal n
  for (size_t it = 0; it < 3; it++) {
    cn += plane[it] * in_xyz[it];
    dn += plane[it] * in_dir[it];
  }

  // The ray must descend to the plane, or else something is not right
  if (dn >= 0.0)
    return false;

  double alpha = -(plane[3] + cn)/dn;

  // The intersection with the plane
  out_xyz = in_xyz + alpha * in_dir;

  return true;
}

// See the .h file for more info
bool snellLaw(Vector3 const& in_xyz, Vector3 const& in_dir,
              std::vector<double> const& plane,
              double refraction_index,
              Vector3 & out_xyz, Vector3 & out_dir) {

  // Find where the ray intersects the plane
  if (!rayPlaneIntersect(in_xyz, in_dir, plane, out_xyz))
    return false;

  // Compute dn for Snell's law calculation
  double dn = 0.0;
  for (size_t it = 0; it < 3; it++)
    dn += plane[it] * in_dir[it];

  // Let n be the plane normal pointing up (the first three components
  // of the plane vector). Let out_dir be the outgoing vector after the
  // ray hits the water, according to Snell's law, with in_dir being the
  // incoming ray. Let a1 be the angles between -in_dir and n, a2 be the
  // angle between out_dir and -n.

  // Then sin(a1) = refraction_index * sin(a2) per Snell's law.
  // Square this. Note that cos^2 (x) + sin^2 (x) = 1.
  // So, 1 - cos(a1)^2 = refraction_index^2 * (1 - cos(a2)^2).
  // But cos(a1) = dot_product(-in_dir, n) = -dn.
  // So, cos(a2)^2 = 1 - (1 - dn^2)/refraction_index^2
  // Call the left-hand value cos_sq.

  double cos_sq = 1.0 - (1.0 - dn * dn)/refraction_index/refraction_index;

  // The outgoing vector out_dir will be a linear combination of -n and d1,
  // normalized to unit length. Let alpha > 0 be the value which will
  // produce the linear combination.  So,
  // out_dir = (-n + alpha * in_dir)/norm(-n + alpha * in_dir)
  // But dot(out_dir, -n) = cos(a2). Hence, if we dot the above with n and square it,
  // we get
  // cos(a2)^2 = (-1 + alpha * dn)^2 / dot( -n + alpha * in_dir, -n + alpha * in_dir)
  // or
  // cos(a2)^2 * dot( -n + alpha * in_dir, -n + alpha * in_dir) = (-1 + alpha * dn)^2
  // or
  // cos_sq * (1 - 2 * alpha * dn + alpha^2) = ( 1 - 2*alpha * dn + alpha^2 * dn^2)
  //
  // Note that we computed cos_sq from Snell's law above.

  // Move everything to the left and find the coefficients of the
  // quadratic equation in alpha, so u * alpha^2 + v * alpha + w = 0.
  double u = cos_sq - dn * dn;  // this is cos(a2)^2 - cos(a1)^2 > 0 as a2 < a1
  double v = -2 * dn * cos_sq + 2.0 * dn;
  double w = cos_sq - 1.0;
  double delta = v * v - 4 * u * w; // discriminant
  if (u <= 0.0 || delta < 0.0)
    return false; // must not happen

  double alpha = (-v + sqrt(delta)) / (2.0 * u); // pick the positive quadratic root

  if (alpha < 0)
    return false; // must not happen

  // The normalized direction after the ray is bent
  out_dir = -Vector3(plane[0], plane[1], plane[2]) + alpha * in_dir;
  out_dir = out_dir / norm_2(out_dir);

  return true;
}

} // namespace vw
