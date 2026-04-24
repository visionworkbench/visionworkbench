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

// Snell's law plane-math primitives. Split out of BathyStereoModel.cc so
// higher-level bathy code and future raster-water-surface code can share
// the same geometric core. See SnellLaw.h for the public interface.

#include <vw/Cartography/SnellLaw.h>

#include <vw/Math/Vector.h>
#include <vw/Math/VectorUtils.h>
#include <vw/Math/LevenbergMarquardt.h>
#include <vw/Cartography/Datum.h>
#include <vw/Core/Exception.h>

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

// Compute the projected coordinates of an ECEF point
Vector3 proj_point(vw::cartography::GeoReference const& projection,
                   Vector3 const& xyz) {
  return projection.geodetic_to_point(projection.datum().cartesian_to_geodetic(xyz));
}

// Reverse this operation
Vector3 unproj_point(vw::cartography::GeoReference const& projection,
                     Vector3 const& proj_pt) {
  return projection.datum().geodetic_to_cartesian(projection.point_to_geodetic(proj_pt));
}

// Test Snell's law in projected and unprojected coordinates
static void testSnellLaw(std::vector<double> const& plane,
                         vw::cartography::GeoReference const& plane_proj,
                         double refraction_index,
                         vw::Vector3 const& out_ecef,
                         vw::Vector3 const& in_ecef_dir, vw::Vector3 const& out_ecef_dir,
                         vw::Vector3 const& out_proj_pt,
                         vw::Vector3 const& in_proj_dir, vw::Vector3 const& out_proj_dir) {

  // 1. In projected coordinates
  Vector3 proj_normal(plane[0], plane[1], plane[2]);
  double sin_in = sin(acos(dot_prod(proj_normal, -in_proj_dir)));
  double sin_out = sin(acos(dot_prod(-proj_normal, out_proj_dir)));

  // 2. In unprojected coordinates
  Vector3 proj_pt_above_normal = out_proj_pt + 1.0 * proj_normal; // go 1 m along the normal
  Vector3 ecef_above_normal = unproj_point(plane_proj,
                                           proj_pt_above_normal);
  Vector3 ecef_normal = ecef_above_normal - out_ecef;
  ecef_normal /= norm_2(ecef_normal); // normalize
  sin_in = sin(acos(dot_prod(ecef_normal, -in_ecef_dir)));
  sin_out = sin(acos(dot_prod(-ecef_normal, out_ecef_dir)));

  // Verify that the incoming ray, outgoing ray, and the
  // normal are in the same plane in projected coordinates

  // 1. In projected coordinates
  Vector3 in_out_normal = vw::math::cross_prod(in_proj_dir, out_proj_dir);
  double plane_error = dot_prod(in_out_normal, proj_normal);

  // 2. In unprojected coordinates
  in_out_normal = vw::math::cross_prod(in_ecef_dir, out_ecef_dir);
  plane_error = dot_prod(in_out_normal, ecef_normal);
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

// Find where a ray (in ECEF) intersects a curved bathy plane. The plane is
// modeled as flat in local stereographic projection. Returns the intersection
// point in ECEF, the intersection point in projected coordinates, and the ray
// direction at that location in projected coordinates.
// Returns false if intersection fails.
bool rayBathyPlaneIntersect(vw::Vector3 const& in_ecef,
                            vw::Vector3 const& in_dir,
                            std::vector<double> const& plane,
                            vw::cartography::GeoReference const& plane_proj,
                            double mean_height,
                            vw::Vector3 & intersect_ecef,
                            vw::Vector3 & intersect_proj_pt,
                            vw::Vector3 & intersect_proj_dir) {

  // Seed the ray-surface intersection from the physical mean water height,
  // supplied by the caller. -plane[3]/plane[2] would be equivalent only when
  // plane_proj is meter-scale (e.g. a local stereographic projection); for
  // geographic plane_proj it can be off by ~100 m.
  //
  // TODO(oalexan1): the per-step refinement below still intersects the ECEF
  // ray with the fitted plane in plane_proj's coords. That step works in any
  // coord system arithmetically, but interpreting the resulting point as
  // "on the water surface" assumes the plane equation describes a meter-
  // scale surface. For a geographic plane_proj the fitted plane is in
  // mixed units and convergence can be loose. See water_surface_notes.sh.
  double major_radius = plane_proj.datum().semi_major_axis() + mean_height;
  double minor_radius = plane_proj.datum().semi_minor_axis() + mean_height;

  // Intersect the ray with the mean water surface, this will give us the
  // initial guess for intersecting with that surface.
  intersect_ecef = vw::cartography::datum_intersection(major_radius, minor_radius,
                                                       in_ecef, in_dir);

  // The fact that we trace a ray below in projected coordinates, even if very
  // close to the bathy plane and very short, can still introduce some small
  // error. So refine intersect_ecef so it is both along the ray in ECEF and on the
  // curved plane.
  for (int pass = 0; pass < 5; pass++) {

    // Move a little up the ray. Move less on later passes.
    vw::Vector3 prev_ecef = intersect_ecef - 1.0 * in_dir / (1.0 + 10.0 * pass);

    // Compute projected entries. These will be exported out of this function.
    intersect_proj_pt = proj_point(plane_proj, intersect_ecef);
    vw::Vector3 prev_proj_pt = proj_point(plane_proj, prev_ecef);
    intersect_proj_dir = intersect_proj_pt - prev_proj_pt;
    intersect_proj_dir /= norm_2(intersect_proj_dir);

    // Stop when we are within 0.1 mm of the plane, while along the ray. Going
    // beyond that seems not useful. This is usually reached on second pass.

    if (std::abs(signed_dist_to_plane(plane, intersect_proj_pt)) < 1e-4)
      break;

    // Intersect the proj ray with the proj plane
    vw::Vector3 refined_intersect_proj_pt;
    if (!rayPlaneIntersect(intersect_proj_pt, intersect_proj_dir, plane,
                           refined_intersect_proj_pt))
      return false;

    // Convert back to ECEF
    intersect_ecef = unproj_point(plane_proj, refined_intersect_proj_pt);

    // Put the point back on the ray. Then it may become slightly off the plane.
    intersect_ecef = in_ecef + dot_prod(intersect_ecef - in_ecef, in_dir) * in_dir;
  }

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

// Consider a stereographic projection and a plane
// a * x + b * y + c * z + d = 0 for (x, y, z) in this projection.
// Intersect it with a ray given in ECEF coordinates.
// If the values a and b are 0, that is the same as intersecting
// the ray with the spheroid of values -d/c above the datum.
// This solver was not used as it was too slow. An approximate
// solution was instead found.
class SolveCurvedPlaneIntersection:
  public vw::math::LeastSquaresModelBase<SolveCurvedPlaneIntersection> {
  vw::Vector3 const& m_ray_pt;
  vw::Vector3 const& m_ray_dir;
  vw::cartography::GeoReference const& m_projection;
  std::vector<double> const& m_proj_plane;
public:

  // This is a one-parameter problem, yet have to use a vector (of size 1)
  // as required by the API.
  typedef vw::Vector<double, 1> result_type;   // residual
  typedef vw::Vector<double, 1> domain_type;   // parameter giving the position on the ray
  typedef vw::Matrix<double>    jacobian_type;

  /// Instantiate the solver with a set of xyz to pixel pairs and a pinhole model
  SolveCurvedPlaneIntersection(vw::Vector3 const& ray_pt, vw::Vector3 const& ray_dir,
                                vw::cartography::GeoReference const& projection,
                                std::vector<double> const& proj_plane):
    m_ray_pt(ray_pt), m_ray_dir(ray_dir), m_projection(projection),
    m_proj_plane(proj_plane) {}

  /// Given the camera, project xyz into it
  inline result_type operator()(domain_type const& t) const {

    // Get the current point along the ray
    Vector3 xyz = m_ray_pt + t[0] * m_ray_dir;

    // Convert to projected coordinates
    Vector3 proj_pt = proj_point(m_projection, xyz);

    result_type ans;
    ans[0] = signed_dist_to_plane(m_proj_plane, proj_pt);
    return ans;
  }
}; // End class SolveCurvedPlaneIntersection

// Given a ray in ECEF and a water surface which is a plane only in a local
// stereographic projection, compute how the ray bends under Snell's law. Use
// the following approximate logic. Find where the ray intersects the datum with
// the mean water height, as then it is close to the water surface, since the
// water surface is almost horizontal in projected coordinates. Find a point on
// that ray 1 m before that. Convert both of these points from ECEF to the
// projected coordinate system. Do Snell's law in that coordinate system for the
// ray going through those two projected points. Find a point on the outgoing
// ray in projected coordinates Find another close point further along it. Undo
// the projection for these two points. That will give the outgoing direction in
// ECEF.
bool curvedSnellLaw(Vector3 const& in_ecef, Vector3 const& in_dir,
                    std::vector<double> const& plane,
                    vw::cartography::GeoReference const& plane_proj,
                    double refraction_index, double mean_height,
                    Vector3 & out_ecef, Vector3 & out_dir) {

  // Intersect the ray with the curved water surface. Return the intersection
  // point in ecef, that point in projected coordinates, and the ray direction
  // at that location in projected coordinates.
  Vector3 intersect_ecef, intersect_proj_pt, intersect_proj_dir;
  if (!rayBathyPlaneIntersect(in_ecef, in_dir, plane, plane_proj, mean_height,
                              intersect_ecef, intersect_proj_pt, intersect_proj_dir))
    return false;

  // Snell's law in projected coordinates
  Vector3 out_proj_pt, out_proj_dir; // in the water
  bool ans = snellLaw(intersect_proj_pt, intersect_proj_dir,
                      plane, refraction_index,
                      out_proj_pt, out_proj_dir);

  // If Snell's law failed to work, exit early
  if (!ans)
    return ans;

  // Move a little on the ray in projected coordinates.
  Vector3 next_proj_pt = out_proj_pt + 1.0 * out_proj_dir;

  // Convert back to ECEF
  out_ecef = unproj_point(plane_proj, out_proj_pt);
  Vector3 next_ecef = unproj_point(plane_proj, next_proj_pt);

  // Finally get the outgoing direction according to Snell's law in ECEF.
  // The assumption here is that at ground level a short vector in ECEF
  // is very close to the same short vector in projected coordinates.
  out_dir = next_ecef - out_ecef;
  out_dir /= norm_2(out_dir);

#if 0
  // Sanity check
  testSnellLaw(plane,
               plane_proj,
               refraction_index,
               out_ecef, in_dir, out_dir,
               out_proj_pt, intersect_proj_dir, out_proj_dir);
#endif

  return true;
}

} // namespace vw
