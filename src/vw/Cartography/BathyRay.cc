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

// Bathy ray-water-surface logic. See BathyRay.h for the public contract.
// The internal helpers (curved-plane intersection iterator, local
// tangent variants, raster neighbor sampling, the private
// curvedSnellLawAgainstPlane that does the actual ray bend) live in
// the anonymous namespace below.

#include <vw/Cartography/BathyRay.h>

#include <vw/Cartography/BathyData.h>
#include <vw/Cartography/Datum.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/SnellLaw.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/Interpolation.h>
#include <vw/Image/PixelMask.h>
#include <vw/Math/Vector.h>
#include <vw/Math/VectorUtils.h>
#include <vw/Core/Exception.h>

#include <cmath>

namespace vw {

namespace {

// True when (px[0], px[1]) is a valid bilinear-sample location inside an
// image of size (cols, rows). The upper bound is inclusive; bilinear
// interpolation at an exact integer coordinate collapses to a direct pixel
// read and does not access an out-of-bounds neighbor.
bool isInBounds(vw::Vector2 const& px, int cols, int rows) {
  return px[0] >= 0.0 && px[0] <= double(cols - 1) &&
         px[1] >= 0.0 && px[1] <= double(rows - 1);
}

// One bilinear-sampled raster pixel: its (possibly fractional) pixel
// coordinates and the height value at that location.
struct LocalRasterSample {
  vw::Vector2 pix;
  double height;
};

// Pick three pixel-aligned raster samples near proj_pt: the center pixel,
// an x-neighbor (+1 with -1 fallback at the high edge), and a y-neighbor
// (+1 with -1 fallback). Bilinear-sample each and return them. Sets
// samples[0] = center, [1] = x-neighbor, [2] = y-neighbor.
//
// proj_pt must be in bp.stereographic_proj's frame (meter-scale). Internally
// converts to bp.plane_proj (the raster's native georef) for the pixel
// lookup, so geographic and stereographic rasters both work.
//
// Returns false if the center pixel is out of bounds, if both +1 and -1
// neighbors are out of bounds on either axis, or if any of the three
// bilinear-interp values is invalid.
// TODO(oalexan1): Make this scheme use centered samples
bool pickLocalRasterSamples(BathyPlane const& bp,
                            vw::Vector3 const& proj_pt,
                            std::vector<LocalRasterSample>& samples) {
  samples.clear();

  vw::ImageView<vw::PixelMask<float>> const& raster = bp.water_surface;
  int cols = raster.cols(), rows = raster.rows();
  if (cols <= 0 || rows <= 0)
    return false;

  // Convert stereographic proj_pt to raster-native coords for pixel lookup.
  vw::Vector3 ecef_at_query = vw::bathyUnprojPoint(bp.stereographic_proj, proj_pt);
  vw::Vector3 raster_proj_pt = vw::bathyProjPoint(bp.plane_proj, ecef_at_query);

  // Pixel at current point
  vw::Vector2 pix0
    = bp.plane_proj.point_to_pixel(vw::Vector2(raster_proj_pt[0], raster_proj_pt[1]));
  if (!isInBounds(pix0, cols, rows))
    return false;

  // Right or left neighbor
  vw::Vector2 pix1 = pix0 + vw::Vector2(1.0, 0.0);
  if (!isInBounds(pix1, cols, rows))
    pix1 = pix0 + vw::Vector2(-1.0, 0.0);
  if (!isInBounds(pix1, cols, rows))
    return false;

  // Top or bottom neighbor
  vw::Vector2 pix2 = pix0 + vw::Vector2(0.0, 1.0);
  if (!isInBounds(pix2, cols, rows))
    pix2 = pix0 + vw::Vector2(0.0, -1.0);
  if (!isInBounds(pix2, cols, rows))
    return false;

  auto interp = vw::BilinearInterpolation::interpolator(raster);
  vw::PixelMask<float> v0 = interp(raster, pix0[0], pix0[1], 0);
  vw::PixelMask<float> v1 = interp(raster, pix1[0], pix1[1], 0);
  vw::PixelMask<float> v2 = interp(raster, pix2[0], pix2[1], 0);
  if (!vw::is_valid(v0) || !vw::is_valid(v1) || !vw::is_valid(v2))
    return false;

  samples.resize(3);
  samples[0].pix = pix0; samples[0].height = v0.child();
  samples[1].pix = pix1; samples[1].height = v1.child();
  samples[2].pix = pix2; samples[2].height = v2.child();
  return true;
}

// Convert a (pixel, height) pair sampled from bp.water_surface into a
// point in bp.stereographic_proj's meter-scale frame, via ECEF.
vw::Vector3 rasterPixelToStereographicPoint(BathyPlane const& bp,
                                            vw::Vector2 const& pix,
                                            double height) {
  vw::Vector2 raster_xy = bp.plane_proj.pixel_to_point(pix);
  vw::Vector3 ecef = vw::bathyUnprojPoint(bp.plane_proj,
                                          vw::Vector3(raster_xy[0], raster_xy[1], height));
  return vw::bathyProjPoint(bp.stereographic_proj, ecef);
}

// Sample 3 points on the projected plane and fit a local ECEF plane
// approximation. Accurate within ~20 m (flat-Earth limit); used by the
// Newton-Raphson refraction solver to iterate entirely in ECEF and avoid
// expensive per-iteration proj/unproj round-trips.
//
// offset_meters is treated as a Euclidean length in plane_proj's axes,
// so those axes must be meter-scale. Callers should pass
// BathyPlane::stereographic_proj; passing a geographic georef would
// interpret 1.0 as 1 degree (~111 km) and produce a bogus tangent.
std::vector<double>
fitLocalEcefPlaneToProjPlane(std::vector<double> const& plane,
                             vw::cartography::GeoReference const& plane_proj,
                             vw::Vector3 const& proj_pt,
                             double offset_meters) {

  if (!plane_proj.is_projected())
    vw_throw(vw::LogicErr() << "fitLocalEcefPlaneToProjPlane requires a "
             << "projected (meter-scale) georef; got geographic instead.\n");

  // Project proj_pt onto the plane to get a point on the plane.
  double dist = signedDistToPlane(plane, proj_pt);
  vw::Vector3 normal(plane[0], plane[1], plane[2]);
  double normal_sq = dot_prod(normal, normal);
  vw::Vector3 proj_on_plane = proj_pt - (dist / normal_sq) * normal;

  // Sample 3 points on the plane by moving along two orthogonal tangent
  // vectors that lie in the plane.
  vw::Vector3 tangent1, tangent2;
  vw::math::computePlaneTangents(normal, offset_meters, tangent1, tangent2);
  vw::Vector3 pt0 = proj_on_plane;
  vw::Vector3 pt1 = proj_on_plane + tangent1;
  vw::Vector3 pt2 = proj_on_plane + tangent2;

  // Convert the 3 projected points to ECEF.
  vw::Vector3 ecef0 = vw::bathyUnprojPoint(plane_proj, pt0);
  vw::Vector3 ecef1 = vw::bathyUnprojPoint(plane_proj, pt1);
  vw::Vector3 ecef2 = vw::bathyUnprojPoint(plane_proj, pt2);

  // Fit a plane through the 3 ECEF points.
  vw::Vector3 v1 = ecef1 - ecef0;
  vw::Vector3 v2 = ecef2 - ecef0;
  vw::Vector3 normal_ecef = normalize(cross_prod(v1, v2));
  double D_ecef = dot_prod(normal_ecef, ecef0);

  std::vector<double> ecef_plane(4);
  ecef_plane[0] = normal_ecef[0];
  ecef_plane[1] = normal_ecef[1];
  ecef_plane[2] = normal_ecef[2];
  ecef_plane[3] = -D_ecef;
  return ecef_plane;
}

// Pick three raster neighbors near proj_pt, convert each to ECEF, and fit
// the exact plane through them. The neighbors are pixel-aligned (not
// proj-meter offsets) so the three samples always span a known scale
// regardless of the raster's CRS. offset_meters is only forwarded to the
// plane-based fallback. proj_pt is in bp.stereographic_proj's meter-scale
// frame; the helper does the stereographic -> raster CRS hop internally.
// Falls back to the plane-based fit if pickLocalRasterSamples reports any
// kind of sampling failure.
std::vector<double> fitLocalEcefPlaneToProjSurface(BathyPlane const& bp,
                                                   vw::Vector3 const& proj_pt,
                                                   double offset_meters) {
  std::vector<LocalRasterSample> samples;
  if (!pickLocalRasterSamples(bp, proj_pt, samples))
    return fitLocalEcefPlaneToProjPlane(bp.bathy_plane, bp.stereographic_proj,
                                        proj_pt, offset_meters);

  // Convert each (pixel, height) to an ECEF point through the raster's
  // native georef. No reprojection of the height; it stays in meters.
  vw::Vector3 ecefs[3];
  for (int i = 0; i < 3; i++) {
    vw::Vector2 raster_xy = bp.plane_proj.pixel_to_point(samples[i].pix);
    ecefs[i] = vw::bathyUnprojPoint(bp.plane_proj,
                                    vw::Vector3(raster_xy[0], raster_xy[1],
                                                samples[i].height));
  }

  // Exact plane through three ECEF points.
  vw::Vector3 v = ecefs[1] - ecefs[0];
  vw::Vector3 w = ecefs[2] - ecefs[0];
  vw::Vector3 normal_ecef = normalize(cross_prod(v, w));

  // Orient normal away from Earth center.
  if (dot_prod(normal_ecef, ecefs[0]) < 0)
    normal_ecef = -normal_ecef;

  double D_ecef = dot_prod(normal_ecef, ecefs[0]);
  std::vector<double> ecef_plane(4);
  ecef_plane[0] = normal_ecef[0];
  ecef_plane[1] = normal_ecef[1];
  ecef_plane[2] = normal_ecef[2];
  ecef_plane[3] = -D_ecef;
  return ecef_plane;
}

// Fit a local plane in stereographic_proj's meter-scale frame from 3
// raster neighbors near proj_pt. Used by the camera-to-ground path to
// replace the global best-fit plane with a local raster-derived plane
// right around the ray-surface hit. Returns false on any sampling
// failure or if the three samples are degenerate; caller falls back to
// the global plane in that case.
bool refineLocalPlaneFromRaster(BathyPlane const& bp,
                                vw::Vector3 const& proj_pt,
                                std::vector<double>& plane) {
  plane.clear();

  std::vector<LocalRasterSample> samples;
  if (!pickLocalRasterSamples(bp, proj_pt, samples))
    return false;

  vw::Vector3 s0 = rasterPixelToStereographicPoint(bp, samples[0].pix, samples[0].height);
  vw::Vector3 s1 = rasterPixelToStereographicPoint(bp, samples[1].pix, samples[1].height);
  vw::Vector3 s2 = rasterPixelToStereographicPoint(bp, samples[2].pix, samples[2].height);

  // Exact plane through the three points in stereographic coords.
  vw::Vector3 u = s1 - s0;
  vw::Vector3 w = s2 - s0;
  vw::Vector3 normal = cross_prod(u, w);
  double normal_norm = vw::math::norm_2(normal);
  if (normal_norm < 1e-12)
    return false; // degenerate triangle, fall back to global plane
  normal /= normal_norm;
  if (normal[2] < 0) normal = -normal; // orient +z

  plane.resize(4);
  plane[0] = normal[0];
  plane[1] = normal[1];
  plane[2] = normal[2];
  plane[3] = -dot_prod(normal, s0);
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
    intersect_proj_pt = vw::bathyProjPoint(plane_proj, intersect_ecef);
    vw::Vector3 prev_proj_pt = vw::bathyProjPoint(plane_proj, prev_ecef);
    intersect_proj_dir = intersect_proj_pt - prev_proj_pt;
    intersect_proj_dir /= norm_2(intersect_proj_dir);

    // Stop when we are within 0.1 mm of the plane, while along the ray. Going
    // beyond that seems not useful. This is usually reached on second pass.

    if (std::abs(signedDistToPlane(plane, intersect_proj_pt)) < 1e-4)
      break;

    // Intersect the proj ray with the proj plane
    vw::Vector3 refined_intersect_proj_pt;
    if (!rayPlaneIntersect(intersect_proj_pt, intersect_proj_dir, plane,
                           refined_intersect_proj_pt))
      return false;

    // Convert back to ECEF
    intersect_ecef = vw::bathyUnprojPoint(plane_proj, refined_intersect_proj_pt);

    // Put the point back on the ray. Then it may become slightly off the plane.
    intersect_ecef = in_ecef + dot_prod(intersect_ecef - in_ecef, in_dir) * in_dir;
  }

  return true;
}

// Test Snell's law in projected and unprojected coordinates
void testSnellLaw(std::vector<double> const& plane,
                  vw::cartography::GeoReference const& plane_proj,
                  double refraction_index,
                  vw::Vector3 const& out_ecef,
                  vw::Vector3 const& in_ecef_dir, vw::Vector3 const& out_ecef_dir,
                  vw::Vector3 const& out_proj_pt,
                  vw::Vector3 const& in_proj_dir, vw::Vector3 const& out_proj_dir) {

  // 1. In projected coordinates
  vw::Vector3 proj_normal(plane[0], plane[1], plane[2]);
  double sin_in = sin(acos(dot_prod(proj_normal, -in_proj_dir)));
  double sin_out = sin(acos(dot_prod(-proj_normal, out_proj_dir)));

  // 2. In unprojected coordinates
  vw::Vector3 proj_pt_above_normal = out_proj_pt + 1.0 * proj_normal; // go 1 m along the normal
  vw::Vector3 ecef_above_normal = vw::bathyUnprojPoint(plane_proj,
                                                       proj_pt_above_normal);
  vw::Vector3 ecef_normal = ecef_above_normal - out_ecef;
  ecef_normal /= norm_2(ecef_normal); // normalize
  sin_in = sin(acos(dot_prod(ecef_normal, -in_ecef_dir)));
  sin_out = sin(acos(dot_prod(-ecef_normal, out_ecef_dir)));

  // Verify that the incoming ray, outgoing ray, and the
  // normal are in the same plane in projected coordinates

  // 1. In projected coordinates
  vw::Vector3 in_out_normal = vw::math::cross_prod(in_proj_dir, out_proj_dir);
  double plane_error = dot_prod(in_out_normal, proj_normal);

  // 2. In unprojected coordinates
  in_out_normal = vw::math::cross_prod(in_ecef_dir, out_ecef_dir);
  plane_error = dot_prod(in_out_normal, ecef_normal);
}

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
bool curvedSnellLawAgainstPlane(vw::Vector3 const& in_ecef, vw::Vector3 const& in_dir,
                                std::vector<double> const& plane,
                                vw::cartography::GeoReference const& plane_proj,
                                double refraction_index, double mean_height,
                                vw::Vector3 & out_ecef, vw::Vector3 & out_dir) {

  // Intersect the ray with the curved water surface. Return the intersection
  // point in ecef, that point in projected coordinates, and the ray direction
  // at that location in projected coordinates.
  vw::Vector3 intersect_ecef, intersect_proj_pt, intersect_proj_dir;
  if (!rayBathyPlaneIntersect(in_ecef, in_dir, plane, plane_proj, mean_height,
                              intersect_ecef, intersect_proj_pt, intersect_proj_dir))
    return false;

  // Snell's law in projected coordinates
  vw::Vector3 out_proj_pt, out_proj_dir; // in the water
  bool ans = snellLaw(intersect_proj_pt, intersect_proj_dir,
                      plane, refraction_index,
                      out_proj_pt, out_proj_dir);

  // If Snell's law failed to work, exit early
  if (!ans)
    return ans;

  // Move a little on the ray in projected coordinates.
  vw::Vector3 next_proj_pt = out_proj_pt + 1.0 * out_proj_dir;

  // Convert back to ECEF
  out_ecef = vw::bathyUnprojPoint(plane_proj, out_proj_pt);
  vw::Vector3 next_ecef = vw::bathyUnprojPoint(plane_proj, next_proj_pt);

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

} // namespace

// Fit local ECEF plane to surface in projected coordinates or to plane in
// projected coordinates.
std::vector<double> fitLocalEcefPlane(BathyPlane const& bp,
                                      vw::Vector3 const& proj_pt,
                                      double offset_meters) {
  if (bp.water_surface.cols() > 0)
    return fitLocalEcefPlaneToProjSurface(bp, proj_pt, offset_meters);
  return fitLocalEcefPlaneToProjPlane(bp.bathy_plane, bp.stereographic_proj,
                                      proj_pt, offset_meters);
}

// Bend a camera ray at the bathy water surface. See BathyRay.h for the
// contract. When bp carries a raster, this routine first locates the
// approximate hit using the global plane, samples three raster neighbors
// there to fit a local plane, and bends with that refined local plane
// (falling back to the global plane on any sampling failure).
//
// TODO(oalexan1): the refined-plane path does two ray-plane intersects
// (one to locate the raster sampling point, one inside curvedSnellLawAgainstPlane).
// Could be collapsed to one by teaching rayBathyPlaneIntersect to take
// an ECEF seed; defer until profiling justifies it.
bool curvedSnellLaw(vw::Vector3 const& in_ecef,
                    vw::Vector3 const& in_dir,
                    BathyPlane const& bp,
                    double refraction_index,
                    vw::Vector3& out_ecef,
                    vw::Vector3& out_dir) {

  // No raster: classic single-plane path.
  if (bp.water_surface.cols() == 0)
    return curvedSnellLawAgainstPlane(in_ecef, in_dir,
                                      bp.bathy_plane, bp.stereographic_proj,
                                      refraction_index, bp.mean_height,
                                      out_ecef, out_dir);

  // Raster path: locate the approximate hit, then refine the plane there.
  vw::Vector3 hit_ecef, hit_proj_pt, hit_proj_dir;
  if (!rayBathyPlaneIntersect(in_ecef, in_dir,
                              bp.bathy_plane, bp.stereographic_proj,
                              bp.mean_height,
                              hit_ecef, hit_proj_pt, hit_proj_dir))
    return false;

  std::vector<double> local_plane;
  bool refined = refineLocalPlaneFromRaster(bp, hit_proj_pt, local_plane);
  std::vector<double> const& active_plane = refined ? local_plane : bp.bathy_plane;

  return curvedSnellLawAgainstPlane(in_ecef, in_dir,
                                    active_plane, bp.stereographic_proj,
                                    refraction_index, bp.mean_height,
                                    out_ecef, out_dir);
}

} // namespace vw
