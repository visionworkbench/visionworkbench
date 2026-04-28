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

// Camera-aware bathymetry helpers. See BathyCamera.h for the contract.

#include <vw/Cartography/BathyCamera.h>

#include <vw/Cartography/BathyData.h>
#include <vw/Cartography/BathyRay.h>
#include <vw/Cartography/SnellLaw.h>
#include <vw/Camera/CameraModel.h>
#include <vw/Math/Vector.h>
#include <vw/Math/VectorUtils.h>
#include <vw/Math/NewtonRaphson.h>
#include <vw/Core/Exception.h>

namespace vw {

namespace {

// A class that helps with finding the projection of an ECEF point through a
// bathy plane in the camera. This plane is in ECEF and is tangent to the true
// bathy plane around the ECEF point. This allows all calculations to be in
// ECEF, which is about 100x faster than having to deal with projection
// operations. Even so this is 10 slower than not having a bathy plane.
class BathyFunctorEcef {
public:
  BathyFunctorEcef(vw::Vector3 const& ecef_point,
                   vw::camera::CameraModel const* cam,
                   std::vector<double> const& ecef_plane,
                   double refraction_index):
  m_cam(cam), m_origin(ecef_point),
  m_ecef_plane(ecef_plane),
  m_refraction_index(refraction_index) {

    // Form local coordinate system using the ECEF point
    vw::math::formBasis(ecef_point, m_x_axis, m_y_axis, m_normal);
  }

  // Intersect a ray with the tangent plane and return 2D coordinates in the tangent
  // plane coordinate system.
  vw::Vector2 rayPlaneIntersect(vw::Vector3 const& ray_pt, vw::Vector3 const& ray_dir) const {

    double denom = vw::math::dot_prod(ray_dir, m_normal);
    if (std::abs(denom) < 1e-10)
      vw::vw_throw(vw::ArgumentErr() << "Ray is parallel to tangent plane.\n");

    double t = vw::math::dot_prod(m_origin - ray_pt, m_normal) / denom;
    vw::Vector3 intersection = ray_pt + t * ray_dir;

    // Project intersection onto tangent plane axes
    vw::Vector3 offset = intersection - m_origin;
    double x = vw::math::dot_prod(offset, m_x_axis);
    double y = vw::math::dot_prod(offset, m_y_axis);

    return vw::Vector2(x, y);
  }

  // Operator for use with Newton-Raphson solver
  // Input: pix is a pixel in the camera
  // Output: 2D coordinates in the tangent plane after ray bending
  vw::Vector2 operator()(vw::Vector2 const& pix) const {
    // Get ray from camera (without bathy correction)
    vw::Vector3 cam_ctr = m_cam->camera_center(pix);
    vw::Vector3 cam_dir = m_cam->pixel_to_vector(pix);

    // Apply Snell's law to bend the ray at the water surface. Everything stays
    // in ECEF.
    vw::Vector3 refracted_pt, refracted_dir;
    bool success = vw::snellLaw(cam_ctr, cam_dir,
                                m_ecef_plane,
                                m_refraction_index,
                                refracted_pt, refracted_dir);

    if (!success)
      vw::vw_throw(vw::ArgumentErr() << "Snell's law refraction failed.\n");

    // Intersect refracted ray with tangent plane and return 2D coordinates
    return rayPlaneIntersect(refracted_pt, refracted_dir);
  }

  vw::camera::CameraModel const* m_cam; // Camera pointer
  vw::Vector3 m_origin;                 // Origin of tangent plane (ECEF point P)
  std::vector<double> m_ecef_plane;     // Local ECEF plane coefficients [A,B,C,-D]
  double m_refraction_index;            // Index of refraction
  vw::Vector3 m_x_axis;                 // East direction (tangent plane X axis)
  vw::Vector3 m_y_axis;                 // North direction (tangent plane Y axis)
  vw::Vector3 m_normal;                 // Up direction (perpendicular to plane)
};

} // namespace

// Intersect a ray from camera center along camera direction with the datum at
// given semi-axes, with optional bathymetry correction. If the ray passes
// through the bathy plane (water surface) before it meets the datum, apply
// Snell's law refraction and continue with the new bent ray until reaching the
// datum. Returns the intersection point, or zero vector on failure.
Vector3 datumBathyIntersection(Vector3 const& cam_ctr,
                               Vector3 const& cam_dir,
                               double major_axis, double minor_axis,
                               BathyPlane const& bathy_plane,
                               double refraction_index) {

  // First, intersect ray with datum
  Vector3 xyz = vw::cartography::datum_intersection(major_axis, minor_axis,
                                                    cam_ctr, cam_dir);

  // If intersection failed, return zero vector
  if (xyz == Vector3(0, 0, 0))
    return Vector3(0, 0, 0);

  // Check signed distance to the water surface (raster-aware if a wl.tif
  // is loaded; falls back to the best-fit plane otherwise).
  double ht_val = signedDistToPlane(bathy_plane, xyz);

  // If point is above water surface, no refraction needed
  if (ht_val >= 0)
    return xyz;

  // Point is below water - need to apply Snell's law refraction. This
  // uses a bathy plane given by a formula or as an raster of values.
  Vector3 out_xyz, out_dir;
  bool success = curvedSnellLaw(cam_ctr, cam_dir,
                                bathy_plane,
                                refraction_index,
                                out_xyz, out_dir);

  // If Snell's law failed, return zero vector
  if (!success)
    return Vector3(0, 0, 0);

  // Continue refracted ray to datum
  Vector3 refracted_xyz = vw::cartography::datum_intersection(major_axis, minor_axis,
                                                              out_xyz, out_dir);

  // If refracted intersection failed, return zero vector
  if (refracted_xyz == Vector3(0, 0, 0))
    return Vector3(0, 0, 0);

  return refracted_xyz;
}

// Project an ECEF point to camera pixel, accounting for bathymetry if the point
// lies below the water surface. Fits a local ECEF tangent plane. Good for
// shallow-water bathymetry so the point is at most meters deep. Degrades for
// deep points and steep off-nadir rays.
vw::Vector2 point_to_pixel(vw::camera::CameraModel const* cam,
                           vw::BathyPlane const& bathy_plane,
                           double refraction_index,
                           vw::Vector3 const& ecef_point) {

  // Get camera pixel as if there was no refraction
  vw::Vector2 pix = cam->point_to_pixel(ecef_point);

  // Check signed distance to the water surface (raster-aware if a wl.tif
  // is loaded; falls back to the best-fit plane otherwise).
  double dist = vw::signedDistToPlane(bathy_plane, ecef_point);

  // If point is above the water surface (positive distance), no refraction
  if (dist >= 0)
    return pix;

  // Project to the stereographic frame for fitLocalEcefPlane below.
  vw::Vector3 proj_pt = vw::bathyProjPoint(bathy_plane.stereographic_proj, ecef_point);

  // Fit local ECEF tangent plane. Newton-Raphson below iterates against
  // this tangent entirely in ECEF - zero per-iteration proj calls. When a
  // raster water surface is present, the tangent reflects the raster near
  // this query; otherwise it reflects the global fitted plane.
  double offset = 1.0;
  std::vector<double> ecef_plane = fitLocalEcefPlane(bathy_plane, proj_pt, offset);

  // Point is below water surface - need to account for refraction
  // Set up the BathyFunctorEcef for Newton-Raphson iteration (uses local ECEF plane)
  // Old version: BathyFunctor bathy_func(ecef_point, cam, bathy_plane, refraction_index);
  BathyFunctorEcef bathy_func(ecef_point, cam, ecef_plane, refraction_index);

  // Set up Newton-Raphson solver with numerical Jacobian
  vw::math::NewtonRaphson nr(bathy_func);

  // TODO(oalexan1): Study what ecef point is between takes in bundle_adjust,
  // so how much ceres moves it for numerical differences. This wil inform the
  // step size and tolerance to use here.

  // Solve: find pixel such that bathy_func(pix) projects to (0, 0) in tangent plane
  // since ecef_point is the origin of the tangent plane
  vw::Vector2 target(0, 0);
  // Use 0.5 pixel step for numerical differentiation to get Jacobian.
  // We want to avoid issues with small steps here, especially that the input
  // are known to within 1e-4 m or so.
  double step = 0.5;
  double tol = 1e-5; // 1e-5 meter tolerance, should be enough.
  pix = nr.solve(pix, target, step, tol);

  return pix;
}

} // namespace vw
