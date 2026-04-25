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

// BathyStereoModel - the stereo triangulator that knows how to bend rays
// at a water surface. The bathy data primitives (BathyPlane, plane fits,
// raster sampling) live in BathyData. The pure ray / plane math lives in
// SnellLaw. The camera-aware helpers (point_to_pixel, datumBathyIntersection)
// live in BathyCamera.

#include <vw/Cartography/BathyStereoModel.h>

#include <vw/Cartography/BathyData.h>
#include <vw/Cartography/BathyRay.h>
#include <vw/Cartography/SnellLaw.h>
#include <vw/Camera/CameraModel.h>
#include <vw/Stereo/StereoModel.h>
#include <vw/Math/Vector.h>
#include <vw/Core/Exception.h>

namespace vw {

// Settings used for bathymetry correction
void BathyStereoModel::set_bathy(double refraction_index,
                                 std::vector<BathyPlane> const& bathy_plane_vec) {

  m_bathy_correct = true;
  m_refraction_index = refraction_index;
  m_bathy_plane_vec = bathy_plane_vec;

  if (m_refraction_index <= 1)
    vw::vw_throw(vw::ArgumentErr() << "The water refraction index must be bigger than 1.");

  if (m_bathy_plane_vec.size() != 2)
    vw::vw_throw(vw::ArgumentErr() << "Expecting two bathy planes (left and right).");

  for (int it = 0; it < 2; it++) {
    if (m_bathy_plane_vec[it].bathy_plane.size() != 4)
      vw::vw_throw(vw::ArgumentErr() << "The bathy plane must have 4 coefficients.");
  }

  // The default behavior is for the left and right bathy planes to be the same.
  // Yet we allow them to be different. Here need to check. Compare the
  // stereographic frame since that is where the plane coefficients live.
  m_single_bathy_plane = true;
  if (m_bathy_plane_vec[0].stereographic_proj.proj4_str()
      != m_bathy_plane_vec[1].stereographic_proj.proj4_str())
      m_single_bathy_plane = false;
  if (m_bathy_plane_vec[0].bathy_plane != m_bathy_plane_vec[1].bathy_plane)
    m_single_bathy_plane = false;
}

// Compute the rays intersection. Note that even if we are in
// bathymetry mode, so m_bathy_correct is true, for this particular
// pair of rays we may have do_bathy false, and then we won't do the
// correction.  Return also a flag saying if we did bathymetry
// correction or not. When the rays intersect above the water surface,
// the correction is not done, but a valid 3D point is still returned.
Vector3 BathyStereoModel::operator()(std::vector<Vector2> const& pixVec,
                                     Vector3& errorVec, bool do_bathy,
                                     bool & did_bathy) const {

  // Initialize the outputs
  did_bathy = false;
  errorVec = Vector3();

  int num_cams = m_cameras.size();
  VW_ASSERT((int)pixVec.size() == num_cams,
            vw::ArgumentErr() << "the number of rays must match "
            << "the number of cameras.\n");

  try {

    std::vector<Vector3> camDirs(num_cams), camCtrs(num_cams);
    camDirs.clear(); camCtrs.clear();

    // Pick the valid rays
    for (int p = 0; p < num_cams; p++) {

      Vector2 pix = pixVec[p];
      if (pix != pix || // i.e., NaN
          pix == camera::CameraModel::invalid_pixel())
        continue;

      camDirs.push_back(m_cameras[p]->pixel_to_vector(pix));
      camCtrs.push_back(m_cameras[p]->camera_center(pix));
    }

    // Not enough valid rays
    if (camDirs.size() < 2)
      return Vector3();

    if (are_nearly_parallel(m_angle_tol, camDirs))
      return Vector3();

    // Determine range by triangulation
    Vector3 uncorr_tri_pt = triangulate_point(camDirs, camCtrs, errorVec);

    // Reflect points that fall behind one of the two cameras.  Do
    // not do this when bathymetry mode is on, as then we surely
    // have satellite images and there is no way a point would be
    // behind the camera.
    if (!m_bathy_correct) {
      bool reflect = false;
      for (int p = 0; p < (int)camCtrs.size(); p++)
        if (dot_prod(uncorr_tri_pt - camCtrs[p], camDirs[p]) < 0)
          reflect = true;
      if (reflect)
        uncorr_tri_pt = -uncorr_tri_pt + 2*camCtrs[0];
    }

    if (!do_bathy || camDirs.size() != 2)
      return uncorr_tri_pt;

    // Continue with bathymetry correction

    if (!m_bathy_correct)
      vw::vw_throw(vw::ArgumentErr()
                    << "Requested to do bathymetry correction while "
                    << "this mode was not set up.");

    // Find the rays after bending, according to Snell's law.
    std::vector<Vector3> waterDirs(2), waterCtrs(2);

    // When there's a single plane, things are simple.
    // Rays get bent, then they intersect, and done.
    if (m_single_bathy_plane) {

      // The water surface is curved. It is however flat (a plane) when we
      // switch to the stereographic frame (where bathy_plane coefs live).
      Vector3 proj_pt = bathyProjPoint(m_bathy_plane_vec[0].stereographic_proj,
                                       uncorr_tri_pt);
      double ht_val = signed_dist_to_plane(m_bathy_plane_vec[0].bathy_plane, proj_pt);
      if (ht_val >= 0) {
        // the rays intersect above the water surface
        did_bathy = false;
        return uncorr_tri_pt;
      }

      for (size_t it = 0; it < 2; it++) {
        // Bend each ray at the surface. Raster-refined plane when available.
        bool ans = curvedSnellLaw(camCtrs[it], camDirs[it],
                                  m_bathy_plane_vec[it],
                                  m_refraction_index,
                                  waterCtrs[it], waterDirs[it]);
        if (!ans) {
          did_bathy = false;
          return uncorr_tri_pt;
        }
      }

      // Re-triangulate with the new rays, after they get bent
      Vector3 corr_tri_pt = triangulate_point(waterDirs, waterCtrs, errorVec);

      did_bathy = true;
      return corr_tri_pt;
    }

    // The case of left and right images having their own bathy planes

    // Bend the rays
    for (size_t it = 0; it < 2; it++) {
      // Bend each ray at the surface. Raster-refined plane when available.
      bool ans = curvedSnellLaw(camCtrs[it], camDirs[it],
                                m_bathy_plane_vec[it],
                                m_refraction_index,
                                waterCtrs[it], waterDirs[it]);
      if (!ans)
        return uncorr_tri_pt;
    }

    // Each ray has two parts: before bending and after it. Two
    // bent rays can intersect on their unbent parts, the bent part
    // of one ray with unbent part of another ray, unbent part of
    // one ray with bent part of another ray, and bent parts of both
    // rays. Handle all these with much care.

    Vector3 err, tri_pt;
    std::vector<double> signed_dists;

    // See if the unbent portions intersect above their planes
    tri_pt = vw::stereo::triangulate_pair(camDirs[0], camCtrs[0], camDirs[1], camCtrs[1], err);
    signed_distances_to_planes(m_bathy_plane_vec, tri_pt, signed_dists);
    if (signed_dists[0] >= 0 && signed_dists[1] >= 0) {
      did_bathy = false; // since the rays did not reach the bathy plane
      errorVec = err;
      return tri_pt;
    }

    // See if the bent portions intersect below their planes
    tri_pt = vw::stereo::triangulate_pair(waterDirs[0], waterCtrs[0],
                                          waterDirs[1], waterCtrs[1], err);
    signed_distances_to_planes(m_bathy_plane_vec, tri_pt, signed_dists);
    if (signed_dists[0] <= 0 && signed_dists[1] <= 0) {
      did_bathy = true; // the resulting point is at least under one plane
      errorVec = err;
      return tri_pt;
    }

    // See if the left unbent portion intersects the right bent portion,
    // above left's water plane and below right's water plane
    tri_pt = vw::stereo::triangulate_pair(camDirs[0], camCtrs[0], waterDirs[1],
                                          waterCtrs[1], err);
    signed_distances_to_planes(m_bathy_plane_vec, tri_pt, signed_dists);
    if (signed_dists[0] >= 0 && signed_dists[1] <= 0) {
      did_bathy = true; // the resulting point is at least under one plane
      errorVec = err;
      return tri_pt;
    }

    // See if the left bent portion intersects the right unbent portion,
    // below left's water plane and above right's water plane
    tri_pt = vw::stereo::triangulate_pair(waterDirs[0], waterCtrs[0],
                                          camDirs[1], camCtrs[1], err);
    signed_distances_to_planes(m_bathy_plane_vec, tri_pt, signed_dists);
    if (signed_dists[0] <= 0 && signed_dists[1] >= 0) {
      did_bathy = true; // the resulting point is at least under one plane
      errorVec = err;
      return tri_pt;
    }

  } catch (const camera::PixelToRayErr& /*e*/) {}

  // We arrive here only when there's bad luck
  did_bathy = false;
  errorVec = vw::Vector3();
  return vw::Vector3();
}

Vector3 BathyStereoModel::operator()(std::vector<Vector2> const& pixVec,
                                      double& error) const {
  vw::vw_throw(vw::NoImplErr() << "Not implemented for BathyStereoModel.");
  return Vector3();
}

Vector3 BathyStereoModel::operator()(Vector2 const& pix1,
                                      Vector2 const& pix2, Vector3& errorVec) const {
  vw::vw_throw(vw::NoImplErr() << "Not implemented for BathyStereoModel.");
  return Vector3();
}

Vector3 BathyStereoModel::operator()(Vector2 const& pix1, Vector2 const& pix2,
                                      double& error) const {
  vw::vw_throw(vw::NoImplErr() << "Not implemented for BathyStereoModel.");
  return Vector3();
}

} // namespace vw
