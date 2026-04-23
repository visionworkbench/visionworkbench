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

#include <vw/Image/ImageView.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/PixelMath.h>
#include <vw/Image/PixelTypeInfo.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Math/Vector.h>
#include <vw/Math/VectorUtils.h>
#include <vw/Math/NewtonRaphson.h>
#include <vw/Camera/CameraModel.h>
#include <vw/Stereo/StereoModel.h>
#include <vw/Cartography/BathyStereoModel.h>
#include <vw/Cartography/SnellLaw.h>
#include <vw/Core/Exception.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/FileIO/DiskImageUtils.h>
#include <vw/Image/MaskViews.h>

#include <iostream>
#include <fstream>
#include <sstream>

namespace vw {

// Functor to invalidate pixels with non-positive values (water pixels)
struct InvalidateNonPositive:
  public vw::ReturnFixedType<vw::PixelMask<float>> {
  vw::PixelMask<float> operator()(vw::PixelMask<float> p) const {
    if (is_valid(p) && p.child() <= 0.0f)
      p.invalidate();
    return p;
  }
};

// Read a bathy mask. Water pixels are those with non-positive values or matching
// the file's nodata value. Both are invalidated in the returned masked image.
// The returned nodata_val is suitable for writing the mask back.
vw::ImageViewRef<vw::PixelMask<float>> read_bathy_mask(std::string const& filename,
                                                       float & nodata_val) {
  float local_nodata = -std::numeric_limits<float>::max();
  bool has_nodata = vw::read_nodata_val(filename, local_nodata);

  // Set a sensible output nodata value for callers that write the mask back
  if (has_nodata)
    nodata_val = std::max(0.0f, local_nodata);
  else
    nodata_val = 0.0f;

  // Create mask using the file's nodata value (if absent, local_nodata is
  // -max_float so nothing matches and nothing is invalidated at this step)
  auto masked = vw::create_mask(vw::DiskImageView<float>(filename), local_nodata);

  // Additionally invalidate pixels with non-positive values (water pixels).
  // This ensures that masks with 0 = water, 1 = land work correctly,
  // regardless of whether the file has a nodata tag.
  return vw::per_pixel_view(masked, InvalidateNonPositive());
}

void read_bathy_masks(std::vector<std::string> const& mask_filenames,
                      std::vector<vw::ImageViewRef<vw::PixelMask<float>>> & bathy_masks) {
  bathy_masks.clear();
  for (size_t i = 0; i < mask_filenames.size(); i++) {
    float nodata_val = -std::numeric_limits<float>::max(); // part of API
    bathy_masks.push_back(read_bathy_mask(mask_filenames[i], nodata_val));
  }
}

// Check if the given left and right pixels are in the masked region (invalid in
// the mask). That will mean bathymetry correction should be applied.
bool areMasked(ImageViewRef<PixelMask<float>> const& left_mask,
               ImageViewRef<PixelMask<float>> const& right_mask,
               Vector2 const& lpix, Vector2 const& rpix) {

  Vector2i ilpix(round(lpix.x()), round(lpix.y()));
  Vector2i irpix(round(rpix.x()), round(rpix.y()));

  if (ilpix.x() < 0 || ilpix.x() >= left_mask.cols() ||
      ilpix.y() < 0 || ilpix.y() >= left_mask.rows()) return false;

  if (irpix.x() < 0 || irpix.x() >= right_mask.cols() ||
      irpix.y() < 0 || irpix.y() >= right_mask.rows()) return false;

  return (!is_valid(left_mask(ilpix.x(), ilpix.y())) &&
          !is_valid(right_mask(irpix.x(), irpix.y())));
}

// Read the bathy plane. It must be in some local stereographic projection which
// is read as well. The file format is: plane coefficients on first line,
// comment line, then projection lat/lon.
void readBathyPlane(std::string const& bathy_plane_file,
                    std::vector<double> & bathy_plane,
                    vw::cartography::GeoReference & plane_proj) {

  vw::vw_out() << "Reading bathy plane: " << bathy_plane_file << "\n";
  std::ifstream handle;
  handle.open(bathy_plane_file.c_str());
  if (handle.fail())
    vw_throw(vw::IOErr() << "Unable to open bathy plane file: " << bathy_plane_file << "\n");

  std::vector<std::string> lines;
  std::string line;
  while (getline(handle, line, '\n'))
    lines.push_back(line);

  if (lines.size() < 3)
    vw_throw(vw::IOErr() << "Bathy plane file must have at least 3 lines: "
              << "plane coefficients, comment, and projection lat/lon.\n");

  // Read the four plane coefficients
  bathy_plane.resize(4);
  {
    std::istringstream iss(lines[0]);
    if (!(iss >> bathy_plane[0] >> bathy_plane[1] >> bathy_plane[2] >> bathy_plane[3]))
      vw_throw(vw::IOErr() << "Could not read four values from first "
                << "line of the bathy plane.\n");
  }

  // The second line must be a comment
  if (lines[1][0] != '#')
    vw_throw(vw::IOErr() << "Second line of bathy plane file must be a comment.\n");

  if (bathy_plane[2] <= 0)
    vw_throw(vw::IOErr() << "For a curved water surface, the third value "
              << "must be positive.\n");

  // Read projection center
  double scale = 1.0;
  double proj_lat, proj_lon;
  {
    std::istringstream iss(lines[2]);
    if (!(iss >> proj_lat >> proj_lon))
      vw_throw(vw::IOErr() << "Could not read the projection latitude and longitude.\n");
  }
  vw::cartography::Datum datum("WGS_1984");
  plane_proj.set_datum(datum);
  plane_proj.set_stereographic(proj_lat, proj_lon, scale);
  vw_out() << "Read projection: " <<  plane_proj.overall_proj4_str() << "\n";
}

// Read the bathy planes and associated data. More often than not they will be
// identical. If there is more than one bathy plane file, they are all kept in
// the same string, separated by space.
void readBathyPlanes(std::string const& bathy_plane_files,
                     int num_images,
                     std::vector<BathyPlane> & bathy_plane_vec) {

  bathy_plane_vec.clear();

  std::string bathy_plane_file;
  std::istringstream iss(bathy_plane_files);
  while (iss >> bathy_plane_file) {
    bathy_plane_vec.push_back(BathyPlane());
    readBathyPlane(bathy_plane_file,
                   // Outputs
                   bathy_plane_vec.back().bathy_plane,
                   bathy_plane_vec.back().plane_proj);
  }

  if (bathy_plane_vec.size() != 1 && bathy_plane_vec.size() != (size_t)num_images)
    vw_throw(vw::ArgumentErr() << "1 or " << num_images << " bathy planes expected.\n");

  // Clone the bathy plane if there's only one
  while (bathy_plane_vec.size() < (size_t)num_images)
    bathy_plane_vec.push_back(bathy_plane_vec[0]);
}

// Given a ECEF point xyz, and two planes, find if xyz is above or below each of
// the plane by finding the signed distances to them. The water surface is curved
// (modeled as a plane in local stereographic projection coordinates).
void signed_distances_to_planes(std::vector<BathyPlane> const& bathy_plane_vec,
                                vw::Vector3 const& xyz,
                                std::vector<double> & distances) {

  if (bathy_plane_vec.size() != 2)
    vw_throw(vw::ArgumentErr() << "Two bathy planes expected.\n");

  distances.resize(2);
  for (size_t it = 0; it < 2; it++) {
    // Convert xyz to projected coordinates, then compute signed distance
    distances[it]
      = signed_dist_to_plane(bathy_plane_vec[it].bathy_plane,
                             proj_point(bathy_plane_vec[it].plane_proj, xyz));
  }
}

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

  // Project the intersection point to local coordinates
  Vector3 proj_pt = proj_point(bathy_plane.plane_proj, xyz);

  // Check signed distance to bathy plane
  double ht_val = signed_dist_to_plane(bathy_plane.bathy_plane, proj_pt);

  // If point is above water surface, no refraction needed
  if (ht_val >= 0)
    return xyz;

  // Point is below water - need to apply Snell's law refraction
  Vector3 out_xyz, out_dir;
  bool success = curvedSnellLaw(cam_ctr, cam_dir,
                                bathy_plane.bathy_plane,
                                bathy_plane.plane_proj,
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

// Project an ECEF point to pixel, accounting for bathymetry if the point
// is below the bathy plane (water surface).
vw::Vector2 point_to_pixel(vw::camera::CameraModel const* cam,
                           vw::BathyPlane const& bathy_plane,
                           double refraction_index,
                           vw::Vector3 const& ecef_point) {

  // Get camera pixel as if there was no refraction
  vw::Vector2 pix = cam->point_to_pixel(ecef_point);

  // Convert ECEF point to projected coordinates for distance check
  vw::Vector3 proj_pt = vw::proj_point(bathy_plane.plane_proj, ecef_point);

  // Check signed distance to bathy plane
  double dist = vw::signed_dist_to_plane(bathy_plane.bathy_plane, proj_pt);

  // If point is above the water surface (positive distance), no refraction
  if (dist >= 0)
    return pix;

  // Fit local ECEF plane approximation to eliminate expensive proj transforms
  double offset = 1.0;
  std::vector<double> ecef_plane = fitLocalEcefPlane(bathy_plane.bathy_plane,
                                                     bathy_plane.plane_proj,
                                                     proj_pt, dist, offset);

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
  // Yet we allow them to be different. Here need to check.
  m_single_bathy_plane = true;
  if (m_bathy_plane_vec[0].plane_proj.proj4_str()
      != m_bathy_plane_vec[1].plane_proj.proj4_str())
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

      // The water surface is curved. It is however flat (a plane) if we
      // switch to proj coordinates.
      Vector3 proj_pt = proj_point(m_bathy_plane_vec[0].plane_proj,
                                   uncorr_tri_pt);
      double ht_val = signed_dist_to_plane(m_bathy_plane_vec[0].bathy_plane, proj_pt);
      if (ht_val >= 0) {
        // the rays intersect above the water surface
        did_bathy = false;
        return uncorr_tri_pt;
      }

      for (size_t it = 0; it < 2; it++) {
        // Bend each ray at the surface according to Snell's law.
        bool ans = curvedSnellLaw(camCtrs[it], camDirs[it],
                                  m_bathy_plane_vec[it].bathy_plane,
                                  m_bathy_plane_vec[it].plane_proj,
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
      // Bend each ray at the surface according to Snell's law.
      bool ans = curvedSnellLaw(camCtrs[it], camDirs[it],
                                m_bathy_plane_vec[it].bathy_plane,
                                m_bathy_plane_vec[it].plane_proj,
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
