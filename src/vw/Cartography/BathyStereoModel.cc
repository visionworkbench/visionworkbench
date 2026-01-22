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
#include <vw/Math/LevenbergMarquardt.h>
#include <vw/Math/NewtonRaphson.h>
#include <vw/Camera/CameraModel.h>
#include <vw/Stereo/StereoModel.h>
#include <vw/Cartography/BathyStereoModel.h>
#include <vw/Core/Exception.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/FileIO/DiskImageUtils.h>
#include <vw/Image/MaskViews.h>

#include <iostream>
#include <fstream>
#include <sstream>

namespace vw {

vw::ImageViewRef<vw::PixelMask<float>> read_bathy_mask(std::string const& filename,
                                                       float & nodata_val) {
  float local_nodata = -std::numeric_limits<float>::max();
  if (!vw::read_nodata_val(filename, local_nodata))
    vw::vw_throw(vw::ArgumentErr() << "Unable to read the nodata value from "
             << filename);
  nodata_val = local_nodata;
  return vw::create_mask(vw::DiskImageView<float>(filename), local_nodata);
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
  vw_out() << "Read projection: " <<  plane_proj.overall_proj4_str()
            << "\n";
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

// Test Snell's law in projected and unprojected coordinates
void testSnellLaw(std::vector<double> const& plane,
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
                            vw::Vector3 & intersect_ecef,
                            vw::Vector3 & intersect_proj_pt,
                            vw::Vector3 & intersect_proj_dir) {
  
  // Find the mean water surface
  double mean_ht = -plane[3] / plane[2];
  double major_radius = plane_proj.datum().semi_major_axis() + mean_ht;
  double minor_radius = plane_proj.datum().semi_minor_axis() + mean_ht;

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
                    double refraction_index,
                    Vector3 & out_ecef, Vector3 & out_dir) {

  // Intersect the ray with the curved water surface. Return the intersection
  // point in ecef, that point in projected coordinates, and the ray direction
  // at that location in projected coordinates.
  Vector3 intersect_ecef, intersect_proj_pt, intersect_proj_dir;
  if (!rayBathyPlaneIntersect(in_ecef, in_dir, plane, plane_proj,
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

// A class to manage a local tangent plane at a given ECEF point. The tangent
// plane is tangent to datum at that point. For an ellipsoid datum, this plane
// is perpendicular to the geodetic normal (not the geocentric radius vector).
// The X and Y axes are East and North respectively.
class BathyFunctor {
public:
  BathyFunctor(vw::Vector3 const& ecef_point,
               vw::camera::CameraModel const* cam,
               vw::BathyPlane const& bathy_plane,
               double refraction_index): m_cam(cam), m_origin(ecef_point),
                                         m_bathy_plane(bathy_plane),
                                         m_refraction_index(refraction_index) {
    // Form local coordinate system using the ECEF point
    vw::math::formBasis(ecef_point, m_x_axis, m_y_axis, m_normal);
  }

  // Intersect a ray with the tangent plane and return 2D coordinates in the tangent
  // plane coordinate system.
  vw::Vector2 rayPlaneIntersect(vw::Vector3 const& ray_pt,
                                 vw::Vector3 const& ray_dir) const {

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

    // Apply Snell's law to bend the ray at the water surface
    vw::Vector3 refracted_pt, refracted_dir;
    bool success = vw::curvedSnellLaw(cam_ctr, cam_dir,
                                      m_bathy_plane.bathy_plane,
                                      m_bathy_plane.plane_proj,
                                      m_refraction_index,
                                      refracted_pt, refracted_dir);

    if (!success)
      vw::vw_throw(vw::ArgumentErr() << "Snell's law refraction failed.\n");

    // Intersect refracted ray with tangent plane and return 2D coordinates
    return rayPlaneIntersect(refracted_pt, refracted_dir);
  }

  vw::camera::CameraModel const* m_cam; // Camera pointer
  vw::Vector3 m_origin;             // Origin of tangent plane (ECEF point P)
  vw::BathyPlane m_bathy_plane;     // Bathy plane for refraction
  double m_refraction_index;        // Index of refraction
  vw::Vector3 m_x_axis;             // East direction (tangent plane X axis)
  vw::Vector3 m_y_axis;             // North direction (tangent plane Y axis)
  vw::Vector3 m_normal;             // Up direction (perpendicular to plane)
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

  // Point is below water surface - need to account for refraction
  // Set up the BathyFunctor for Newton-Raphson iteration
  BathyFunctor bathy_func(ecef_point, cam, bathy_plane, refraction_index);

  // Set up Newton-Raphson solver with numerical Jacobian
  vw::math::NewtonRaphson nr(bathy_func);

  // Solve: find pixel such that bathy_func(pix) projects to (0, 0) in tangent plane
  // since ecef_point is the origin of the tangent plane
  vw::Vector2 target(0, 0);
  double step = 0.5;    // 0.5 meter step for numerical differentiation
  double tol = 1e-6;    // 1e-6 pixels tolerance

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
