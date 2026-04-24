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

// Implementations of the bathy data free functions declared in
// vw/Cartography/BathyData.h: areMasked, read_bathy_mask(s), and
// readBathyPlanes. These are I/O and predicates on the bathy data
// types; the camera-aware refraction operations (point_to_pixel,
// datumBathyIntersection) and the BathyStereoModel triangulator live
// in BathyStereoModel.cc.

#include <vw/Cartography/BathyData.h>

#include <vw/Image/ImageView.h>
#include <vw/Image/Interpolation.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/PixelMath.h>
#include <vw/Image/PixelTypeInfo.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/MaskViews.h>
#include <vw/Math/Vector.h>
#include <vw/Math/VectorUtils.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/LinearAlgebra.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/SnellLaw.h>
#include <vw/Core/Exception.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/FileIO/DiskImageUtils.h>

#include <cmath>
#include <fstream>
#include <sstream>
#include <limits>

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

} // namespace

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

// Read a plane from the legacy 3-line text format: four plane coefficients
// on line 1, comment on line 2, projection center lat/lon on line 3.
void readBathyPlaneFromText(std::string const& bathy_plane_file, BathyPlane & bp) {
  std::vector<double> & bathy_plane = bp.bathy_plane;
  vw::cartography::GeoReference & plane_proj = bp.plane_proj;

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

  // Text input is already stereographic by construction; stereographic_proj
  // is just a copy of plane_proj.
  bp.stereographic_proj = plane_proj;

  // Guard the invariant: downstream plane math (signed_dist_to_plane,
  // Snell's law in proj coords, tangent-fallback offsets) assumes
  // stereographic_proj has meter-scale axes.
  if (!bp.stereographic_proj.is_projected())
    vw_throw(vw::LogicErr() << "BathyPlane::stereographic_proj must be a "
              << "projected (meter-scale) CRS; got geographic instead.\n");

  // The text format's plane_proj is stereographic meters by construction,
  // so -d/c is the mean water-surface height in meters above the datum.
  bp.mean_height = -bathy_plane[3] / bathy_plane[2];

  vw_out() << "Read projection: " <<  plane_proj.overall_proj4_str() << "\n";
}

// Best-fit plane through a set of 3D points. Returns the four coefficients
// (a, b, c, d) such that a*x + b*y + c*z + d = 0, with (a, b, c) unit length
// and oriented so c > 0. No outlier rejection - use on clean samples.
//
// Formulation: ordinary least squares treating z as the dependent variable,
// minimizing Sum (z_i - (alpha*x_i + beta*y_i + gamma))^2 via the 3x3 normal
// equations. This is insensitive to the relative numerical scale of (x, y)
// vs z, so the fit works across input coordinate systems - including
// geographic (lon-lat-h), where the lon/lat numerical span can be smaller
// than the h span and a centered-points SVD would mistakenly pick a nearly-
// horizontal direction as the plane normal. Least-squares-on-height treats
// height as a function of the horizontal axes by construction.
void fitPlaneToPoints(std::vector<vw::Vector3> const& points,
                      std::vector<double> & bathy_plane) {

  if (points.size() < 3)
    vw_throw(vw::ArgumentErr() << "fitPlaneToPoints: need at least 3 points.\n");

  // Accumulate sums for the normal equations
  //   M * [alpha, beta, gamma]^T = rhs
  // where M is the 3x3 normal-equations matrix (symmetric positive definite
  // for any non-collinear point set).
  double Sx = 0, Sy = 0, Sz = 0;
  double Sxx = 0, Syy = 0, Sxy = 0;
  double Sxz = 0, Syz = 0;
  for (auto const& p : points) {
    Sx += p[0]; Sy += p[1]; Sz += p[2];
    Sxx += p[0]*p[0];
    Syy += p[1]*p[1];
    Sxy += p[0]*p[1];
    Sxz += p[0]*p[2];
    Syz += p[1]*p[2];
  }

  vw::Matrix<double> M(3, 3);
  M(0, 0) = Sxx; M(0, 1) = Sxy; M(0, 2) = Sx;
  M(1, 0) = Sxy; M(1, 1) = Syy; M(1, 2) = Sy;
  M(2, 0) = Sx;  M(2, 1) = Sy;  M(2, 2) = double(points.size());

  vw::Vector<double> rhs(3);
  rhs[0] = Sxz; rhs[1] = Syz; rhs[2] = Sz;

  vw::Vector<double> coef = vw::math::solve_symmetric(M, rhs);
  double alpha = coef[0], beta = coef[1], gamma = coef[2];

  // The fitted equation z = alpha*x + beta*y + gamma rewrites as
  //   alpha*x + beta*y + (-1)*z + gamma = 0.
  // Flip signs to make c = +1, then normalize to a unit normal.
  double a = -alpha, b = -beta, c = 1.0, d = -gamma;
  double norm = std::sqrt(a*a + b*b + c*c);
  a /= norm; b /= norm; c /= norm; d /= norm;

  bathy_plane.resize(4);
  bathy_plane[0] = a;
  bathy_plane[1] = b;
  bathy_plane[2] = c;
  bathy_plane[3] = d;
}

// Sample 3 points on the projected plane and fit a local ECEF plane
// approximation. Accurate within ~20 m (flat-Earth limit); used by the
// Newton-Raphson refraction solver to iterate entirely in ECEF and avoid
// expensive per-iteration proj/unproj round-trips. This assumes 
// the projection is in meters, not degrees.
std::vector<double>
fitLocalEcefPlaneToProjPlane(std::vector<double> const& plane,
                             vw::cartography::GeoReference const& plane_proj,
                             vw::Vector3 const& proj_pt,
                             double offset_meters) {

  // offset_meters is treated as a Euclidean length in plane_proj's axes, so
  // those axes must be meter-scale. Callers should pass
  // BathyPlane::stereographic_proj; passing a geographic georef would
  // interpret 1.0 as 1 degree (~111 km) and produce a bogus tangent.
  if (!plane_proj.is_projected())
    vw_throw(vw::LogicErr() << "fitLocalEcefPlaneToProjPlane requires a "
              << "projected (meter-scale) georef; got geographic instead.\n");

  // Project proj_pt onto the plane to get a point on the plane.
  double dist = signed_dist_to_plane(plane, proj_pt);
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
  vw::Vector3 ecef0 = vw::unproj_point(plane_proj, pt0);
  vw::Vector3 ecef1 = vw::unproj_point(plane_proj, pt1);
  vw::Vector3 ecef2 = vw::unproj_point(plane_proj, pt2);

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

// Bilinear-sample three raster heights one pixel apart, convert to ECEF, fit
// a plane. The neighbors are pixel-aligned (not proj-meter offsets) so the
// three samples always span a known scale regardless of the raster's CRS:
// in a 1-degree-per-unit geographic proj, a 1-pixel offset is still one
// pixel, not 1 degree. offset_meters is only forwarded to the plane-based
// fallback. Falls back to the plane-based fit if the center pixel is out of
// bounds, if both +1 and -1 neighbors are out of bounds in either axis, or
// if any sample value is invalid.
//
// Input proj_pt is in bp.stereographic_proj's frame (meter-scale). The
// raster is sampled through bp.plane_proj (its native georef, possibly
// geographic), so we convert stereographic -> ECEF -> raster coords
// before looking up the pixel.
std::vector<double> fitLocalEcefPlaneToProjSurface(BathyPlane const& bp,
                                                   vw::Vector3 const& proj_pt,
                                                   double offset_meters) {

  vw::ImageView<vw::PixelMask<float>> const& raster = bp.water_surface;
  int cols = raster.cols(), rows = raster.rows();

  // Convert stereographic proj_pt to raster-native coords for pixel lookup.
  vw::Vector3 ecef_at_query = vw::unproj_point(bp.stereographic_proj, proj_pt);
  vw::Vector3 raster_proj_pt = vw::proj_point(bp.plane_proj, ecef_at_query);

  // Center sample (pixel coords, possibly fractional).
  vw::Vector2 pix0 = bp.plane_proj.point_to_pixel(
                       vw::Vector2(raster_proj_pt[0], raster_proj_pt[1]));
  if (!isInBounds(pix0, cols, rows))
    return fitLocalEcefPlaneToProjPlane(bp.bathy_plane, bp.stereographic_proj,
                                        proj_pt, offset_meters);

  // x-neighbor: +1 pixel, falling back to -1 pixel at the high edge.
  vw::Vector2 pix1 = pix0 + vw::Vector2(1.0, 0.0);
  if (!isInBounds(pix1, cols, rows))
    pix1 = pix0 + vw::Vector2(-1.0, 0.0);
  if (!isInBounds(pix1, cols, rows))
    return fitLocalEcefPlaneToProjPlane(bp.bathy_plane, bp.stereographic_proj,
                                        proj_pt, offset_meters);

  // y-neighbor: +1 pixel, falling back to -1 pixel.
  vw::Vector2 pix2 = pix0 + vw::Vector2(0.0, 1.0);
  if (!isInBounds(pix2, cols, rows))
    pix2 = pix0 + vw::Vector2(0.0, -1.0);
  if (!isInBounds(pix2, cols, rows))
    return fitLocalEcefPlaneToProjPlane(bp.bathy_plane, bp.stereographic_proj,
                                        proj_pt, offset_meters);

  auto interp = vw::BilinearInterpolation::interpolator(raster);
  vw::PixelMask<float> v0 = interp(raster, pix0[0], pix0[1], 0);
  vw::PixelMask<float> v1 = interp(raster, pix1[0], pix1[1], 0);
  vw::PixelMask<float> v2 = interp(raster, pix2[0], pix2[1], 0);

  if (!vw::is_valid(v0) || !vw::is_valid(v1) || !vw::is_valid(v2))
    return fitLocalEcefPlaneToProjPlane(bp.bathy_plane, bp.stereographic_proj,
                                        proj_pt, offset_meters);

  // Convert (pixel -> proj_xy, sampled height) -> ECEF for each sample.
  vw::Vector2 p0 = bp.plane_proj.pixel_to_point(pix0);
  vw::Vector2 p1 = bp.plane_proj.pixel_to_point(pix1);
  vw::Vector2 p2 = bp.plane_proj.pixel_to_point(pix2);
  vw::Vector3 ecef0
    = vw::unproj_point(bp.plane_proj, vw::Vector3(p0[0], p0[1], v0.child()));
  vw::Vector3 ecef1
    = vw::unproj_point(bp.plane_proj, vw::Vector3(p1[0], p1[1], v1.child()));
  vw::Vector3 ecef2
    = vw::unproj_point(bp.plane_proj, vw::Vector3(p2[0], p2[1], v2.child()));

  // Exact plane through three ECEF points.
  vw::Vector3 v = ecef1 - ecef0;
  vw::Vector3 w = ecef2 - ecef0;
  vw::Vector3 normal_ecef = normalize(cross_prod(v, w));

  // Orient normal away from Earth center.
  if (dot_prod(normal_ecef, ecef0) < 0)
    normal_ecef = -normal_ecef;

  double D_ecef = dot_prod(normal_ecef, ecef0);
  std::vector<double> ecef_plane(4);
  ecef_plane[0] = normal_ecef[0];
  ecef_plane[1] = normal_ecef[1];
  ecef_plane[2] = normal_ecef[2];
  ecef_plane[3] = -D_ecef;
  return ecef_plane;
}

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

// Fit a local plane in stereographic_proj's meter-scale frame from 3
// raster neighbors near proj_pt. Mirrors fitLocalEcefPlaneToProjSurface's
// sampling logic but keeps the output plane in stereographic coords for
// use with rayBathyPlaneIntersect / Snell's law in proj coords. Returns
// an empty vector on any sampling failure.
std::vector<double> refineLocalPlaneFromRaster(
    BathyPlane const& bp,
    vw::Vector3 const& proj_pt) {

  vw::ImageView<vw::PixelMask<float>> const& raster = bp.water_surface;
  int cols = raster.cols(), rows = raster.rows();
  if (cols <= 0 || rows <= 0)
    return {};

  // Convert stereographic proj_pt to raster-native coords for pixel lookup.
  vw::Vector3 ecef_at_query = vw::unproj_point(bp.stereographic_proj, proj_pt);
  vw::Vector3 raster_proj_pt = vw::proj_point(bp.plane_proj, ecef_at_query);

  vw::Vector2 pix0 = bp.plane_proj.point_to_pixel(
                       vw::Vector2(raster_proj_pt[0], raster_proj_pt[1]));
  if (!isInBounds(pix0, cols, rows))
    return {};

  vw::Vector2 pix1 = pix0 + vw::Vector2(1.0, 0.0);
  if (!isInBounds(pix1, cols, rows))
    pix1 = pix0 + vw::Vector2(-1.0, 0.0);
  if (!isInBounds(pix1, cols, rows))
    return {};

  vw::Vector2 pix2 = pix0 + vw::Vector2(0.0, 1.0);
  if (!isInBounds(pix2, cols, rows))
    pix2 = pix0 + vw::Vector2(0.0, -1.0);
  if (!isInBounds(pix2, cols, rows))
    return {};

  auto interp = vw::BilinearInterpolation::interpolator(raster);
  vw::PixelMask<float> v0 = interp(raster, pix0[0], pix0[1], 0);
  vw::PixelMask<float> v1 = interp(raster, pix1[0], pix1[1], 0);
  vw::PixelMask<float> v2 = interp(raster, pix2[0], pix2[1], 0);
  if (!vw::is_valid(v0) || !vw::is_valid(v1) || !vw::is_valid(v2))
    return {};

  // Convert (pix_i, h_i) to stereographic coords via ECEF round-trip.
  auto toStereo = [&](vw::Vector2 const& pix, double h) -> vw::Vector3 {
    vw::Vector2 rp = bp.plane_proj.pixel_to_point(pix);
    vw::Vector3 ecef = vw::unproj_point(bp.plane_proj,
                                        vw::Vector3(rp[0], rp[1], h));
    return vw::proj_point(bp.stereographic_proj, ecef);
  };
  vw::Vector3 s0 = toStereo(pix0, v0.child());
  vw::Vector3 s1 = toStereo(pix1, v1.child());
  vw::Vector3 s2 = toStereo(pix2, v2.child());

  // Exact plane through the three points in stereographic coords.
  vw::Vector3 u = s1 - s0;
  vw::Vector3 w = s2 - s0;
  vw::Vector3 normal = cross_prod(u, w);
  double normal_norm = vw::math::norm_2(normal);
  if (normal_norm < 1e-12)
    return {};  // degenerate triangle, fall back to global plane
  normal /= normal_norm;
  if (normal[2] < 0) normal = -normal;   // orient +z

  std::vector<double> plane(4);
  plane[0] = normal[0];
  plane[1] = normal[1];
  plane[2] = normal[2];
  plane[3] = -dot_prod(normal, s0);
  return plane;
}

// Read a georeferenced raster of water-surface heights. The raster's own
// georef is stored as plane_proj (used only for pixel lookups). A local
// stereographic is derived at the valid-pixel centroid and stored as
// stereographic_proj. Valid pixels are reprojected to that stereographic
// frame via ECEF and passed to fitPlaneToPoints, so the four plane
// coefficients live in meter-scale coords regardless of the raster CRS.
// Downstream consumers pair bathy_plane with stereographic_proj.
static void readBathyPlaneFromRaster(std::string const& bathy_plane_file,
                                     BathyPlane & bp) {
  std::vector<double> & bathy_plane = bp.bathy_plane;
  vw::cartography::GeoReference & plane_proj = bp.plane_proj;

  // Read pixel data with nodata-aware mask. Use raw DiskImageView rather than
  // read_bathy_mask(), which also invalidates non-positive pixels - wrong for
  // water-surface heights (typically negative relative to WGS84 ellipsoid).
  float nodata = -std::numeric_limits<float>::max();
  vw::read_nodata_val(bathy_plane_file, nodata);
  vw::ImageView<vw::PixelMask<float>> raster
    = vw::create_mask(vw::DiskImageView<float>(bathy_plane_file), nodata);

  // Walk valid pixels, collect points in the raster's projection coordinates.
  // Also accumulate the sum of valid heights to compute the physical mean,
  // and the mean proj_xy to use as the stereographic origin.
  std::vector<vw::Vector3> raster_proj_pts;
  raster_proj_pts.reserve(size_t(raster.cols()) * size_t(raster.rows()));
  double height_sum = 0.0;
  vw::Vector2 proj_xy_sum(0, 0);
  for (int row = 0; row < raster.rows(); row++) {
    for (int col = 0; col < raster.cols(); col++) {
      vw::PixelMask<float> pix = raster(col, row);
      if (!vw::is_valid(pix))
        continue;
      vw::Vector2 proj_xy = plane_proj.pixel_to_point(vw::Vector2(col, row));
      raster_proj_pts.push_back(vw::Vector3(proj_xy.x(), proj_xy.y(), pix.child()));
      height_sum += pix.child();
      proj_xy_sum += proj_xy;
    }
  }
  if (raster_proj_pts.size() < 3)
    vw_throw(vw::IOErr() << "Bathy plane raster " << bathy_plane_file
              << " has fewer than 3 valid pixels; cannot fit a plane.\n");
  bp.mean_height = height_sum / double(raster_proj_pts.size());

  // Downstream plane math needs meter-scale axes. If the raster's native
  // georef is already stereographic, reuse it verbatim so the fitted
  // coefficients stay in that frame. Otherwise derive a local stereographic
  // centered at the valid-pixel centroid (covers geographic rasters like
  // Monica's EPSG:4326 wl.tifs).
  if (plane_proj.overall_proj4_str().find("+proj=stere") != std::string::npos) {
    bp.stereographic_proj = plane_proj;
  } else {
    vw::Vector2 center_proj_xy = proj_xy_sum / double(raster_proj_pts.size());
    vw::Vector2 center_lonlat = plane_proj.point_to_lonlat(center_proj_xy);
    bp.stereographic_proj.set_datum(plane_proj.datum());
    bp.stereographic_proj.set_stereographic(center_lonlat[1], center_lonlat[0], 1.0);
  }

  // Guard the invariant: stereographic_proj must be meter-scale regardless
  // of input path.
  if (!bp.stereographic_proj.is_projected())
    vw_throw(vw::LogicErr() << "BathyPlane::stereographic_proj must be a "
              << "projected (meter-scale) CRS; derivation from "
              << bathy_plane_file << " did not produce one.\n");

  // Re-express each sample in the stereographic frame via ECEF, so the
  // fitted bathy_plane coefficients live in meter-scale coords regardless
  // of the raster's native CRS. signed_dist_to_plane, rayPlaneIntersect,
  // and Snell's law in proj coords all then work unchanged downstream.
  std::vector<vw::Vector3> stereo_pts;
  stereo_pts.reserve(raster_proj_pts.size());
  for (auto const& p : raster_proj_pts) {
    vw::Vector3 ecef = vw::unproj_point(plane_proj, p);
    stereo_pts.push_back(vw::proj_point(bp.stereographic_proj, ecef));
  }
  fitPlaneToPoints(stereo_pts, bathy_plane);

  // Report fit residuals so the user can see whether a plane is a good model.
  double max_abs_residual = 0.0, sum_sq_residual = 0.0;
  for (auto const& p : stereo_pts) {
    double r = bathy_plane[0]*p[0] + bathy_plane[1]*p[1]
             + bathy_plane[2]*p[2] + bathy_plane[3];
    max_abs_residual = std::max(max_abs_residual, std::abs(r));
    sum_sq_residual += r * r;
  }
  double rms = std::sqrt(sum_sq_residual / double(stereo_pts.size()));

  vw_out() << "Fitted bathy plane from raster: "
           << bathy_plane[0] << " " << bathy_plane[1] << " "
           << bathy_plane[2] << " " << bathy_plane[3] << "\n";
  vw_out() << "Plane-fit residuals (in stereographic frame, meters): max "
           << max_abs_residual << ", RMS " << rms << "\n";
  vw_out() << "Projection: " << plane_proj.overall_proj4_str() << "\n";
  if (bp.stereographic_proj.overall_proj4_str() != plane_proj.overall_proj4_str())
    vw_out() << "Stereographic companion: "
             << bp.stereographic_proj.overall_proj4_str() << "\n";
}

// Dispatch: try to read as a georeferenced raster first; on failure, fall
// back to the legacy 3-line text format.
void readBathyPlane(std::string const& bathy_plane_file, BathyPlane & bp) {

  vw::vw_out() << "Reading bathy plane: " << bathy_plane_file << "\n";

  bool is_raster = false;
  try {
    is_raster = vw::cartography::read_georeference(bp.plane_proj, bathy_plane_file);
  } catch (...) {
    is_raster = false;
  }

  if (is_raster)
    readBathyPlaneFromRaster(bathy_plane_file, bp);
  else
    readBathyPlaneFromText(bathy_plane_file, bp);
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
    readBathyPlane(bathy_plane_file, bathy_plane_vec.back());
  }

  if (bathy_plane_vec.size() != 1 && bathy_plane_vec.size() != (size_t)num_images)
    vw_throw(vw::ArgumentErr() << "1 or " << num_images << " bathy planes expected.\n");

  // Clone the bathy plane if there's only one
  while (bathy_plane_vec.size() < (size_t)num_images)
    bathy_plane_vec.push_back(bathy_plane_vec[0]);
}

} // namespace vw
