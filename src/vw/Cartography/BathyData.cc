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
// vw/Cartography/BathyData.h: areMasked, read_bathy_mask(s),
// readBathyPlanes, and the ECEF / projected-coords primitives
// bathyProjPoint and bathyUnprojPoint. The actual ray-bending logic
// (curvedSnellLaw, rayBathyPlaneIntersect, local tangent fits, raster
// neighbor sampling) lives in BathyRay.cc.

#include <vw/Cartography/BathyData.h>

#include <vw/Cartography/Datum.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/SnellLaw.h>
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
#include <vw/Core/Exception.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/FileIO/DiskImageUtils.h>

#include <cmath>
#include <fstream>
#include <sstream>
#include <limits>

namespace vw {

namespace {

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
                      std::vector<double>& bathy_plane) {

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

// Compute the projected coordinates of an ECEF point.
Vector3 bathyProjPoint(vw::cartography::GeoReference const& projection,
                       Vector3 const& xyz) {
  return projection.geodetic_to_point(projection.datum().cartesian_to_geodetic(xyz));
}

// Reverse this operation.
Vector3 bathyUnprojPoint(vw::cartography::GeoReference const& projection,
                         Vector3 const& proj_pt) {
  return projection.datum().geodetic_to_cartesian(projection.point_to_geodetic(proj_pt));
}

// Given an ECEF point xyz and two bathy planes, find if xyz is above or below
// each plane by signed distance in each plane's stereographic frame.
void signedDistToPlanes(std::vector<BathyPlane> const& bathy_plane_vec,
                        vw::Vector3 const& xyz,
                        std::vector<double>& distances) {

  if (bathy_plane_vec.size() != 2)
    vw_throw(vw::ArgumentErr() << "Two bathy planes expected.\n");

  distances.resize(2);
  for (size_t it = 0; it < 2; it++) {
    // bathy_plane coefs live in stereographic_proj's frame, so project xyz
    // into that frame (meter-scale) before evaluating the signed distance.
    distances[it]
      = signedDistToPlane(bathy_plane_vec[it].bathy_plane,
                          bathyProjPoint(bathy_plane_vec[it].stereographic_proj, xyz));
  }
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

  // Guard the invariant: downstream plane math (signedDistToPlane,
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
  // centered at the valid-pixel centroid (covers geographic rasters).
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
  // of the raster's native CRS. signedDistToPlane, rayPlaneIntersect,
  // and Snell's law in proj coords all then work unchanged downstream.
  std::vector<vw::Vector3> stereo_pts;
  stereo_pts.reserve(raster_proj_pts.size());
  for (auto const& p : raster_proj_pts) {
    vw::Vector3 ecef = bathyUnprojPoint(plane_proj, p);
    stereo_pts.push_back(bathyProjPoint(bp.stereographic_proj, ecef));
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
