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
#include <vw/Image/PixelMask.h>
#include <vw/Image/PixelMath.h>
#include <vw/Image/PixelTypeInfo.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/MaskViews.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/LinearAlgebra.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Core/Exception.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/FileIO/DiskImageUtils.h>

#include <cmath>
#include <fstream>
#include <sstream>
#include <limits>

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

// Read a plane from the legacy 3-line text format: four plane coefficients
// on line 1, comment on line 2, projection center lat/lon on line 3.
static void readBathyPlaneFromText(std::string const& bathy_plane_file,
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
  vw_out() << "Read projection: " <<  plane_proj.overall_proj4_str() << "\n";
}

// Best-fit plane through a set of 3D points via SVD. Returns the four
// coefficients (a, b, c, d) such that a*x + b*y + c*z + d = 0, with the
// normal (a, b, c) of unit length and oriented so c > 0. No outlier
// rejection - use on clean samples.
void fitPlaneToPoints(std::vector<vw::Vector3> const& points,
                      std::vector<double> & bathy_plane) {

  if (points.size() < 3)
    vw_throw(vw::ArgumentErr() << "fitPlaneToPoints: need at least 3 points.\n");

  // Centroid.
  vw::Vector3 centroid(0, 0, 0);
  for (auto const& p : points) centroid += p;
  centroid /= double(points.size());

  // SVD on the centered point matrix. The plane normal is the right-singular
  // vector for the smallest singular value, i.e. the last row of V^T.
  vw::Matrix<double> A(points.size(), 3);
  for (size_t i = 0; i < points.size(); i++) {
    A(i, 0) = points[i][0] - centroid[0];
    A(i, 1) = points[i][1] - centroid[1];
    A(i, 2) = points[i][2] - centroid[2];
  }
  vw::Matrix<double> U, VT;
  vw::Vector<double> S;
  vw::math::svd(A, U, S, VT);

  vw::Vector3 normal(VT(2, 0), VT(2, 1), VT(2, 2));
  if (normal[2] < 0) normal = -normal;
  double d = -dot_prod(normal, centroid);

  bathy_plane.resize(4);
  bathy_plane[0] = normal[0];
  bathy_plane[1] = normal[1];
  bathy_plane[2] = normal[2];
  bathy_plane[3] = d;
}

// Read a georeferenced raster of water-surface heights. The raster's own
// georef becomes plane_proj - no reprojection. Valid pixels are converted
// to (proj_x, proj_y, height) in that georef's coordinates and passed to
// fitPlaneToPoints to derive the four plane coefficients.
static void readBathyPlaneFromRaster(std::string const& bathy_plane_file,
                                     std::vector<double> & bathy_plane,
                                     vw::cartography::GeoReference & plane_proj) {

  // Read pixel data with nodata-aware mask. Use raw DiskImageView rather than
  // read_bathy_mask(), which also invalidates non-positive pixels - wrong for
  // water-surface heights (typically negative relative to WGS84 ellipsoid).
  float nodata = -std::numeric_limits<float>::max();
  vw::read_nodata_val(bathy_plane_file, nodata);
  vw::ImageView<vw::PixelMask<float>> raster
    = vw::create_mask(vw::DiskImageView<float>(bathy_plane_file), nodata);

  // Walk valid pixels, collect points in the raster's projection coordinates,
  // then find the best-fit plane we will use as an inital guess in solvers.
  std::vector<vw::Vector3> proj_pts;
  proj_pts.reserve(size_t(raster.cols()) * size_t(raster.rows()));
  for (int row = 0; row < raster.rows(); row++) {
    for (int col = 0; col < raster.cols(); col++) {
      vw::PixelMask<float> pix = raster(col, row);
      if (!vw::is_valid(pix))
        continue;
      vw::Vector2 proj_xy = plane_proj.pixel_to_point(vw::Vector2(col, row));
      proj_pts.push_back(vw::Vector3(proj_xy.x(), proj_xy.y(), pix.child()));
    }
  }
  if (proj_pts.size() < 3)
    vw_throw(vw::IOErr() << "Bathy plane raster " << bathy_plane_file
              << " has fewer than 3 valid pixels; cannot fit a plane.\n");
  fitPlaneToPoints(proj_pts, bathy_plane);

  // Report fit residuals so the user can see whether a plane is a good model.
  double max_abs_residual = 0.0, sum_sq_residual = 0.0;
  for (auto const& p : proj_pts) {
    double r = bathy_plane[0]*p[0] + bathy_plane[1]*p[1]
             + bathy_plane[2]*p[2] + bathy_plane[3];
    max_abs_residual = std::max(max_abs_residual, std::abs(r));
    sum_sq_residual += r * r;
  }
  double rms = std::sqrt(sum_sq_residual / double(proj_pts.size()));

  vw_out() << "Fitted bathy plane from raster: "
           << bathy_plane[0] << " " << bathy_plane[1] << " "
           << bathy_plane[2] << " " << bathy_plane[3] << "\n";
  vw_out() << "Plane-fit residuals: max " << max_abs_residual
           << ", RMS " << rms << "\n";
  vw_out() << "Projection: " << plane_proj.overall_proj4_str() << "\n";
}

// Dispatch: try to read as a georeferenced raster first; on failure, fall
// back to the legacy 3-line text format.
void readBathyPlane(std::string const& bathy_plane_file,
                    std::vector<double> & bathy_plane,
                    vw::cartography::GeoReference & plane_proj) {

  vw::vw_out() << "Reading bathy plane: " << bathy_plane_file << "\n";

  bool is_raster = false;
  try {
    is_raster = vw::cartography::read_georeference(plane_proj, bathy_plane_file);
  } catch (...) {
    is_raster = false;
  }

  if (is_raster)
    readBathyPlaneFromRaster(bathy_plane_file, bathy_plane, plane_proj);
  else
    readBathyPlaneFromText(bathy_plane_file, bathy_plane, plane_proj);
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

} // namespace vw
