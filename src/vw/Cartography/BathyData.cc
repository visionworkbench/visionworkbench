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
#include <vw/Core/Exception.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/FileIO/DiskImageUtils.h>

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

} // namespace vw
