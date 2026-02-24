// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
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

// Very small standalone utilities for working with georeferences and images

#include <vw/Cartography/GeoReferenceBaseUtils.h>
#include <vw/Image/MaskViews.h>
#include <vw/FileIO/DiskImageUtils.h>
#include <vw/FileIO/DiskImageResource.h>

namespace vw {
namespace cartography {

// Convert from projected coordinates to ECEF
vw::Vector3 projToEcef(GeoReference const& georef,
                       vw::Vector3  const& proj) {
  vw::Vector3 llh = georef.point_to_geodetic(proj);
  vw::Vector3 ecef = georef.datum().geodetic_to_cartesian(llh);
  return ecef;
}

// Convert from ECEF to projected coordinates
vw::Vector3 ecefToProj(GeoReference const& georef,
                       vw::Vector3  const& ecef) {
  vw::Vector3 llh = georef.datum().cartesian_to_geodetic(ecef);
  vw::Vector3 proj = georef.geodetic_to_point(llh);
  return proj;
}

// A function that will read a geo-referenced image, its nodata value,
// and the georeference, and will return a PixelMasked image, the nodata
// value, and the georeference.
void readGeorefImage(std::string const& image_file, 
  float & nodata_val, vw::cartography::GeoReference & georef,
  vw::ImageViewRef<vw::PixelMask<float>> & masked_image) {

  // Initial value, in case the image has no nodata field
  nodata_val = std::numeric_limits<float>::quiet_NaN();
  if (!vw::read_nodata_val(image_file, nodata_val))
        vw::vw_out() << "Warning: Could not read the nodata value for: "
                      << image_file << "\nUsing: " << nodata_val << ".\n";

  // Read the image
  vw::vw_out() << "Reading: " << image_file << std::endl;
  vw::DiskImageView<float> image(image_file);
  // Create the masked image
  masked_image = vw::create_mask(image, nodata_val);

  // Read the georeference, and throw an exception if it is missing
  bool has_georef = vw::cartography::read_georeference(georef, image_file);
  if (!has_georef)
    vw::vw_throw(vw::ArgumentErr() << "Missing georeference in: "
                                    << image_file << ".\n");
}

// Given a georeferenced image and a point in ECEF, return the closest
// pixel value in the image. If the point is outside the image, return
// an invalid pixel.
vw::PixelMask<float> closestPixelVal(vw::ImageViewRef<vw::PixelMask<float>> const& image,
                                     vw::cartography::GeoReference const& georef,
                                     vw::Vector3 const& ecef) {

  // Convert the point to lon-lat-height
  vw::Vector3 llh = georef.datum().cartesian_to_geodetic(ecef);
  
  // Convert to pixel
  vw::Vector2 pix = georef.lonlat_to_pixel(subvector(llh, 0, 2));
  
  // Find col and row by rounding
  int col = round(pix[0]);
  int row = round(pix[1]);
  
  // If outside bounds return invalid pixel
  if (col < 0 || row < 0 || col >= image.cols() || row >= image.rows())
    return vw::PixelMask<float>();
    
  return image(col, row);  
}

}} // vw::cartography
