// __BEGIN_LICENSE__
//  Copyright (c) 2006-2025, United States Government as represented by the
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


#include <vw/Math/Matrix.h>
#include <vw/Core/Exception.h>
#include <vw/Cartography/GeoReferenceResourceGDAL.h>
#include <vw/FileIO/DiskImageResourceGDAL.h>
#include <vw/Cartography/GeoReference.h>

// TODO(oalexan1): Can probably trim down these headers
#include "gdal.h"
#include "gdal_priv.h"
#include "cpl_string.h"
#include "ogr_spatialref.h"
#include "ogr_api.h"

// For boost::split
#include <boost/algorithm/string.hpp>

namespace vw {
namespace cartography {

  // This is to ensure a given warning for a given file is printed only once.
  static std::set<std::string> pix_sign_warning, non_normal_georef_warning;
  static Mutex warning_mutex; // protect as shared resource

  bool read_gdal_georeference(GeoReference& georef,
                              DiskImageResourceGDAL const& resource) {
    boost::shared_ptr<GDALDataset> dataset = resource.get_dataset_ptr();
    if (!dataset)
      vw_throw(LogicErr() << "read_gdal_georeference: Could not read georeference. "
                << "No file has been opened.");
    // Pull the projection and datum information out of the file if available
    if (dataset->GetProjectionRef() != NULL)
      georef.set_wkt(dataset->GetProjectionRef());

    double geo_transform[6];
    Matrix<double,3,3> transform;
    if (dataset->GetGeoTransform(geo_transform) == CE_None) {
      transform.set_identity();
      transform(0,0) = geo_transform[1];
      transform(0,1) = geo_transform[2];
      transform(0,2) = geo_transform[0];
      transform(1,0) = geo_transform[4];
      transform(1,1) = geo_transform[5];
      transform(1,2) = geo_transform[3];

      // It is highly unusual for a georeference to have the y axis go up.
      // This breaks some assumptions in the code. Do not throw a fatal
      // error as sometimes the georeference is not used.
      if (transform(1, 1) > 0) {
        // Ensure the warning is printed once for each file.
        Mutex::Lock lock(warning_mutex);
        std::string fileName = resource.filename();
        if (pix_sign_warning.find(fileName) == pix_sign_warning.end()) {
          pix_sign_warning.insert(fileName);
          vw_out(WarningMessage)
            << "Found a georeference with a positive value of the y pixel component "
            << "in file: " << fileName << ". This is not standard. Incorrect results "
            << "may be produced. Check the pixel size with gdalinfo.\n";
        }
      }

      georef.set_transform(transform);

      // Determine the pixel interpretation for the image.  See the
      // comments in GeoReference for more information.  If nothing is
      // encoded in the file, the default is to assume PixelAsArea.
      georef.set_pixel_interpretation(GeoReference::PixelAsArea);
      char **metadata = dataset->GetMetadata();
      if (CSLCount(metadata) > 0) {
        for (int i = 0; metadata[i] != NULL; i++) {
          std::vector<std::string> split_vec;
          boost::split(split_vec, metadata[i], boost::is_any_of("="));
          if (split_vec[0] == GDALMD_AREA_OR_POINT && split_vec.size() >= 2)
            if (boost::trim_copy(split_vec[1]) == GDALMD_AOP_POINT)
              georef.set_pixel_interpretation(GeoReference::PixelAsPoint);
        }
      }
    } else {
      return false;
    }

    // Update the lon-lat georef bbox using the image size. This greatly
    // helps disambiguate the longitude when it comes to a 360 degree offset.
    int cols = resource.format().cols;
    int rows = resource.format().rows;
    georef.ll_box_from_pix_box(vw::BBox2(0, 0, cols, rows));

    // Georeference functions need not be invertible. When we perform a reverse
    // lookup (e.g. during a geotransform) we rely on PROJ.4 to pick one
    // possible value.  However, the georeference might actually place the image
    // at another (equivalent) location in the projected space.  This is hard to
    // recover from at run-time, and we currently have no generic solution to
    // this problem.  In the mean time, we at least test whether the
    // georeference is likely to run into this problem (and warn the user) by
    // checking whether forward- and reverse-projecting the four corner pixels
    // lands us back at the same pixel.
    std::vector<Vector2> test_pixels(4);
    test_pixels[0] = Vector2(0, 0);
    test_pixels[1] = Vector2(0, rows-1);
    test_pixels[2] = Vector2(cols-1, 0);
    test_pixels[3] = Vector2(cols-1, rows-1);
    for (int i = 0; i < 4; i++) {
      double error = 9999;
      bool have_error = false;
      try {
        error = georef.test_pixel_reprojection_error(test_pixels[i]);
      } catch (std::exception &e) {
        vw_out(WarningMessage) << e.what() << std::endl;
        have_error = true;
      }
      if (error > 0.1 || have_error) {
        // Print the warning just once per file
        Mutex::Lock lock(warning_mutex);
        std::string fileName = resource.filename();
        if (non_normal_georef_warning.find(fileName) == non_normal_georef_warning.end()) {
          non_normal_georef_warning.insert(fileName);
          vw_out(WarningMessage) << "read_gdal_georeference(): WARNING! Resource file " 
            << fileName << " contains a non-normal georeference.\n";
        }
      }
    }

    return true;
  }

  void write_gdal_georeference(DiskImageResourceGDAL& resource,
                               GeoReference const& georef) {

    boost::shared_ptr<GDALDataset> dataset = resource.get_dataset_ptr();
    if (!dataset)
      vw_throw(LogicErr() << "GeoReferenceHelperGDAL: Could not write georeference. "
               << "No file has been opened.");

    // Store the transform matrix
    double geo_transform[6] = { georef.transform()(0,2), georef.transform()(0,0),
                                georef.transform()(0,1), georef.transform()(1,2),
                                georef.transform()(1,0), georef.transform()(1,1) };
    dataset->SetGeoTransform(geo_transform);

    // This is a little ridiculous, but GDAL can't write geotiffs
    // without a string of Well Known Text (WKT), so we must conjure
    // up an OGRSpatialReference here and use it to convert from proj.4 to WKT.
    std::string wkt_str = georef.get_wkt();
    dataset->SetProjection(wkt_str.c_str());

    // Set the pixel interpretation for the image.  See the comments
    // in GeoReference for more information.
    if (georef.pixel_interpretation() == GeoReference::PixelAsArea)
      dataset->SetMetadataItem(GDALMD_AREA_OR_POINT, GDALMD_AOP_AREA);
    else // PixelAsPoint
      dataset->SetMetadataItem(GDALMD_AREA_OR_POINT, GDALMD_AOP_POINT);
  }

  // Read an arbitrary name = value pair from the geoheader.
  bool read_gdal_string(DiskImageResourceGDAL const& resource,
                        std::string const& str_name,
                        std::string & str_val) {

    // Initialize the output string
    str_val = "";

    // Call read_gdal_strings and then extract the value of the desired key.
    std::map<std::string, std::string> value_pairs;
    read_gdal_strings(resource, value_pairs);
    auto it = value_pairs.find(str_name);
    if (it != value_pairs.end()) {
      str_val = it->second;
      return true;
    }

    return false;
  }

  bool read_gdal_strings(DiskImageResourceGDAL const& resource,
                         std::map<std::string, std::string>& value_pairs) {

    // Wipe the output
    value_pairs.clear();

    boost::shared_ptr<GDALDataset>dataset = resource.get_dataset_ptr();
    if (!dataset)
      vw_throw(LogicErr() << "read_gdal_string: Could not read string. "
                << "No file has been opened.");

    char **metadata = dataset->GetMetadata();
    if (CSLCount(metadata) == 0)
      return false;

    for (int i = 0; metadata[i] != NULL; i++) {

      // Find the location of the first equal sign
      std::string line = metadata[i];
      auto loc = line.find("=");
      if (loc == std::string::npos)
        continue;
      // The part before the equal sign is the key
      std::string local_key = line.substr(0, loc);
      // The part after the equal sign is the value
      std::string local_val = line.substr(loc+1, line.size()-loc-1);
      value_pairs[local_key] = local_val;
    }
    return true;
  }

  // Write an arbitrary name=value pair in the geoheader.
  void write_gdal_string(DiskImageResourceGDAL& resource,
                          std::string const& str_name,
                          std::string const& str_val) {

    boost::shared_ptr<GDALDataset> dataset = resource.get_dataset_ptr();
    if (!dataset)
      vw_throw(LogicErr() << "write_gdal_string: Could not write string. "
                           << "No file has been opened.");

    dataset->SetMetadataItem(str_name.c_str(), str_val.c_str());
  }

}} // namespace vw::cartography
