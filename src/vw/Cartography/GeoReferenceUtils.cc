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

#include <vw/Cartography/GeoReferenceUtils.h>
#include <vw/Cartography/GeoTransform.h>
#include <gdal_priv.h>
#include <gdal_utils.h>
#include <cpl_string.h>
#include <boost/filesystem.hpp>

namespace vw {
namespace cartography {

GeoReference crop(GeoReference const& input,
                  double upper_left_x, double upper_left_y,
                  double /*width*/, double /*height*/) {
  Vector2 top_left_ll;
  if (input.pixel_interpretation() == GeoReference::PixelAsArea) {
    top_left_ll = input.pixel_to_point(Vector2(upper_left_x, upper_left_y) - Vector2(0.5,0.5));
  } else {
    top_left_ll = input.pixel_to_point(Vector2(upper_left_x, upper_left_y));
  }
  GeoReference output = input;      // Start with copy of current transform
  Matrix3x3 T = output.transform();
  T(0,2) = top_left_ll[0];          // Shift the translation to the crop region
  T(1,2) = top_left_ll[1];          //  (don't need to worry about width/height)
  output.set_transform(T);
  // TODO(oalexan1): Need to update m_image_ll_box here! Must also use width and
  // height.
  return output;
}

GeoReference crop(GeoReference const& input, BBox2 const& bbox) {
  // Redirect to the other georeference crop call
  // TODO(oalexan1): Need to update m_image_ll_box here!
  return crop(input, bbox.min().x(), bbox.min().y(),
              bbox.width(), bbox.height());
}

GeoReference resample(GeoReference const& input, double scale_x, double scale_y) {
  GeoReference output = input;
  Matrix3x3 T = output.transform();
  T(0,0) /= scale_x;
  T(1,1) /= scale_y;
  if (input.pixel_interpretation() == GeoReference::PixelAsArea) {
    Vector2 top_left_ll = input.pixel_to_point(-Vector2(0.5 / scale_x, 0.5 / scale_y));
    T(0,2) = top_left_ll[0];
    T(1,2) = top_left_ll[1];
  }
  output.set_transform(T);
  return output;
}

GeoReference resample(GeoReference const& input, double scale) {
  return resample(input, scale, scale);
}

double haversine_circle_distance(Vector2 a, Vector2 b, double radius) {
  const double DEG_TO_RAD = M_PI / 180.0;
  double dlon = fabs(a[0] - b[0]);
  double dlat = fabs(a[1] - b[1]);
  double term1 = pow(sin(DEG_TO_RAD* dlat/2), 2.0);
  double term2 = cos(DEG_TO_RAD* a[1])*cos(DEG_TO_RAD* b[1])*pow(sin(DEG_TO_RAD* dlon/2), 2.0);
  double central_angle = 2.0 * asin(sqrt(term1 + term2)); // In radians
  return central_angle * radius;
}

double get_image_meters_per_pixel(int width, int height, GeoReference const& georef) {
  // Get lonlat coordinates of the four corners
  Vector2 lonlat_tl(georef.pixel_to_lonlat(Vector2(0,       0)));
  Vector2 lonlat_tr(georef.pixel_to_lonlat(Vector2(width-1, 0)));
  Vector2 lonlat_br(georef.pixel_to_lonlat(Vector2(width-1, height-1)));
  Vector2 lonlat_bl(georef.pixel_to_lonlat(Vector2(0,       height-1)));
  // Take the mean of the two diagonal distances.  If this does not prove robust
  // enough (say at the poles) we can add some more work here.
  double diag_pixel_dist = sqrt(height*height + width*width);
  double radius = georef.datum().semi_major_axis();
  double d1 = haversine_circle_distance(lonlat_tl, lonlat_br, radius) / diag_pixel_dist;
  double d2 = haversine_circle_distance(lonlat_tr, lonlat_bl, radius) / diag_pixel_dist;
  return (d1 + d2) / 2.0;
}

void convert_to_cog(const std::string& filename, GdalWriteOptions const& opt) {
  
  vw_out() << "Converting to cloud-optimized GeoTIFF.\n";
  
  // Open file to check data type
  GDALDataset* check_dataset = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
  if (check_dataset == NULL) {
    vw_out(WarningMessage) << "Failed to open file to check data type for COG conversion.\n";
    return;
  }
  GDALDataType dtype = check_dataset->GetRasterBand(1)->GetRasterDataType();
  GDALClose(check_dataset);
  
  // Determine appropriate predictor
  std::string predictor_opt;
  if (dtype == GDT_Float32 || dtype == GDT_Float64)
    predictor_opt = "PREDICTOR=3";  // Float predictor
  else
    predictor_opt = "PREDICTOR=2";
  
  // Create temp filename
  std::string temp_file = filename + ".tmp_cog.tif";
  
  // Move original to temp
  try {
    boost::filesystem::rename(filename, temp_file);
  } catch (const boost::filesystem::filesystem_error& e) {
    vw_out(WarningMessage) << "Failed to rename file for COG conversion: " << e.what() << "\n";
    return;
  }
  
  // Open source file
  GDALDataset* src_dataset = (GDALDataset*)GDALOpen(temp_file.c_str(), GA_ReadOnly);
  if (src_dataset == NULL) {
    vw_out(WarningMessage) << "Failed to open source file for COG conversion.\n";
    boost::filesystem::rename(temp_file, filename);
    return;
  }
  
  // Get BIGTIFF setting from user options
  std::string bigtiff = "IF_SAFER";
  auto it = opt.gdal_options.find("BIGTIFF");
  if (it != opt.gdal_options.end())
    bigtiff = it->second;
  
  // Set up GDALTranslate options for COG
  std::string bigtiff_opt = "BIGTIFF=" + bigtiff;
  char* argv[] = {
    const_cast<char*>("-of"), const_cast<char*>("COG"),
    const_cast<char*>("-co"), const_cast<char*>("COMPRESS=DEFLATE"),
    const_cast<char*>("-co"), const_cast<char*>(predictor_opt.c_str()),
    const_cast<char*>("-co"), const_cast<char*>("BLOCKSIZE=512"),
    const_cast<char*>("-co"), const_cast<char*>(bigtiff_opt.c_str()),
    const_cast<char*>("-co"), const_cast<char*>("NUM_THREADS=ALL_CPUS"),
    nullptr
  };
  
  vw_out() << "COG options: COMPRESS=DEFLATE, " << predictor_opt 
           << ", BLOCKSIZE=512, " << bigtiff_opt << "\n";
  
  GDALTranslateOptions* opts = GDALTranslateOptionsNew(argv, nullptr);
  if (opts == NULL) {
    vw_out(WarningMessage) << "Failed to create GDALTranslate options.\n";
    GDALClose(src_dataset);
    boost::filesystem::rename(temp_file, filename);
    return;
  }
  
  // Perform translation to COG
  GDALDatasetH dst_dataset = GDALTranslate(filename.c_str(), 
                                           (GDALDatasetH)src_dataset, 
                                           opts, nullptr);
  
  GDALTranslateOptionsFree(opts);
  GDALClose(src_dataset);
  
  if (dst_dataset == NULL) {
    vw_out(WarningMessage) << "COG conversion failed. Restoring original file.\n";
    boost::filesystem::rename(temp_file, filename);
    return;
  }
  
  GDALClose(dst_dataset);
  
  // Remove temp file
  try {
    boost::filesystem::remove(temp_file);
  } catch (const boost::filesystem::filesystem_error& e) {
    vw_out(WarningMessage) << "Failed to remove temporary file: " << e.what() << "\n";
  }
}

}} // vw::cartography
