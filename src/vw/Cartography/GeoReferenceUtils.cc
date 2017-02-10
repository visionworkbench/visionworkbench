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

namespace vw {
namespace cartography {

GdalWriteOptions::GdalWriteOptions() {
#if defined(VW_HAS_BIGTIFF) && VW_HAS_BIGTIFF == 1
  gdal_options["BIGTIFF"]  = "YES";
#else
  gdal_options["BIGTIFF"] = "NO";
#endif
  gdal_options["COMPRESS"] = "LZW";
  raster_tile_size = Vector2i(vw_settings().default_tile_size(),
                              vw_settings().default_tile_size());
  num_threads = vw_settings().default_num_threads();
}

GdalWriteOptionsDescription::GdalWriteOptionsDescription( GdalWriteOptions& opt ) {
  namespace po = boost::program_options;
  (*this).add_options()
    ("threads",      po::value(&opt.num_threads)->default_value(0),
        "Select the number of processors (threads) to use.")
    ("tile-size",  po::value(&opt.raster_tile_size)->default_value
     (Vector2i(vw_settings().default_tile_size(),
               vw_settings().default_tile_size()),"256, 256"),
        "Image tile size used for multi-threaded processing.")
    ("no-bigtiff",   "Tell GDAL to not create bigtiffs.")  // gets stored in vm.count("no-bigtiff")
    ("tif-compress", po::value(&opt.tif_compress)->default_value("LZW"),
        "TIFF Compression method. [None, LZW, Deflate, Packbits]")
    ("version,v",    "Display the version of software.")
    ("help,h",       "Display this help message.");
}

GeoReference crop( GeoReference const& input,
                   double upper_left_x, double upper_left_y,
                   double /*width*/, double /*height*/ ) {
  Vector2 top_left_ll;
  if ( input.pixel_interpretation() == GeoReference::PixelAsArea ) {
    top_left_ll = input.pixel_to_point( Vector2(upper_left_x, upper_left_y ) - Vector2(0.5,0.5) );
  } else {
    top_left_ll = input.pixel_to_point( Vector2(upper_left_x, upper_left_y ) );
  }
  GeoReference output = input;      // Start with copy of current transform
  Matrix3x3 T = output.transform();
  T(0,2) = top_left_ll[0];          // Shift the translation to the crop region
  T(1,2) = top_left_ll[1];          //  (don't need to worry about width/height)
  output.set_transform(T);
  return output;
}

GeoReference crop( GeoReference const& input,
                   BBox2 const& bbox ) {
  // Redirect to the other georeference crop call
  return crop(input, bbox.min().x(), bbox.min().y(),
              bbox.width(), bbox.height());
}

GeoReference resample( GeoReference const& input, double scale_x, double scale_y ) {
  GeoReference output = input;
  Matrix3x3 T = output.transform();
  T(0,0) /= scale_x;
  T(1,1) /= scale_y;
  if ( input.pixel_interpretation() == GeoReference::PixelAsArea ) {
    Vector2 top_left_ll = input.pixel_to_point( -Vector2(0.5 / scale_x, 0.5 / scale_y) );
    T(0,2) = top_left_ll[0];
    T(1,2) = top_left_ll[1];
  }
  output.set_transform(T);
  return output;
}

GeoReference resample( GeoReference const& input, double scale ) {
  return resample(input, scale, scale );
}

double haversine_circle_distance(Vector2 a, Vector2 b) {
  const double DEG_TO_RAD = M_PI / 180.0; // TODO: Where is this?
  double dlon = fabs(a[0] - b[0]);
  double dlat = fabs(a[1] - b[1]);
  double term1 = pow(sin(DEG_TO_RAD* dlat/2), 2.0);
  double term2 = cos(DEG_TO_RAD* a[1])*cos(DEG_TO_RAD* b[1])*pow(sin(DEG_TO_RAD* dlon/2), 2.0);
  double central_angle = 2.0 * asin(sqrt(term1 + term2)); // In radians
  double APPROX_RADIUS = 6371008.8; // In meters TODO: Collect constants somewhere!
  return central_angle * APPROX_RADIUS;
}


double get_image_meters_per_pixel(int width, int height, GeoReference const& georef) {
  // Get lonlat coordinates of the four corners
  Vector2 lonlat_tl(georef.pixel_to_lonlat(Vector2(0,       0       )));
  Vector2 lonlat_tr(georef.pixel_to_lonlat(Vector2(width-1, 0       )));
  Vector2 lonlat_br(georef.pixel_to_lonlat(Vector2(width-1, height-1)));
  Vector2 lonlat_bl(georef.pixel_to_lonlat(Vector2(0,       height-1)));
  // Take the mean of the two diagonal distances.  If this does not prove robust
  // enough (say at the poles) we can add some more work here.
  double diag_pixel_dist = sqrt(height*height + width*width);
  double d1 = haversine_circle_distance(lonlat_tl, lonlat_br) / diag_pixel_dist;
  double d2 = haversine_circle_distance(lonlat_tr, lonlat_bl) / diag_pixel_dist;
  return (d1 + d2) / 2.0;
}

}} // vw::cartography

#undef CHECK_PROJ_ERROR
