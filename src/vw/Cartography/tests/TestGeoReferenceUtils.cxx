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


// TestGeoReference.h
#include <gtest/gtest_VW.h>
#include <test/Helpers.h>

#include <vw/Cartography/GeoReferenceUtils.h>

#include <iomanip>

using namespace vw;
using namespace vw::cartography;
using namespace vw::test;

TEST(GeoReferenceUtils, CropAndResample) {
  Matrix3x3 affine;
  affine(0,0) = 0.01; // 100 pix/degree
  affine(1,1) = -0.01; // 100 pix/degree
  affine(2,2) = 1;
  affine(0,2) = 30;   // 30 deg east
  affine(1,2) = -35;  // 35 deg south
  GeoReference input_pa( Datum("WGS84"), affine, GeoReference::PixelAsArea );
  GeoReference crop_pa = crop( input_pa, 200, 400 );
  GeoReference resample_pa = resample( input_pa, 0.5, 2 );
  EXPECT_VECTOR_NEAR( input_pa.pixel_to_lonlat(Vector2(200, 400)),
                      crop_pa.pixel_to_lonlat(Vector2()), 1e-7 );
  EXPECT_VECTOR_NEAR( input_pa.pixel_to_lonlat(Vector2()),
                      crop_pa.pixel_to_lonlat(Vector2(-200,-400)), 1e-7 );
  EXPECT_VECTOR_NEAR( input_pa.pixel_to_lonlat(Vector2()),
                      resample_pa.pixel_to_lonlat(Vector2()), 1e-7 );
  EXPECT_VECTOR_NEAR( input_pa.pixel_to_lonlat(Vector2(100,100)),
                      resample_pa.pixel_to_lonlat(Vector2(50,200)), 1e-7 );

  GeoReference input_pp( Datum("WGS84"), affine, GeoReference::PixelAsPoint );
  GeoReference crop_pp = crop( input_pp, 200, 400 );
  GeoReference resample_pp = resample( input_pp, 0.5, 2 );
  EXPECT_VECTOR_NEAR( input_pp.pixel_to_lonlat(Vector2(200, 400)),
                      crop_pp.pixel_to_lonlat(Vector2()), 1e-7 );
  EXPECT_VECTOR_NEAR( input_pp.pixel_to_lonlat(Vector2()),
                      crop_pp.pixel_to_lonlat(Vector2(-200,-400)), 1e-7 );
  EXPECT_VECTOR_NEAR( input_pp.pixel_to_lonlat(Vector2()),
                      resample_pp.pixel_to_lonlat(Vector2()), 1e-7 );
  EXPECT_VECTOR_NEAR( input_pp.pixel_to_lonlat(Vector2(100,100)),
                      resample_pp.pixel_to_lonlat(Vector2(50,200)), 1e-7 );

  // Test with a projection that doesn't use degrees.
  input_pp.set_equirectangular( 40, 40, 50, 0, 0 );
  crop_pp =  crop( input_pp, 200, 400 );
  resample_pp = resample( input_pp, 0.5, 2 );
  EXPECT_VECTOR_NEAR( input_pp.pixel_to_lonlat(Vector2(200, 400)),
                      crop_pp.pixel_to_lonlat(Vector2()), 1e-7 );
  EXPECT_VECTOR_NEAR( input_pp.pixel_to_lonlat(Vector2()),
                      crop_pp.pixel_to_lonlat(Vector2(-200,-400)), 1e-7 );
  EXPECT_VECTOR_NEAR( input_pp.pixel_to_lonlat(Vector2()),
                      resample_pp.pixel_to_lonlat(Vector2()), 1e-7 );
  EXPECT_VECTOR_NEAR( input_pp.pixel_to_lonlat(Vector2(100,100)),
                      resample_pp.pixel_to_lonlat(Vector2(50,200)), 1e-7 );
}


TEST( GeoReferenceUtils, gdal_write_checks ) {

  // Init a georef, the numbers are pretty arbitrary, it just must be valid
  cartography::GeoReference georef;
  georef.set_geographic();
  georef.set_proj4_projection_str("+proj=longlat +a=3396190 +b=3396190 +no_defs ");
  georef.set_well_known_geogcs("D_MARS");
  Matrix3x3 affine;
  affine(0,0) = 0.01; // 100 pix/degree
  affine(1,1) = -0.01; // 100 pix/degree
  affine(2,2) = 1;
  affine(0,2) = 30;   // 30 deg east
  affine(1,2) = -35;  // 35 deg south
  georef.set_transform(affine);

  // For the test below to pass, the files must be present.
  ImageView<float> dem(100, 100);
  double nodata = -1000;
  bool has_nodata = true, has_georef = true;
  TerminalProgressCallback tpc("vw", "");
  GdalWriteOptions opt;

  block_write_gdal_image("dem.tif", dem, has_georef, georef, has_nodata, nodata, opt, tpc);
} 

TEST( GeoReferenceUtils, gdal_read_checks) {
  // Verify that our GDAL read function assigns the -180 to 180
  // longitude range to this image that could go either way.
  GeoReference georef;
  EXPECT_TRUE(read_georeference(georef, "dem.tif"));
  EXPECT_TRUE(georef.is_lon_center_around_zero());
  EXPECT_TRUE(georef.proj4_str().find("+over") == std::string::npos);
  
  // Check handling of a 0-360 image
  //EXPECT_TRUE(read_georeference(georef, "dem360.tif"));
  //EXPECT_FALSE(georef.is_lon_center_around_zero());
  //EXPECT_TRUE(georef.proj4_str().find("+over") != std::string::npos);

  // Handle an image that is larger than 360 degrees!
  //EXPECT_TRUE(read_georeference(georef, "/home/smcmich1/repo/BinaryBuilder/build_root/install/share/geoids/egm96-5.jp2"));
  //EXPECT_FALSE(georef.is_lon_center_around_zero());
  //EXPECT_TRUE(georef.proj4_str().find("+over") != std::string::npos);
}


TEST( GeoReferenceUtils, haversine_distance) {
  // Simple check like what we use to compute image meters per pixel
  
  double d1    = haversine_circle_distance(Vector2(-83.1074910, 39.7140823), Vector2(-79.8817279, 37.8158211));
  double d2    = haversine_circle_distance(Vector2(-83.1074910, 37.8158211), Vector2(-79.8817279, 39.7140823));
  double denom = sqrt(30977*30977 + 18229*18229);

  EXPECT_NEAR(9.74744, d1/denom, 0.0001);
  EXPECT_NEAR(9.74744, d2/denom, 0.0001);
}

