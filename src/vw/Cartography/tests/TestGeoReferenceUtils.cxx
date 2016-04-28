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

TEST(GeoReferenceUtils, RefToRef) {
  
  // Set up a pair of GeoReference objects which overlap but
  //  do not use the same longitude convention.
  
  Matrix3x3 affine;
  Datum d("WGS72");

  affine(0,0) =  0.00015; // meters per pixel
  affine(1,1) = -0.00015; // meters per pixel
  affine(2,2) = 1;
  affine(0,2) = 246.738; // Degrees
  affine(1,2) = -74.936;
  
  GeoReference georef1(d, affine);

  affine(0,2) = -113.2547; // Degrees
  affine(1,2) = -74.9393;
  
  GeoReference georef2(d, affine);
  
  georef1.set_proj4_projection_str("+proj=longlat +ellps=WGS72 +no_defs");
  georef2.set_proj4_projection_str("+proj=longlat +ellps=WGS72 +no_defs");

  // Test a few GeoReference utility functions and verify the results are ok.
  // - These numbers were not carefully checked but they look ok.

  Vector2 pixel1(130,96);
  BBox2   pixel_bbox_1(130, 96, 11, 9);
  BBox2   point_bbox_1(246.75772, -74.952049, 0.001, 0.001);
  
  Vector2 point1 = georef1.pixel_to_point(pixel1);
  //std::cout << "pixel1 = "       << pixel1 << std::endl;
  //std::cout << "point1 = "       << point1 << std::endl;
  
  Vector2 lonlat = georef1.point_to_lonlat(point1);
  Vector2 point2 = georef2.lonlat_to_point(lonlat);
  //Vector2 pixel2 = georef2.point_to_pixel(point2);
  //std::cout << "lonlat = " << lonlat << std::endl;
  //std::cout << "point2 = " << point2 << std::endl;
  //std::cout << "pixel2 = " << pixel2 << std::endl;
  
  Vector2 pixel2       = georef_point_to_georef_pixel     (point1,       georef1, georef2);
  BBox2   pixel_bbox_2 = georef_point_to_georef_pixel_bbox(point_bbox_1, georef1, georef2);
  BBox2   point_bbox_2 = georef_pixel_to_georef_point_bbox(pixel_bbox_1, georef1, georef2);
  
  //std::cout << "pixel2 = "       << pixel2 << std::endl;
  //std::cout << "pixel_bbox_2 = " << pixel_bbox_2 << std::endl;
  //std::cout << "point_bbox_2 = " << std::setprecision(12) << point_bbox_2 << std::endl;
  
  const double EPS = 1e-4;

  EXPECT_VECTOR_NEAR(pixel2, Vector2(81.3333,74), EPS);
  EXPECT_VECTOR_NEAR(pixel_bbox_2.min(),  Vector2(82,77), EPS);
  EXPECT_VECTOR_NEAR(pixel_bbox_2.size(), Vector2(7,8), EPS);
  EXPECT_VECTOR_NEAR(point_bbox_2.min(),  Vector2(-113.242425, -74.951825), EPS);
  EXPECT_VECTOR_NEAR(point_bbox_2.size(), Vector2(0.0015, 0.00135), EPS);

}




