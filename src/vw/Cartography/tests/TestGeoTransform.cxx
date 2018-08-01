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


// TestGeoTransform.h
#include <gtest/gtest_VW.h>
#include <test/Helpers.h>

#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoTransform.h>
#include <vw/Cartography/GeoReferenceUtils.h>

using namespace vw;
using namespace vw::cartography;
using namespace vw::test;


TEST(GeoTransform, safe_set_center) {
  // This test will fail if the correct grid files are not installed!
  Matrix3x3 a; // Set an affine matrix for centering in the 0-360 range.
  a(0,0) = 6.44839e-06;  a(0,1) =  0;            a(0,2) = 243.745;
  a(1,0) = 0;            a(1,1) = -6.44839e-06;  a(1,2) = 35.2553;
  a(2,0) = 0;            a(2,1) = 0;             a(2,2) = 1;

  GeoReference in_georef, out_georef;
  in_georef.set_well_known_geogcs("NAD27");
  out_georef.set_well_known_geogcs("WGS84");
  in_georef.set_transform(a);
  out_georef.set_transform(a);
  in_georef.safe_set_lon_center(true); // Move center to -180 to 180 range.
  out_georef.safe_set_lon_center(true);
  
  GeoTransform tf(in_georef, out_georef);
  
  Vector2 pixel = tf.reverse(Vector2(100,100));
  EXPECT_VECTOR_NEAR(Vector2(232.2344,96.8230), pixel, 1e-4);
}


/*
// --> This test fails on some machines where the absolute path to
//     the .ct2 file is too long!
/// Try running a pixel-to-pixel transform with a +nadgrids file.
TEST(GeoTransform, nadgrids) {

  const std::string proj_str_in  = "+proj=longlat +datum=WGS84";
  const std::string proj_str_out = "+proj=longlat +datum=NAD83 +no_defs +nadgrids="TEST_SRCDIR"/wgs84_to_nad83.ct2";

  Matrix3x3 af;
  Datum d_in, d_out;
  d_in.set_datum_from_proj_str(proj_str_in);
  d_out.set_datum_from_proj_str(proj_str_out);

  af(0,0) = 6.44839e-06; af(0,1) =  0;           af(0,2) = 243.745;
  af(1,0) = 0;           af(1,1) = -6.44839e-06; af(1,2) = 35.2553;
  af(2,0) = 0;           af(2,1) = 0;            af(2,2) = 1;

  GeoReference georef1(d_in,  af);
  GeoReference georef2(d_out, af);

  GeoTransform trans12(georef1, georef2);

  Vector2 pixel(87, 72);
  Vector2 p1 = trans12.pixel_to_pixel(pixel);
 
  EXPECT_VECTOR_NEAR(p1,  Vector2(84.3167187348,71.0064820955), 1e-4);
}
*/

TEST( GeoTransform, BasicTransform ) {
  GeoReference src_georef;
  src_georef.set_well_known_geogcs("WGS84");

  GeoReference dst_georef;
  dst_georef.set_well_known_geogcs("WGS84");

  GeoTransform geotx(src_georef,dst_georef);

  Vector2 fwd = geotx.forward(Vector2(25,25)),
    rev = geotx.reverse(Vector2(25,25));

  EXPECT_VECTOR_NEAR( fwd, Vector2(25,25), 1e-16 );
  EXPECT_VECTOR_NEAR( rev, Vector2(25,25), 1e-16 );
}

TEST( GeoTransform, UTMFarZone ) {
  // This tests for a bug where forward_bbox calls latlon_to_* for a latlon
  // that is invalid for a utm zone.

  std::vector<Vector2> utm(4);
  utm[0] = Vector2(419832.648, 5184829.285); // UL
  utm[1] = Vector2(419832.648, 5129329.285); // LL
  utm[2] = Vector2(469332.648, 5184829.285); // UR
  utm[3] = Vector2(469332.648, 5129329.285); // LR

  Vector2i size(3300,3700);

  GeoReference ll_georef, utm_georef;

  Matrix3x3 utm_map;
  utm_map(0,0) =  (utm[2][0] - utm[0][0]) / static_cast<double>(size(0));
  utm_map(1,1) = -(utm[0][1] - utm[1][1]) / static_cast<double>(size(1));
  utm_map(0,2) = utm[0][0];
  utm_map(1,2) = utm[0][1];
  utm_map(2,2) = 1;
  utm_georef.set_transform(utm_map);
  utm_georef.set_pixel_interpretation(GeoReference::PixelAsPoint);
  // pick one far away from the poles
  utm_georef.set_UTM(59, false); // 59S

  Matrix3x3 ll_map;
  ll_map(0,0) =  360.0;
  ll_map(1,1) = -180.0;
  ll_map(0,2) = -180;
  ll_map(1,2) = 90;
  ll_map(2,2) = 1;
  ll_georef.set_transform(ll_map);

  GeoTransform geotx(utm_georef, ll_georef);

  BBox2i bbox, bbox2;

  // This should work
  EXPECT_NO_THROW( bbox = geotx.forward_bbox(BBox2i(0,0,size[0],size[1])) );

  // This should not throw an exception, even though we encounter a
  // proj.4 error because it tries to look up a lonlat that is outside
  // the utm zone, because this routine catches all exceptions.
  EXPECT_NO_THROW( bbox2 = geotx.reverse_bbox(BBox2i(0,0,size[0],size[1])));
}

TEST( GeoTransform, StereographicSingularity ) {
  // Test forward bbox actually hits the limit of platecarre
  GeoReference ll_georef, stereo_georef;

  Matrix3x3 transform = math::identity_matrix<3>();
  transform(1,1) = -1; transform(1,2) = 90;
  ll_georef.set_transform(transform);

  transform = math::identity_matrix<3>();
  transform(0,0) = transform(1,1) = 1e4;
  transform(0,2) = transform(1,2) = -1e6;
  stereo_georef.set_transform(transform);
  stereo_georef.set_stereographic(90,0,1);

  GeoTransform geotx(stereo_georef, ll_georef);
  BBox2i output, input(50,50,100,100);
  EXPECT_NO_THROW( output = geotx.forward_bbox(input) );
  EXPECT_NEAR( 0, output.min()[1], 2 );
}

TEST(GeoTransform, skipProjectionTest) {

  // Test GeoTransform functions when the skip_map_projection flag is set.
  Matrix3x3 affine;
  Datum d("WGS72");

  affine(0,0) =  0.01; // degrees per pixel
  affine(1,1) = -0.01; // degrees per pixel
  affine(2,2) = 1;
  affine(0,2) = 100.0; // Degrees
  affine(1,2) = 30.0;

  GeoReference georef1(d, affine, GeoReference::PixelAsPoint);

  affine(0,0) =  0.02; // degrees per pixel
  affine(1,1) = -0.02; // degrees per pixel
  
  GeoReference georef2(d, affine, GeoReference::PixelAsPoint);

  georef1.set_proj4_projection_str("+proj=longlat +ellps=WGS72 +no_defs");
  georef2.set_proj4_projection_str("+proj=longlat +ellps=WGS72 +no_defs");

  GeoTransform trans(georef1, georef2);

  // georef2 is half the pixel resolution as georef1, so the pixel coordinates
  // should be halved but the projected coordinates should be the same.

  Vector2 pixel1 = Vector2(10.0,20.0);
  Vector2 point1 = georef1.pixel_to_point(pixel1);
  
  Vector2 pixel2  = trans.point_to_pixel(point1);
  Vector2 point2  = trans.pixel_to_point(pixel1);
  Vector2 point2b = trans.point_to_point(point1); // Should be same result

  const double EPS = 1e-4;
  EXPECT_VECTOR_NEAR(pixel2, Vector2(5.0, 10.0), EPS);
  EXPECT_VECTOR_NEAR(point1, point2,  EPS);
  EXPECT_VECTOR_NEAR(point2, point2b, EPS);
}

TEST(GeoTransform, RefToRef) {
  
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

  GeoTransform trans(georef1, georef2);

  // Test a few GeoReference utility functions and verify the results are ok.
  // - These numbers were not carefully checked but they look ok.

  Vector2 pixel1(130,96);
  BBox2   pixel_bbox_1(130, 96, 11, 9);
  BBox2   point_bbox_1(246.75772, -74.952049, 0.001, 0.001);
  
  Vector2 point1 = georef1.pixel_to_point(pixel1);

  //std::cout << "pixel1 = "       << pixel1 << std::endl;
  //std::cout << "point1 = "       << point1 << std::endl;
  
  //Vector2 lonlat = georef1.point_to_lonlat(point1);
  //Vector2 point2 = georef2.lonlat_to_point(lonlat);
  //Vector2 pixel2 = georef2.point_to_pixel(point2);
  //std::cout << "lonlat = " << lonlat << std::endl;
  //std::cout << "point2 = " << point2 << std::endl;
  //std::cout << "pixel2 = " << pixel2 << std::endl;
  
  Vector2 point2       = trans.pixel_to_point     (pixel1      );
  Vector2 pixel2       = trans.point_to_pixel     (point1      );
  BBox2   pixel_bbox_2 = trans.point_to_pixel_bbox(point_bbox_1);
  BBox2   point_bbox_2 = trans.pixel_to_point_bbox(pixel_bbox_1);
  
  //std::cout << "point2 = "       << point2 << std::endl;
  //std::cout << "pixel2 = "       << pixel2 << std::endl;
  //std::cout << "pixel_bbox_2 = " << pixel_bbox_2 << std::endl;
  //std::cout << "point_bbox_2 = " << std::setprecision(12) << point_bbox_2 << std::endl;
  //std::cout << "point2 = " << std::setprecision(12) << point2 << std::endl;
  
  const double EPS = 1e-4;

  EXPECT_VECTOR_NEAR(point2, Vector2(-113.242425,-74.950475), EPS);
  EXPECT_VECTOR_NEAR(pixel2, Vector2(81.3333,74), EPS);
  EXPECT_VECTOR_NEAR(pixel_bbox_2.min(),  Vector2(82, 77),             EPS);
  EXPECT_VECTOR_NEAR(pixel_bbox_2.size(), Vector2(7,8),           EPS);
  EXPECT_VECTOR_NEAR(point_bbox_2.min(),  Vector2(-113.242425, -74.9517),   EPS);
  EXPECT_VECTOR_NEAR(point_bbox_2.size(), Vector2(0.00149999999999, 0.0012), EPS);
  EXPECT_VECTOR_NEAR(point_bbox_2.min(),  Vector2(-113.242425,-74.951675),    EPS);
  EXPECT_VECTOR_NEAR(point_bbox_2.size(), Vector2(0.0015,0.0012),             EPS);
}

TEST(GeoTransform, bboxTest) {

  Matrix3x3 affine;
  Datum d("D_MOON");

  affine(0,0) =  0.02; // degrees per pixel
  affine(1,1) = -0.02; // degrees per pixel
  affine(2,2) = 1;
  affine(0,2) = -1.0; // Degrees
  affine(1,2) = 20.0;
  
  // This georef will cross the 0 degree line
  GeoReference georef1(d, affine);

  // This one will not.
  affine(0,2) = 50.0; // Degrees
  affine(1,2) = 20.0; // 
  GeoReference georef2(d, affine);

  // This one is in the 0-360 space.
  affine(0,2) = 300.0; // Degrees
  affine(1,2) =  20.0; // 
  GeoReference georef3(d, affine);
  
  georef1.set_proj4_projection_str("+proj=longlat");
  georef2.set_proj4_projection_str("+proj=longlat");
  georef3.set_proj4_projection_str("+proj=longlat");

  BBox2 image_bbox(0, 0, 1000, 1000);  
  GeoTransform trans12(georef1, georef2, image_bbox, image_bbox);
  GeoTransform trans13(georef1, georef3, image_bbox, image_bbox);
  GeoTransform trans21(georef2, georef1, image_bbox, image_bbox);
  GeoTransform trans31(georef3, georef1, image_bbox, image_bbox);

  // The check_bbox_wraparound function should only trigger when the 
  //  lon center changes in the source georef.
  EXPECT_FALSE(trans12.check_bbox_wraparound());
  EXPECT_TRUE( trans13.check_bbox_wraparound());
  EXPECT_FALSE(trans21.check_bbox_wraparound());
  EXPECT_FALSE(trans31.check_bbox_wraparound());
}

/// Set up two sinusoidal GeoRefs and see if we can go back-forth to the same bbox
TEST(GeoTransform, sinusoidalTest) {

  Matrix3x3 affine;
  Datum d("D_EUROPA", "EUROPA", "EUROPA", 1564130, 1564130, 0);

  affine(0,0) =  233.0515264;
  affine(1,1) = -233.0515264;
  affine(2,2) = 1;
  affine(0,2) = -120254.588;
  affine(1,2) = 1676805.732;  
  GeoReference georef1(d, affine);

  affine(0,0) =  227.994294037;
  affine(1,1) = -227.994294037;
  affine(0,2) = -102825.427;
  affine(1,2) =  725933.832;
  GeoReference georef2(d, affine);
  
  georef1.set_proj4_projection_str("+proj=sinu +lon_0=132.64203131646 +x_0=0 +y_0=0 +a=1564130 +b=1564130 +units=m +no_defs");
  georef2.set_proj4_projection_str("+proj=sinu +lon_0=137.63572555431 +x_0=0 +y_0=0 +a=1564130 +b=1564130 +units=m +no_defs");

  BBox2 image_bbox(0,0,924,914);
  GeoTransform trans12(georef1, georef2, image_bbox, image_bbox);
  GeoTransform trans21(georef2, georef1, image_bbox, image_bbox);

  BBox2 point1 = trans21.pixel_to_point_bbox(image_bbox);
  BBox2 pixel2 = trans12.point_to_pixel_bbox(point1);
  
  double eps = 35.0; // Yes, our back-forth accuracy is this bad!
  EXPECT_NEAR(pixel2.min().x(),    0.0, eps);
  EXPECT_NEAR(pixel2.min().y(),    0.0, eps);
  EXPECT_NEAR(pixel2.max().x(), 924.0, eps);
  EXPECT_NEAR(pixel2.max().y(), 914.0, eps);  
}
