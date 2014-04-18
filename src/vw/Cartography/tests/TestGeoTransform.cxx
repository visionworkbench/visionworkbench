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

using namespace vw;
using namespace vw::cartography;
using namespace vw::test;

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

  // This should throw a proj.4 error because it tries to look up a lonlat
  // that is outside the utm zone
  EXPECT_THROW( bbox2 = geotx.reverse_bbox(BBox2i(0,0,size[0],size[1])) ,
                vw::ArgumentErr );
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

/*
TEST( GeoTransform, LargeEqcTransform ) {

  // Construct two GeoReference objects
  // - They are identical except for their transform offsets.
  GeoReference georef1, georef2;
  Matrix3x3 affine1, affine2;
  georef1.set_equirectangular(0, 0, 24.1, 0, 0);
  georef1.set_well_known_geogcs("D_MOON");
  georef2.set_equirectangular(0, 0, 24.1, 0, 0);
  georef2.set_well_known_geogcs("D_MOON");

  // Set up a transform so the projected coordinates
  affine1(0,0) =  1.0;
  affine1(1,1) = -1.0;
  affine1(2,2) = 1;
  affine1(0,2) = -1.35197e+06;
  affine1(1,2) = 731470;
  georef1.set_transform(affine1);

  affine2(0,0) =  1.0;
  affine2(1,1) = -1.0;
  affine2(2,2) = 1;
  affine2(0,2) = -1.35176e+06;
  affine2(1,2) = 731352;
  georef2.set_transform(affine2);

  std::cout << georef1 << std::endl << std::endl;
  std::cout << georef2 << std::endl;

  GeoTransform geotx(georef1, georef2);
  std::cout << geotx.reverse(Vector2(0,    0)) << std::endl;
  std::cout << geotx.reverse(Vector2(0,  256)) << std::endl;
  std::cout << geotx.reverse(Vector2(256,512)) << std::endl;
  EXPECT_TRUE(false);


  left georef   = -- Proj.4 Geospatial Reference Object --
    Transform  : Matrix3x3((1,0,-1.35197e+06)(0,-1,731470)(0,0,1))
    Geodeditic Datum --> Name: unknown  Spheroid: unnamed  Semi-major: 1.7374e+06  Semi-minor: 1.7374e+06  Meridian: Greenwich  at 0
    Proj.4 String: +proj=eqc +lat_ts=24.1 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +units=m
    Pixel Interpretation: pixel as area

  right georef  = -- Proj.4 Geospatial Reference Object --
    Transform  : Matrix3x3((1,0,-1.35176e+06)(0,-1,731352)(0,0,1))
    Geodeditic Datum --> Name: unknown  Spheroid: unnamed  Semi-major: 1.7374e+06  Semi-minor: 1.7374e+06  Meridian: Greenwich  at 0
    Proj.4 String: +proj=eqc +lat_ts=24.1 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +units=m
    Pixel Interpretation: pixel as area


  $1 = (const vw::BBox2i &) @0x7fffe41cc860: {<vw::math::BBoxBase<vw::math::BBox<int, 2ul>, int, 2ul>> = {m_min =
      {<vw::math::VectorBase<vw::math::Vector<int, 2ul> >> = {<No data fields>}, core_ = {elems = {0, 256}}}, m_max =
      {<vw::math::VectorBase<vw::math::Vector<int, 2ul> >> = {<No data fields>}, core_ = {elems = {256, 512}}}}, <No data fields>}

  (gdb) info locals
  transformed_bbox = {<vw::math::BBoxBase<vw::math::BBox<int, 2ul>, int, 2ul>> = {m_min =
      {<vw::math::VectorBase<vw::math::Vector<int, 2ul> >> = {<No data fields>}, core_ = {elems = {0, 138}}}, m_max =
      {<vw::math::VectorBase<vw::math::Vector<int, 2ul> >> = {<No data fields>}, core_ = {elems = {9964869, 395}}}}, <No data fields>}


}
*/




