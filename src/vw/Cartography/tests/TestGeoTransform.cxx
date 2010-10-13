// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestGeoTransform.h
#include <gtest/gtest.h>
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
