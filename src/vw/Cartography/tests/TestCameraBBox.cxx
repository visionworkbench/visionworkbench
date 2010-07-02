// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestCameraBBox.h
#include <gtest/gtest.h>
#include <test/Helpers.h>

#include <vw/Cartography/CameraBBox.h>
#include <vw/Camera/PinholeModel.h>

using namespace vw;
using namespace vw::cartography;
using namespace vw::test;
using namespace vw::camera;

#if defined(VW_HAVE_PKG_PROTOBUF) && VW_HAVE_PKG_PROTOBUF==1
TEST( CameraBBox, GeospatialIntersectDatum ) {
  boost::shared_ptr<CameraModel> apollo( new PinholeModel("apollo.pinhole") );
  GeoReference moon;
  moon.set_well_known_geogcs("D_MOON");

  for ( unsigned i = 0; i < 10; i++ ) {
    Vector2 input_image( rand()%4096, rand()%4096 );
    bool did_intersect;
    Vector2 lonlat =
      geospatial_intersect( input_image, moon,
                            apollo, 1, did_intersect );
    ASSERT_TRUE( did_intersect );
    double radius = moon.datum().radius( lonlat[0], lonlat[1] );
    EXPECT_NEAR( radius, 1737400, 1e-3 );
    Vector3 llr( lonlat[0], lonlat[1], 0 );
    Vector3 ecef = moon.datum().geodetic_to_cartesian(llr);
    Vector3 llr2 = moon.datum().cartesian_to_geodetic(ecef);
    EXPECT_VECTOR_NEAR( llr2, llr, 1e-4 );

    Vector2 retrn_image = apollo->point_to_pixel( ecef );
    EXPECT_VECTOR_NEAR( retrn_image, input_image, 1e-3 );
  }
}

TEST( CameraBBox, CameraBBoxDatum ) {
  boost::shared_ptr<CameraModel> apollo( new PinholeModel("apollo.pinhole") );
  GeoReference moon;
  moon.set_well_known_geogcs("D_MOON");

  float scale;
  BBox2 image_bbox = camera_bbox( moon, apollo, 4096, 4096, scale );
  EXPECT_VECTOR_NEAR( image_bbox.min(), Vector2(86,-1), 2 );
  EXPECT_VECTOR_NEAR( image_bbox.max(), Vector2(95,7), 2 );
}
#endif
