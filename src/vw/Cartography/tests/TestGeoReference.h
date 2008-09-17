// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
//
// Copyright 2006 Carnegie Mellon University. All rights reserved.
//
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
//
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
//
// __END_LICENSE__

// TestGeoReference.h
#include <cxxtest/TestSuite.h>

#include <vw/Cartography/GeoReference.h>

using namespace std;
using namespace vw;
using namespace vw::cartography;

class TestGeoReference : public CxxTest::TestSuite
{
public:

  void test_basic_geographic()
  {
    GeoReference georef;
    georef.set_well_known_geogcs("WGS84");

    // Start with some basic test of a geographic (unprojected)
    // coordinate system.
    Vector2 loc = georef.pixel_to_point(Vector2(0,0));
    Vector2 pix = georef.point_to_pixel(Vector2(0,0));
    TS_ASSERT_EQUALS(loc[0], 0);
    TS_ASSERT_EQUALS(loc[1], 0);
    TS_ASSERT_EQUALS(pix[0], 0);
    TS_ASSERT_EQUALS(pix[1], 0);

    loc = georef.pixel_to_point(Vector2(90,90));
    pix = georef.point_to_pixel(Vector2(90,90));
    TS_ASSERT_EQUALS(loc[0], 90);
    TS_ASSERT_EQUALS(loc[1], 90);
    TS_ASSERT_EQUALS(pix[0], 90);
    TS_ASSERT_EQUALS(pix[1], 90);

    loc = georef.pixel_to_point(Vector2(-90,-90));
    pix = georef.point_to_pixel(Vector2(-90,-90));
    TS_ASSERT_EQUALS(loc[0], -90);
    TS_ASSERT_EQUALS(loc[1], -90);
    TS_ASSERT_EQUALS(pix[0], -90);
    TS_ASSERT_EQUALS(pix[1], -90);
  }

  void test_basic_affine_transform()
  {
    GeoReference georef;
    georef.set_well_known_geogcs("WGS84");
    Matrix3x3 affine;
    affine(0,0) = 0.01; // 100 pix/degree
    affine(1,1) = -0.01; // 100 pix/degree
    affine(2,2) = 1;
    affine(0,2) = 30;   // 30 deg east
    affine(1,2) = -35;  // 35 deg south
    georef.set_transform(affine);

    // Test that, given a pixel location, you can recover the position
    // in the projected space. (remember, the VW assumes that the
    // position (0,0) appears at the center of the upper left hand
    // pixel in an image).
    georef.set_pixel_interpretation(GeoReference::PixelAsPoint);
    Vector2 loc = georef.pixel_to_point(Vector2(0,0));
    TS_ASSERT_EQUALS(loc[0], 30);
    TS_ASSERT_EQUALS(loc[1], -35);

    // Here, we are essentially querying for location (-0.5,-0.5)
    georef.set_pixel_interpretation(GeoReference::PixelAsArea);
    loc = georef.pixel_to_point(Vector2(-0.5,-0.5));
    TS_ASSERT_EQUALS(loc[0], 30 );
    TS_ASSERT_EQUALS(loc[1], -35 );

    // Test that, given a point in the projected space, you can
    // recover the proper pixel location.  (remember, the VW assumes
    // that the position (0,0) appears at the center of the upper left
    // hand pixel in an image).
    georef.set_pixel_interpretation(GeoReference::PixelAsPoint);
    Vector2 pix = georef.point_to_pixel(Vector2(30,-35));
    TS_ASSERT_EQUALS(pix[0], 0);
    TS_ASSERT_EQUALS(pix[1], 0);

    georef.set_pixel_interpretation(GeoReference::PixelAsArea);
    pix = georef.point_to_pixel(Vector2(30,-35));
    TS_ASSERT_EQUALS(pix[0], -0.5);
    TS_ASSERT_EQUALS(pix[1], -0.5);
  }

  void test_basic_projections()
  {
    // Set up a spherical datum for testing purposes.
    GeoReference georef;
    Datum d = georef.datum();
    d.set_semi_minor_axis(d.semi_major_axis());
    georef.set_datum(d);

    // -------------------------------------
    // Test an unprojected coordinate system
    // -------------------------------------
    Vector2 loc = georef.pixel_to_point(Vector2(500,500));
    TS_ASSERT_EQUALS(loc[0],500);
    TS_ASSERT_EQUALS(loc[1],500);

    Vector2 pix = georef.point_to_pixel(Vector2(500,500));
    TS_ASSERT_EQUALS(pix[0],500);
    TS_ASSERT_EQUALS(pix[1],500);

    Vector2 lon_lat = georef.point_to_lonlat(Vector2(140,89));
    TS_ASSERT_DELTA(lon_lat[0],140,1e-7);
    TS_ASSERT_DELTA(lon_lat[1],89,1e-7);

    loc = georef.lonlat_to_point(lon_lat);
    TS_ASSERT_DELTA(loc[0],140,1e-7);
    TS_ASSERT_DELTA(loc[1],89,1e-7);

    // ----------------------------------
    // Test a projected coordinate system
    // ----------------------------------
    georef.set_sinusoidal(0,0,0); // center longitude = 0

    loc = georef.pixel_to_point(Vector2(500,500));
    TS_ASSERT_EQUALS(loc[0],500);
    TS_ASSERT_EQUALS(loc[1],500);

    pix = georef.point_to_pixel(Vector2(500,500));
    TS_ASSERT_EQUALS(pix[0],500);
    TS_ASSERT_EQUALS(pix[1],500);

    // Whip up a home-brew sinusoidal projection for testing against
    // proj.4
    double lon = 20;
    double lat = 15;
    double meters_per_degree = d.semi_major_axis() * M_PI * 2.0/360.0;
    double x = lon*cos(lat*M_PI/180) * meters_per_degree; // sinusoidal projection
    double y = lat * meters_per_degree;

    lon_lat = Vector2(lon,lat);
    loc = georef.lonlat_to_point(lon_lat);
    TS_ASSERT_DELTA(loc[0],x,2);
    TS_ASSERT_DELTA(loc[1],y,2);

    lon_lat = georef.point_to_lonlat(loc);
    TS_ASSERT_DELTA(lon_lat[0],lon,1e-7);
    TS_ASSERT_DELTA(lon_lat[1],lat,1e-7);

    double x_1 = 7000/meters_per_degree;
    double y_1 = 3000/meters_per_degree;
    double lon_1 = x_1/cos(y_1*M_PI/180);
    double lat_1 = y_1;

    lon_lat = georef.point_to_lonlat(Vector2(7000,3000));
    TS_ASSERT_DELTA(lon_lat[0],lon_1,1e-7);
    TS_ASSERT_DELTA(lon_lat[1],lat_1,1e-7);

    loc = georef.lonlat_to_point(lon_lat);
    TS_ASSERT_DELTA(loc[0],7000,1e-7);
    TS_ASSERT_DELTA(loc[1],3000,1e-7);
  }


};
