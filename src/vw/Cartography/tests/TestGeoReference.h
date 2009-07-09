// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
    georef.set_pixel_interpretation(GeoReferenceBase::PixelAsPoint);
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
    georef.set_pixel_interpretation(GeoReferenceBase::PixelAsPoint);

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

  void test_utm_lonlat() {
    vector<Vector2> ll(5), utm(5), px(5);

    ll[0] = Vector2(170.008619281089, -43.4851542659474); // UL
    ll[1] = Vector2(170.000341300606, -43.9847965111766); // LL
    ll[2] = Vector2(170.620731735563, -43.4888289623561); // UR
    ll[3] = Vector2(170.617564317696, -43.9885355997305); // LR
    ll[4] = Vector2(170.311817924037, -43.7372482005704); // Center

    utm[0] = Vector2(419832.648, 5184829.285); // UL
    utm[1] = Vector2(419832.648, 5129329.285); // LL
    utm[2] = Vector2(469332.648, 5184829.285); // UR
    utm[3] = Vector2(469332.648, 5129329.285); // LR
    utm[4] = Vector2(444582.648, 5157079.285); // Center

    // Fake image size
    Vector2 size(3300,3700);
    px[0] = Vector2(0,0);
    px[1] = Vector2(0,size(1));
    px[2] = Vector2(size(0),0);
    px[3] = size;
    px[4] = size/2;

    GeoReference ll_georef, utm_georef;

    Matrix3x3 ll_map;
    ll_map(0,0) =  360.0;
    ll_map(1,1) = -180.0;
    ll_map(0,2) = -180;
    ll_map(1,2) = 90;
    ll_map(2,2) = 1;
    ll_georef.set_transform(ll_map);

    Matrix3x3 utm_map;
    utm_map(0,0) =  (utm[2][0] - utm[0][0]) / size(0);
    utm_map(1,1) = -(utm[0][1] - utm[1][1]) / size(1);
    utm_map(0,2) = utm[0][0];
    utm_map(1,2) = utm[0][1];
    utm_map(2,2) = 1;
    utm_georef.set_transform(utm_map);
    utm_georef.set_pixel_interpretation(GeoReference::PixelAsPoint);
    utm_georef.set_UTM(59, false); // 59S

    for (int i = 0; i < ll.size(); ++i) {
      Vector2 check_point  = utm_georef.pixel_to_point(px[i]);
      Vector2 check_lonlat = utm_georef.pixel_to_lonlat(px[i]);

      // At this latitude, 1 degree is ~111km lat and ~81km lon
      TS_ASSERT_DELTA(check_point(0),  utm[i](0), 1);    // Accurate within 1m
      TS_ASSERT_DELTA(check_point(1),  utm[i](1), 1);    // Accurate within 1m
      TS_ASSERT_DELTA(check_lonlat(0), ll[i](0),  1e-5); // Accurate within 0.8m
      TS_ASSERT_DELTA(check_lonlat(1), ll[i](1),  1e-5); // Accurate within 1.1m
    }
  }

};
