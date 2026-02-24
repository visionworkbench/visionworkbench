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

#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoReferenceUtils.h>
#include <vw/FileIO/GdalWriteOptions.h>

using namespace vw;
using namespace vw::cartography;
using namespace vw::test;

TEST( GeoReference, Core) {
  GeoReference georef;
  georef.set_pixel_interpretation(GeoReference::PixelAsPoint);
  georef.set_well_known_geogcs("WGS84");

  Matrix3x3 affine;
  affine(0,0) = 0.01; // 100 pix/degree
  affine(1,1) = -0.01; // 100 pix/degree
  affine(2,2) = 1;
  affine(0,2) = 30;   // 30 deg east
  affine(1,2) = -35;  // 35 deg south
  georef.set_transform(affine);

  Vector2 pix = georef.point_to_pixel(Vector2(30,-35));
  EXPECT_VECTOR_DOUBLE_EQ( pix, Vector2() );

  georef.set_pixel_interpretation(GeoReference::PixelAsArea);
  pix = georef.point_to_pixel(Vector2(30,-35));
  EXPECT_VECTOR_DOUBLE_EQ( pix, Vector2(-0.5,-0.5) );
}

TEST( GeoReference, BasicGeographic ) {
  GeoReference georef;
  georef.set_pixel_interpretation(GeoReference::PixelAsPoint);
  georef.set_well_known_geogcs("WGS84");

  // Start with some basic test of a geographic (unprojected)
  // coordinate system.
  Vector2 loc = georef.pixel_to_point(Vector2(0,0));
  Vector2 pix = georef.point_to_pixel(Vector2(0,0));
  EXPECT_VECTOR_DOUBLE_EQ( loc, Vector2() );
  EXPECT_VECTOR_DOUBLE_EQ( pix, Vector2() );

  loc = georef.pixel_to_point(Vector2(90,90));
  pix = georef.point_to_pixel(Vector2(90,90));
  EXPECT_VECTOR_DOUBLE_EQ( loc, Vector2(90,90) );
  EXPECT_VECTOR_DOUBLE_EQ( pix, Vector2(90,90) );

  loc = georef.pixel_to_point(Vector2(-90,-90));
  pix = georef.point_to_pixel(Vector2(-90,-90));
  EXPECT_VECTOR_DOUBLE_EQ( loc, Vector2(-90,-90) );
  EXPECT_VECTOR_DOUBLE_EQ( pix, Vector2(-90,-90) );
}

TEST( GeoReference, AffineTransform ) {
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
  EXPECT_VECTOR_DOUBLE_EQ( loc, Vector2(30,-35) );

  // Here, we are essentially querying for location (-0.5,-0.5)
  georef.set_pixel_interpretation(GeoReference::PixelAsArea);
  loc = georef.pixel_to_point(Vector2(-0.5,-0.5));
  EXPECT_VECTOR_DOUBLE_EQ( loc, Vector2(30,-35) );

  // Test that, given a point in the projected space, you can
  // recover the proper pixel location.  (remember, the VW assumes
  // that the position (0,0) appears at the center of the upper left
  // hand pixel in an image).
  georef.set_pixel_interpretation(GeoReference::PixelAsPoint);
  Vector2 pix = georef.point_to_pixel(Vector2(30,-35));
  EXPECT_VECTOR_DOUBLE_EQ( pix, Vector2() );

  georef.set_pixel_interpretation(GeoReference::PixelAsArea);
  pix = georef.point_to_pixel(Vector2(30,-35));
  EXPECT_VECTOR_DOUBLE_EQ( pix, Vector2(-0.5,-0.5) );
}

TEST( GeoReference, BasicProjections ) {
  // Set up a spherical datum for testing purposes.
  GeoReference georef;
  georef.set_pixel_interpretation(GeoReference::PixelAsPoint);

  Datum d = georef.datum();
  d.set_well_known_datum("D_MOON");
  georef.set_datum(d);

  // -------------------------------------
  // Test an unprojected coordinate system
  // -------------------------------------
  Vector2 loc = georef.pixel_to_point(Vector2(500,500));
  EXPECT_VECTOR_DOUBLE_EQ( loc, Vector2(500,500) );

  Vector2 pix = georef.point_to_pixel(Vector2(500,500));
  EXPECT_VECTOR_DOUBLE_EQ( pix, Vector2(500,500) );

  Vector2 lon_lat = georef.point_to_lonlat(Vector2(140,89));
  EXPECT_VECTOR_NEAR( lon_lat, Vector2(140,89), 1e-7 );

  loc = georef.lonlat_to_point(lon_lat);
  EXPECT_VECTOR_NEAR( loc, Vector2(140,89), 1e-7 );

  // ----------------------------------
  // Test a projected coordinate system
  // ----------------------------------
  georef.set_sinusoidal(0,0,0); // center longitude = 0

  loc = georef.pixel_to_point(Vector2(500,500));
  EXPECT_VECTOR_DOUBLE_EQ( loc, Vector2(500,500) );

  pix = georef.point_to_pixel(Vector2(500,500));
  EXPECT_VECTOR_DOUBLE_EQ( pix, Vector2(500,500) );

  // Whip up a home-brew sinusoidal projection for testing against proj.4
  double lon = 20;
  double lat = 15;
  double meters_per_degree = d.semi_major_axis() * M_PI * 2.0/360.0;
  double x = lon*cos(lat*M_PI/180) * meters_per_degree; // sinusoidal projection
  double y = lat * meters_per_degree;

  lon_lat = Vector2(lon,lat);
  loc = georef.lonlat_to_point(lon_lat);
  EXPECT_VECTOR_NEAR( loc, Vector2(x,y), 2 );

  lon_lat = georef.point_to_lonlat(loc);
  EXPECT_VECTOR_NEAR( lon_lat, Vector2(lon,lat), 1e-7 );

  double x_1 = 7000/meters_per_degree;
  double y_1 = 3000/meters_per_degree;
  double lon_1 = x_1/cos(y_1*M_PI/180);
  double lat_1 = y_1;

  lon_lat = georef.point_to_lonlat(Vector2(7000,3000));
  EXPECT_VECTOR_NEAR( lon_lat, Vector2(lon_1,lat_1), 1e-7 );

  loc = georef.lonlat_to_point(lon_lat);
  EXPECT_VECTOR_NEAR( loc, Vector2(7000,3000), 1e-7 );
}

TEST( GeoReference, UTM_to_LonLat ) {
  std::vector<Vector2> ll(5), utm(5), px(5);

  // lon, lat coordinates
  ll[0] = Vector2(170.008619281089, -43.4851542659474); // UL
  ll[1] = Vector2(170.000341300606, -43.9847965111766); // LL
  ll[2] = Vector2(170.620731735563, -43.4888289623561); // UR
  ll[3] = Vector2(170.617564317696, -43.9885355997305); // LR
  ll[4] = Vector2(170.311817924037, -43.7372482005704); // Center

  // UTM coordinates
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

  for (unsigned i = 0; i < ll.size(); ++i) {
    Vector2 check_point  = utm_georef.pixel_to_point(px[i]);
    Vector2 check_lonlat = utm_georef.pixel_to_lonlat(px[i]);

    // At this latitude, 1 degree is ~111km lat and ~81km lon
    EXPECT_VECTOR_NEAR( check_point, utm[i], 1 ); // Accurate within 1m
    EXPECT_VECTOR_NEAR( check_lonlat, ll[i], 1e-5 );
  }

}

TEST( GeoReference, LonLat_to_UTM ) {

  // Set up a UTM object. The "within-zone" part can be constant.
  
  Datum d("NAD83");
  Matrix3x3 affine;
  affine(0,0) =  1.0; // meters per pixel
  affine(1,1) = 1.0; // meters per pixel
  affine(2,2) = 1;
  affine(0,2) = 0; // Degrees
  affine(1,2) = 0;
  
  GeoReference georef(d, affine);

  // Loop through all 60 UTM zones and make sure we get a valid pixel for each one.
  const double ZONE_WIDTH = 6.0;
  double longitude = -180 + 3.012345;
  for (int zone=1; zone<=60; ++zone) {
  
    // Re-init the georef object for each UTM zone.
    std::stringstream s;
    s << "+proj=utm +zone=" << zone << " +units=m";
    std::string proj_string = s.str();
    georef.set_proj4_projection_str(proj_string );
    
    // Test lonlat to pixel conversion.
    Vector2 pixel = georef.lonlat_to_pixel(Vector2(longitude, 44.059883787682551));
    //std::cout << std::setprecision(10) << "pixel = " << pixel << std::endl;
    const double EPS = 1e-2;
    EXPECT_VECTOR_NEAR( pixel, Vector2(500988.2534,4878523.61), EPS);

    longitude += ZONE_WIDTH;
  }
}

// TODO(oalexan1): Fix this test.
#if 0
TEST( GeoReference, IOLoop ) {
  // Make a dummy image and projection transform
  ImageView<PixelRGB<float> > test_image(2,2);
  test_image(0,0) = PixelRGB<float>(1,2,3);
  test_image(0,1) = PixelRGB<float>(4,1,4);
  test_image(1,0) = PixelRGB<float>(7,7,2);
  test_image(1,1) = PixelRGB<float>(8,9,2);

  Matrix3x3 test_transform;
  test_transform(0,0) = 1.2; test_transform(0,1) = 3.2;
  test_transform(1,1) = 4.5; test_transform(1,2) = 6.3;
  test_transform(2,2) = 1; test_transform(2,1) = 0;

  {
    // Pack into dummy GeoReference
    Datum test_datum( "monkey", "dog", "cow", 7800, 6600, 3 );
    GeoReference test_georeference( test_datum, test_transform );

    // Write it to a temporary file
    UnlinkName test_filename("georeference_test.tif");
    ASSERT_NO_THROW(
      write_gdal_image(test_filename, test_image,
                             test_georeference, GdalWriteOptions()));

    // Reading back in and compare
    GeoReference retn_georeference;
    ImageView<PixelRGB<float>> retn_image;
    // Note: The function below got wiped. Read the image and georef separately.
    EXPECT_TRUE(read_georeferenced_image(retn_image, retn_georeference,
                                         test_filename));

    // Verify that the pixels are identical
    typedef ImageView<PixelRGB<float> >::iterator iterator;
    for ( iterator test = test_image.begin(), retn = retn_image.begin();
          test != test_image.end(); test++, retn++ )
      EXPECT_PIXEL_EQ( *retn, *test );

    // Check the proj4 strings and transform matrices
    EXPECT_STREQ( boost::trim_copy(retn_georeference.proj4_str()).c_str(),
                  boost::trim_copy(test_georeference.proj4_str()).c_str() );
    EXPECT_MATRIX_DOUBLE_EQ( retn_georeference.transform(),
                             test_georeference.transform() );
    // TODO: Fix this test!
/*
    std::cout << "INPUT: " << test_georeference << std::endl;
    std::cout << "OUTPUT: " << retn_georeference << std::endl;

    std::ostringstream retn_ostr, test_ostr;
    retn_ostr << retn_georeference;
    test_ostr << test_georeference;
    EXPECT_STREQ( boost::erase_all_copy(retn_ostr.str()," ").c_str(),
                  boost::erase_all_copy(test_ostr.str()," ").c_str() );*/
  }

  { // Test that spherical is handled correctly
    Datum test_datum("D_MOON");
    GeoReference test_georef( test_datum, test_transform );

    UnlinkName test_filename( "georeference_test.tif" );
    ASSERT_NO_THROW(
      write_gdal_image(test_filename, test_image,
                                 test_georef, GdalWriteOptions()));
    // Read back in and compare
    GeoReference retn_georef;
    ImageView<PixelRGB<float> > retn_image;
    // Note: The function below got wiped. Read the image and georef separately.
    EXPECT_TRUE( read_georeferenced_image( retn_image, retn_georef,
                                           test_filename ) );

    EXPECT_STREQ( boost::trim_copy(retn_georef.proj4_str()).c_str(),
                  boost::trim_copy(test_georef.proj4_str()).c_str() );
    EXPECT_MATRIX_DOUBLE_EQ( retn_georef.transform(),
                             test_georef.transform() );

  }

  { // Test that georef png can't be written
    UnlinkName fail_filename("bird.png");
    Datum test_datum( "monkey", "dog", "cow", 7800, 6600, 3 );
    GeoReference test_georeference( test_datum, test_transform );
    EXPECT_THROW( write_gdal_image(fail_filename, test_image,
                                    test_georeference, GdalWriteOptions()),
                  NoImplErr );
  }
}
#endif

TEST(GeoReference, BoundingBoxNoProj) {
  GeoReference georef;
  georef.set_pixel_interpretation(GeoReference::PixelAsPoint);
  georef.set_well_known_geogcs("WGS84");

  Matrix3x3 affine;
  affine(0,0) = 0.01; // 100 pix/degree
  affine(1,1) = -0.01; // 100 pix/degree
  affine(2,2) = 1;
  affine(0,2) = 30;   // 30 deg east
  affine(1,2) = -35;  // 35 deg south
  georef.set_transform(affine);

  BBox2i pixel_bbox(400, 300, 200, 100);
  BBox2 lonlat_bbox(georef.pixel_to_lonlat_bbox(pixel_bbox));
  BBox2 lonlat_bbox2(georef.point_to_lonlat_bbox(georef.pixel_to_point_bbox(pixel_bbox)));
  BBox2i pixel_bbox2(georef.lonlat_to_pixel_bbox(lonlat_bbox));
  BBox2i pixel_bbox3(georef.point_to_pixel_bbox(georef.lonlat_to_point_bbox(lonlat_bbox)));

  // Verify the operation of lonlat_to_pixel transforms
  EXPECT_VECTOR_NEAR(pixel_bbox.min(), pixel_bbox2.min(), 1);
  EXPECT_VECTOR_NEAR(pixel_bbox.max(), pixel_bbox2.max(), 1);
  // Verify the operation of lonlat_to_point, point_to_pixel
  EXPECT_VECTOR_NEAR(pixel_bbox.min(), pixel_bbox3.min(), 1);
  EXPECT_VECTOR_NEAR(pixel_bbox.max(), pixel_bbox3.max(), 1);
  // Verify the operation of pixel_to_point, point_to_lonlat
  EXPECT_VECTOR_NEAR(lonlat_bbox.min(), lonlat_bbox2.min(), .01);
  EXPECT_VECTOR_NEAR(lonlat_bbox.max(), lonlat_bbox2.max(), .01);
}

TEST(GeoReference, BoundingBox) {
  // TODO: I'm not familar with projections, so this is the best I'm going to do for now.
  // TODO: Test different types of projections
  Matrix3x3 affine;
  affine(0,0) = 0.01; // 100 pix/degree
  affine(1,1) = -0.01; // 100 pix/degree
  affine(2,2) = 1;
  affine(0,2) = 30;   // 30 deg east
  affine(1,2) = -35;  // 35 deg south
  GeoReference georef(Datum("WGS84"), affine, GeoReference::PixelAsPoint);
  georef.set_equirectangular(0.0, 0.0, 1.0, 0.0, 0.0);

  BBox2i pixel_bbox(400, 300, 200, 100);
  BBox2 lonlat_bbox(georef.pixel_to_lonlat_bbox(pixel_bbox));
  BBox2 lonlat_bbox2(georef.point_to_lonlat_bbox(georef.pixel_to_point_bbox(pixel_bbox)));
  BBox2i pixel_bbox2(georef.lonlat_to_pixel_bbox(lonlat_bbox));
  BBox2i pixel_bbox3(georef.point_to_pixel_bbox(georef.lonlat_to_point_bbox(lonlat_bbox)));

  // Verify the operation of lonlat_to_pixel transforms
  EXPECT_VECTOR_NEAR(pixel_bbox.min(), pixel_bbox2.min(), 1);
  EXPECT_VECTOR_NEAR(pixel_bbox.max(), pixel_bbox2.max(), 1);
  // Verify the operation of lonlat_to_point, point_to_pixel
  EXPECT_VECTOR_NEAR(pixel_bbox.min(), pixel_bbox3.min(), 1);
  EXPECT_VECTOR_NEAR(pixel_bbox.max(), pixel_bbox3.max(), 1);
  // Verify the operation of pixel_to_point, point_to_lonlat
  EXPECT_VECTOR_NEAR(lonlat_bbox.min(), lonlat_bbox2.min(), .01);
  EXPECT_VECTOR_NEAR(lonlat_bbox.max(), lonlat_bbox2.max(), .01);
}

TEST(GeoReference, NED_MATRIX) {

  // Test the lonlat_to_ned_matrix() function. Create a Cartesian
  // vector which is a combination of vectors pointing North, East,
  // and down. See if this matrix can find correctly the combination.

  Matrix3x3 test_transform;
  test_transform(0,0) = 1.2; test_transform(0,1) = 3.2;
  test_transform(1,1) = 4.5; test_transform(1,2) = 6.3;
  test_transform(2,2) = 1;   test_transform(2,1) = 0;

  Datum test_datum("D_MOON");
  GeoReference georef(test_datum, test_transform);

  Vector3 xyz(100, 232, -47);
  Vector3 G = georef.datum().cartesian_to_geodetic(xyz);

  Vector3 xyz_eps_p, xyz_eps_m, xyz_lat, xyz_lon, xyz_rad;
  double eps = 1e-6;

  // Vector pointing North
  xyz_eps_p = georef.datum().geodetic_to_cartesian(G + Vector3(0, eps, 0));
  xyz_eps_m = georef.datum().geodetic_to_cartesian(G - Vector3(0, eps, 0));
  xyz_lat   = (xyz_eps_p - xyz_eps_m)/eps;
  xyz_lat   = xyz_lat/norm_2(xyz_lat);

  // Vector pointing East
  xyz_eps_p = georef.datum().geodetic_to_cartesian(G + Vector3(eps, 0, 0));
  xyz_eps_m = georef.datum().geodetic_to_cartesian(G - Vector3(eps, 0, 0));
  xyz_lon   = (xyz_eps_p - xyz_eps_m)/eps;
  xyz_lon   = xyz_lon/norm_2(xyz_lon);

  // Vector pointing down
  xyz_rad = -xyz/norm_2(xyz);

  // The Matrix should be the same as the three vectors side by side.
  Matrix3x3 M = georef.datum().lonlat_to_ned_matrix(G);
  for (int i=0; i<3; ++i) {
    EXPECT_NEAR(M(i, 0), xyz_lat[i], eps);
    EXPECT_NEAR(M(i, 1), xyz_lon[i], eps);
    EXPECT_NEAR(M(i, 2), xyz_rad[i], eps);
  }
}

/// Loop through a bunch of pixels in an image and
///  make sure we can go from and back to the same pixel.
void georefMatchTest(const GeoReference &georef)
{

  const double MAX_PIXEL_DIFF  = 0.01;
  const double MAX_POINT_DIFF  = 0.01;
  const double MAX_DEGREE_DIFF = 0.00001;
  const double MAX_POINT_CHANGE_DIFF   = 0.1;
  const double MAX_DEGREE_CHANGE_DIFF  = 0.01;
  const int    PIXEL_SPACING = 50;

  for (int r=0; r<5; ++r) // Check a few rows
  {
    double lastX       = 0;
    double lastLon     = 0;
    double lastDiffX   = 0;
    double lastDiffLon = 0;
    for (int c=0; c<5; ++c) // Check several columns
    {
      // Make a test pixel
      Vector2 pixel1(c*PIXEL_SPACING,r*PIXEL_SPACING);
      
      // Go to a point and lonlat and back to the pixel
      Vector2 point1  = georef.pixel_to_point(pixel1);
      Vector2 lonlat1 = georef.point_to_lonlat(point1);
      Vector2 point2  = georef.lonlat_to_point(lonlat1);
      Vector2 pixel2  = georef.point_to_pixel(point2);

      // Also check equivalent shifted longitude values.
      Vector2 lonlat2 = lonlat1;
      if (lonlat1[0] > 180)
        lonlat2[0] -= 360;
      else
        lonlat2[0] += 360;

      //Vector2 point3  = georef.lonlat_to_point(lonlat2);
      Vector2 pixel3  = georef.lonlat_to_pixel(lonlat2);

      // Also check that we can go from lon to lon (at least in a limited range)
      Vector2 lonlat3 = georef.point_to_lonlat(point2);

      EXPECT_LT(fabs(pixel1 [0] - pixel2 [0]), MAX_PIXEL_DIFF ); // Did we reproject back to the same pixel?
      EXPECT_LT(fabs(pixel1 [1] - pixel2 [1]), MAX_PIXEL_DIFF );
      EXPECT_LT(fabs(point1 [0] - point2 [0]), MAX_POINT_DIFF ); // Did we go back and forth through the same point?
      EXPECT_LT(fabs(point1 [1] - point2 [1]), MAX_POINT_DIFF );
      EXPECT_LT(fabs(pixel1 [0] - pixel3 [0]), MAX_PIXEL_DIFF ); // Can we handle longitudes offset by 360 degrees?
      EXPECT_LT(fabs(pixel1 [1] - pixel3 [1]), MAX_PIXEL_DIFF );
      EXPECT_LT(fabs(lonlat1[0] - lonlat3[0]), MAX_DEGREE_DIFF); // Can we go back and forth to the same longitude?
      EXPECT_LT(fabs(lonlat1[1] - lonlat3[1]), MAX_DEGREE_DIFF);

      // How much did point and lon location change since last time?
      double diffX   = point2[0] - lastX;
      double diffLon = lonlat1[0] - lastLon;

      // Make sure that the rate of change horizontally does not vary much as we move through pixels.
      if (c > 1) {
        EXPECT_LT(fabs(diffX   - lastDiffX  ), MAX_POINT_CHANGE_DIFF);
        EXPECT_LT(fabs(diffLon - lastDiffLon), MAX_DEGREE_CHANGE_DIFF);
      }

      // Update records
      lastDiffX   = diffX;
      lastDiffLon = diffLon;
      lastX       = point2[0];
      lastLon     = lonlat1[0];

    }
  }

}

/// Check back-and-forth conversions for several equirectangular test cases
/// - Our code can't handle all test cases, but it should be able to handle
///   anything that is likely to be seen.
TEST( GeoReference, eqcReverseTest) {

  //std::cout << "Default init.\n";
  GeoReference georef;
  Matrix3x3 affine;

  // Test #1

  georef.set_well_known_geogcs("D_MOON");

  // Offset these points near where proj4 can mess with the conversions
  double pi     = 3.141592653589793238462643383279;
  double radius = georef.datum().semi_major_axis();
  double bound  = radius * pi;

  // Set up a transform that needs to be handled in the 0-360 range.
  // - This image could wrap around the 360 line, but that would not be a valid georeference.
  affine(0,0) =  1.0; // 1 meters per pixel
  affine(1,1) = -1.0; // 1 meters per pixel
  affine(2,2) = 1;
  affine(0,2) = 2*bound - 5000; // This is short of 360 degrees in proj4 space
  affine(1,2) = 50000;  // Some value in the northern hemisphere
  georef.set_transform(affine);

  georef.set_well_known_geogcs("D_MOON");
  georef.set_equirectangular(0, 0, 0, 0, 0); // Default EQC parameters
  //std::cout << georef << std::endl;
  //std::cout << "Finished initializing georef!\n";

  georefMatchTest(georef); // Run a set of tests on the georef

  // Test #2

  // An arbitrary weird set of parameters
  affine(0,0) =  10.0; // 10 meters per pixel
  affine(1,1) = -8.0; // 8 meters per pixel
  affine(2,2) = 1;
  affine(0,2) = 600; //
  affine(1,2) = 50000;  // Some value in the northern hemisphere
  georef.set_transform(affine);

  // center_latitude, center_longitude, latitude_of_true_scale, false_easting, false_northing
  georef.set_equirectangular(-18, 13, -45, -500, -4000);
  georef.set_well_known_geogcs("D_MOON");
  //std::cout << georef << std::endl;
  //std::cout << "Finished initializing georef!\n";

  georefMatchTest(georef); // Run a set of tests on the georef

  // Test #3

  // The center lon and the projection offset put this near the -180 line
  // - This georef needs to be handled in the -180 to 180 range.
  affine(0,0) =  -3.0; // 3 meters per pixel
  affine(1,1) = 11.0; // 11 meters per pixel
  affine(2,2) = 1;
  affine(0,2) = -bound/2;
  affine(1,2) = -50000;  // Some value in the southern hemisphere
  georef.set_transform(affine);

  // center_latitude, center_longitude, latitude_of_true_scale, false_easting, double false_northing
  georef.set_equirectangular(-45, -15, 13, 14, 8.5);
  georef.set_well_known_geogcs("D_MOON");
  //std::cout << georef << std::endl;
  //std::cout << "Finished initializing georef!\n";

  georefMatchTest(georef); // Run a set of tests on the georef
}

TEST( GeoReference, orthoTest) {

  Matrix3x3 affine;
  Datum d;
  std::string proj_str = "+proj=ortho +lat_0=37 +lon_0=350 +x_0=0 +y_0=0 +a=606000 +b=606000 +units=m +no_defs";
  d.set_datum_from_proj_str(proj_str);

  affine(0,0) =  1200; // meters per pixel
  affine(1,1) = -1200; // meters per pixel
  affine(2,2) = 1;
  affine(0,2) = -606000.000; // Degrees
  affine(1,2) = 606000.000;

  GeoReference georef(d, affine); 
  georef.set_proj4_projection_str(proj_str);
  //std::cout << georef << std::endl;

  Vector2 lonlat = georef.pixel_to_lonlat(Vector2(400,600));
  EXPECT_VECTOR_NEAR(lonlat, Vector2(-23.2265,25.2555), 10e-4);
}


TEST( GeoReference, albersTestNegLon) {
  Matrix3x3 affine;
  Datum d;
  std::string proj_str = "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs";
  d.set_datum_from_proj_str(proj_str);

  affine(0,0) =  50; // meters per pixel
  affine(1,1) = -50; // meters per pixel
  affine(2,2) = 1;
  affine(0,2) = -379823.0;
  affine(1,2) = 2025169.0;

  GeoReference georef(d, affine); 
  georef.set_proj4_projection_str(proj_str); 
  //std::cout << georef << std::endl;

  Vector2 lonlat = georef.pixel_to_lonlat(Vector2(30,30));
  EXPECT_VECTOR_NEAR(lonlat, Vector2(-162.981, 67.9422), 10e-3);
  Vector2 pix = georef.lonlat_to_pixel(lonlat);
  EXPECT_VECTOR_NEAR(pix, Vector2(30,30), 10e-2);
}

TEST( GeoReference, albersTestHighLon) {
  Matrix3x3 affine;
  Datum d;
  std::string proj_str = "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=200 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs";
  d.set_datum_from_proj_str(proj_str);

  affine(0,0) =  50; // meters per pixel
  affine(1,1) = -50; // meters per pixel
  affine(2,2) = 1;
  affine(0,2) = -379823.0;
  affine(1,2) = 2025169.0;

  GeoReference georef(d, affine); 
  georef.set_proj4_projection_str(proj_str); 
  //std::cout << georef << std::endl;

  Vector2 lonlat = georef.pixel_to_lonlat(Vector2(30,30));
  EXPECT_VECTOR_NEAR(lonlat, Vector2(191.019, 67.9422), 10e-3);
  Vector2 pix = georef.lonlat_to_pixel(lonlat);
  EXPECT_VECTOR_NEAR(pix, Vector2(30,30), 10e-2);
}

//TEST( GeoReference, read_strings) {
//
//  std::string file = "/home/smcmich1/data/icebridge_bulk/AN_2011_10_18/ortho/DMS_1281710_00355_20111018_14300236.tif";
//  boost::shared_ptr<vw::DiskImageResource> resource(vw::DiskImageResource::open(file));
//  std::map<std::string, std::string> entries;
//  read_header_strings(*resource.get(), entries);
//  std::map<std::string, std::string>::const_iterator iter;
//  for (iter = entries.begin(); iter!=entries.end(); ++iter) {
//    std::cout << iter->first << ", " << iter->second << std::endl;
//  }
//  EXPECT_TRUE(false);
//}



