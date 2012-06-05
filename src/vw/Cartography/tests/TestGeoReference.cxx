// __BEGIN_LICENSE__
//  Copyright (c) 2006-2012, United States Government as represented by the
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

using namespace vw;
using namespace vw::cartography;
using namespace vw::test;


#if defined(VW_HAVE_PKG_PROTOBUF) && VW_HAVE_PKG_PROTOBUF==1
TEST( GeoReference, GeoReferenceDesc) {
  GeoReference georef;
  georef.set_pixel_interpretation(GeoReferenceBase::PixelAsPoint);
  georef.set_well_known_geogcs("WGS84");

  Matrix3x3 affine;
  affine(0,0) = 0.01; // 100 pix/degree
  affine(1,1) = -0.01; // 100 pix/degree
  affine(2,2) = 1;
  affine(0,2) = 30;   // 30 deg east
  affine(1,2) = -35;  // 35 deg south
  georef.set_transform(affine);

  GeoReference desc = georef.build_desc();

  GeoReference georef2(desc);

  EXPECT_EQ(georef.build_desc().DebugString(), georef2.build_desc().DebugString());

  Vector2 pix = georef.point_to_pixel(Vector2(30,-35));
  EXPECT_VECTOR_DOUBLE_EQ( pix, Vector2() );
  pix = georef2.point_to_pixel(Vector2(30,-35));
  EXPECT_VECTOR_DOUBLE_EQ( pix, Vector2() );


  georef.set_pixel_interpretation(GeoReference::PixelAsArea);
  georef2.set_pixel_interpretation(GeoReference::PixelAsArea);
  pix = georef.point_to_pixel(Vector2(30,-35));
  EXPECT_VECTOR_DOUBLE_EQ( pix, Vector2(-0.5,-0.5) );
  pix = georef2.point_to_pixel(Vector2(30,-35));
  EXPECT_VECTOR_DOUBLE_EQ( pix, Vector2(-0.5,-0.5) );
}
#endif

TEST( GeoReference, BasicGeographic ) {
  GeoReference georef;
  georef.set_pixel_interpretation(GeoReferenceBase::PixelAsPoint);
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
  georef.set_pixel_interpretation(GeoReferenceBase::PixelAsPoint);

  Datum d = georef.datum();
  d.set_semi_minor_axis(d.semi_major_axis());
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

  // Whip up a home-brew sinusoidal projection for testing against
  // proj.4
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

  for (unsigned i = 0; i < ll.size(); ++i) {
    Vector2 check_point  = utm_georef.pixel_to_point(px[i]);
    Vector2 check_lonlat = utm_georef.pixel_to_lonlat(px[i]);

    // At this latitude, 1 degree is ~111km lat and ~81km lon
    EXPECT_VECTOR_NEAR( check_point, utm[i], 1 ); // Accurate within 1m
    EXPECT_VECTOR_NEAR( check_lonlat, ll[i], 1e-5 );
  }
}

TEST( GeoReference, IOLoop ) {
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
    Datum test_datum( "monkey", "dog", "cow", 7800, 6600, 3 );
    GeoReference test_georeference( test_datum, test_transform );

    UnlinkName test_filename( "georeference_test.tif" );
    ASSERT_NO_THROW(
      write_georeferenced_image( test_filename, test_image,
                             test_georeference ));

    // Reading back in and compare
    GeoReference retn_georeference;
    ImageView<PixelRGB<float> > retn_image;
    EXPECT_TRUE( read_georeferenced_image( retn_image, retn_georeference,
                                           test_filename ) );

    typedef ImageView<PixelRGB<float> >::iterator iterator;
    for ( iterator test = test_image.begin(), retn = retn_image.begin();
          test != test_image.end(); test++, retn++ )
      EXPECT_PIXEL_EQ( *retn, *test );

    EXPECT_STREQ( boost::trim_copy(retn_georeference.proj4_str()).c_str(),
                  boost::trim_copy(test_georeference.proj4_str()).c_str() );
    EXPECT_MATRIX_DOUBLE_EQ( retn_georeference.transform(),
                             test_georeference.transform() );

    EXPECT_STREQ( retn_georeference.gml_str().c_str(),
                  test_georeference.gml_str().c_str() );

    std::ostringstream retn_ostr, test_ostr;
    retn_ostr << retn_georeference;
    test_ostr << test_georeference;
    EXPECT_STREQ( boost::erase_all_copy(retn_ostr.str()," ").c_str(),
                  boost::erase_all_copy(test_ostr.str()," ").c_str() );
  }

  { // Test that spherical is handled correctly
    Datum test_datum("D_MOON");
    GeoReference test_georef( test_datum, test_transform );

    UnlinkName test_filename( "georeference_test.tif" );
    ASSERT_NO_THROW(
      write_georeferenced_image( test_filename, test_image,
                                 test_georef ));
    // Read back in and compare
    GeoReference retn_georef;
    ImageView<PixelRGB<float> > retn_image;
    EXPECT_TRUE( read_georeferenced_image( retn_image, retn_georef,
                                           test_filename ) );

    EXPECT_STREQ( boost::trim_copy(retn_georef.proj4_str()).c_str(),
                  boost::trim_copy(test_georef.proj4_str()).c_str() );
    EXPECT_MATRIX_DOUBLE_EQ( retn_georef.transform(),
                             test_georef.transform() );

    EXPECT_STREQ( retn_georef.gml_str().c_str(),
                  test_georef.gml_str().c_str() );
  }

  { // Test that georef png can't be written
    UnlinkName fail_filename("bird.png");
    Datum test_datum( "monkey", "dog", "cow", 7800, 6600, 3 );
    GeoReference test_georeference( test_datum, test_transform );
    EXPECT_THROW( write_georeferenced_image( fail_filename, test_image,
                                             test_georeference ),
                  NoImplErr );
  }
}

TEST(GeoReference, BoundingBoxNoProj) {
  GeoReference georef;
  georef.set_pixel_interpretation(GeoReferenceBase::PixelAsPoint);
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
  // TODO: I'm not familar with projections, so this is the best I'm going to
  // do for now.
  // TODO: Test different types of projections
  Matrix3x3 affine;
  affine(0,0) = 0.01; // 100 pix/degree
  affine(1,1) = -0.01; // 100 pix/degree
  affine(2,2) = 1;
  affine(0,2) = 30;   // 30 deg east
  affine(1,2) = -35;  // 35 deg south
  GeoReference georef(Datum("WGS84"), affine, GeoReferenceBase::PixelAsPoint);
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

TEST(GeoReference, CropAndResample) {
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
}
