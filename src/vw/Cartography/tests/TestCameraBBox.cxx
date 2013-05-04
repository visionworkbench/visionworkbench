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


// TestCameraBBox.h
#include <gtest/gtest_VW.h>
#include <test/Helpers.h>

#include <vw/Cartography/CameraBBox.h>
#include <vw/Camera/PinholeModel.h>

// Must have protobuf to be able to read camera
#if defined(VW_HAVE_PKG_PROTOBUF) && VW_HAVE_PKG_PROTOBUF==1 && defined(VW_HAVE_PKG_CAMERA) && VW_HAVE_PKG_CAMERA

using namespace vw;
using namespace vw::cartography;
using namespace vw::test;
using namespace vw::camera;

class CameraBBoxTest :  public ::testing::Test {
protected:
  virtual void SetUp() {
    apollo_camera = boost::shared_ptr<CameraModel>(new PinholeModel(TEST_SRCDIR"/apollo.pinhole"));
    moon_georef.set_well_known_geogcs("D_MOON");
    DEM.set_size(20,20); // DEM covering lat {10,-10} long {80,100}
    for ( int32 i = 0; i < DEM.cols(); i++ )
      for ( int32 j = 0; j <DEM.rows(); j++ )
        DEM(i,j) = 1000 - 10*(pow(DEM.cols()/2.0-i,2.0)+pow(DEM.rows()/2.0-j,2.0));
  }

  boost::shared_ptr<CameraModel> apollo_camera;
  GeoReference moon_georef;
  ImageView<float> DEM;
};

TEST_F( CameraBBoxTest, GeospatialIntersectDatum ) {
  for ( unsigned i = 0; i < 10; i++ ) {
    Vector2 input_image( rand()%4096, rand()%4096 );
    bool did_intersect;
    Vector2 lonlat =
      geospatial_intersect( moon_georef,
                            apollo_camera->camera_center(input_image),
                            apollo_camera->pixel_to_vector(input_image),
                            did_intersect );
    ASSERT_TRUE( did_intersect );
    double radius = moon_georef.datum().radius( lonlat[0], lonlat[1] );
    EXPECT_NEAR( radius, 1737400, 1e-3 );
    Vector3 llr( lonlat[0], lonlat[1], 0 );
    Vector3 ecef = moon_georef.datum().geodetic_to_cartesian(llr);
    Vector3 llr2 = moon_georef.datum().cartesian_to_geodetic(ecef);
    EXPECT_VECTOR_NEAR( llr2, llr, 1e-4 );

    Vector2 retrn_image = apollo_camera->point_to_pixel( ecef );
    EXPECT_VECTOR_NEAR( retrn_image, input_image, 1e-3 );
  }
}

TEST_F( CameraBBoxTest, CameraBBoxDatum ) {
  float scale; // degrees per pixel
  BBox2 image_bbox = camera_bbox( moon_georef, apollo_camera, 4096, 4096, scale );
  EXPECT_VECTOR_NEAR( image_bbox.min(), Vector2(86,-1), 2 );
  EXPECT_VECTOR_NEAR( image_bbox.max(), Vector2(95,7), 2 );
  EXPECT_NEAR( scale, (95-86.)/sqrt(4096*4096*2), 1e-3 ); // Cam is rotated
}

TEST_F( CameraBBoxTest, CameraBBoxDEM ) {
  Matrix<double> geotrans = vw::math::identity_matrix<3>();
  geotrans(0,2) = 80;
  geotrans(1,1) = -1;
  geotrans(1,2) = 10;
  moon_georef.set_transform(geotrans);

  BBox2 image_bbox = camera_bbox( DEM, moon_georef, apollo_camera, 4096, 4096 );
  EXPECT_VECTOR_NEAR( image_bbox.min(), Vector2(87,0), 2 );
  EXPECT_VECTOR_NEAR( image_bbox.max(), Vector2(94,6), 2 );
}

TEST_F( CameraBBoxTest, CameraPixelToXYZ) {
  Matrix<double> geotrans = vw::math::identity_matrix<3>();
  geotrans(0,2) = 80;
  geotrans(1,1) = -1;
  geotrans(1,2) = 10;
  moon_georef.set_transform(geotrans);

  Vector2 input_pixel(50,10);
  bool treat_nodata_as_zero = false;
  bool has_intersection;
  Vector3 xyz = camera_pixel_to_dem_xyz(apollo_camera->camera_center(input_pixel),
                                        apollo_camera->pixel_to_vector(input_pixel),
                                        DEM, moon_georef,
                                        treat_nodata_as_zero,
                                        has_intersection
                                        );

  // Verification
  Vector2 output_pixel = apollo_camera->point_to_pixel(xyz);

  EXPECT_EQ(has_intersection, 1);
  EXPECT_VECTOR_NEAR(input_pixel, output_pixel, 1e-8);
}

#endif
