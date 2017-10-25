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
#include <vw/FileIO/DiskImageView.h>

#if defined(VW_HAVE_PKG_CAMERA) && VW_HAVE_PKG_CAMERA

using namespace vw;
using namespace vw::cartography;
using namespace vw::test;
using namespace vw::camera;

// Load a pinhole camera and an associated low res DEM in Antarctica.
class CameraBBoxTest :  public ::testing::Test {
protected:
  virtual void SetUp() {
    pinhole_camera = boost::shared_ptr<CameraModel>(new PinholeModel("pinhole_AN.tsai"));
    earth_georef.set_well_known_geogcs("WGS84");
    DEM = create_mask(DiskImageView<short>("tinyDemAN.tif"), -32768);
  }

  boost::shared_ptr<CameraModel> pinhole_camera;
  GeoReference earth_georef;
  ImageView< PixelMask<float> > DEM;
};

// Test datum intersection
TEST_F( CameraBBoxTest, GeospatialIntersectDatum ) {
  for ( unsigned i = 0; i < 10; i++ ) {
    Vector2 input_image( rand()%5616, rand()%3744 );
    bool did_intersect;
    Vector2 lonlat =
      geospatial_intersect( earth_georef,
                            pinhole_camera->camera_center(input_image),
                            pinhole_camera->pixel_to_vector(input_image),
                            did_intersect );
    ASSERT_TRUE( did_intersect );
    double radius = earth_georef.datum().radius( lonlat[0], lonlat[1] );
    EXPECT_NEAR( radius, 6360070, 10);
    Vector3 llr( lonlat[0], lonlat[1], 0 );
    Vector3 ecef = earth_georef.datum().geodetic_to_cartesian(llr);
    Vector3 llr2 = earth_georef.datum().cartesian_to_geodetic(ecef);
    EXPECT_VECTOR_NEAR( llr2, llr, 1e-4 );

    Vector2 retrn_image = pinhole_camera->point_to_pixel( ecef );
    EXPECT_VECTOR_NEAR( retrn_image, input_image, 1e-3 );
  }
}

// Test the datum version of camera_bbox 
TEST_F( CameraBBoxTest, CameraBBoxDatum ) {
  float scale; // degrees per pixel
  BBox2 image_bbox = camera_bbox( earth_georef, pinhole_camera, 5616, 3744, scale );
  EXPECT_VECTOR_NEAR( image_bbox.min(), Vector2(-64.3442, -66.7563), 0.0004 );
  EXPECT_VECTOR_NEAR( image_bbox.max(), Vector2(-64.2840, -66.7254), 0.0004 );
  EXPECT_NEAR( scale, 2.8e-6, 1e-7 ); // Cam is rotated
}

// Test the dem version of camera_bbox 
TEST_F( CameraBBoxTest, CameraBBoxDEM ) {
  GeoReference dem_georef;
  read_georeference(dem_georef, "tinyDemAN.tif");
  BBox2 image_bbox = camera_bbox( DEM, dem_georef, dem_georef, pinhole_camera, 5616, 3744 );
  EXPECT_VECTOR_NEAR( image_bbox.min(), Vector2(-2.30944e+06,1.10892e+06), 5 );
  EXPECT_VECTOR_NEAR( image_bbox.max(), Vector2(-2.30635e+06,1.11089e+06), 5 );
}

// Try projecting a pixel to and from the DEM using the bbox tool.
TEST_F( CameraBBoxTest, CameraPixelToXYZ) {
  GeoReference dem_georef;
  read_georeference(dem_georef, "tinyDemAN.tif");
  Vector2 input_pixel(50,10);
  bool treat_nodata_as_zero = false;
  bool has_intersection;
  Vector3 xyz = camera_pixel_to_dem_xyz(pinhole_camera->camera_center(input_pixel),
                                        pinhole_camera->pixel_to_vector(input_pixel),
                                        DEM, dem_georef,
                                        treat_nodata_as_zero,
                                        has_intersection
                                        );

  // Verification
  Vector2 output_pixel = pinhole_camera->point_to_pixel(xyz);

  EXPECT_EQ(has_intersection, 1);
  EXPECT_VECTOR_NEAR(input_pixel, output_pixel, 1e-4);
}


/* // Make sure camera_bbox works properly in a non-toy case.
TEST( CameraBBox, CameraBBoxDEM11 ) {
  
  std::string demPath = "/home/smcmich1/data/icebridge_bulk/reference_dems/ramp200dem_wgs_v2.tif";
  GeoReference georef;
  if (!read_georeference(georef, demPath))
    vw_throw( ArgumentErr() << "Missing georef.\n");


  float dem_nodata_val = -std::numeric_limits<float>::max(); 
  vw::read_nodata_val(demPath, dem_nodata_val);
  ImageViewRef< PixelMask<double> > dem = create_mask
    (channel_cast<double>(DiskImageView<float>(demPath)), dem_nodata_val);

  boost::shared_ptr<vw::camera::CameraModel> camPtr(new PinholeModel("/home/smcmich1/data/icebridge_bulk/AN_20111114/nav_camera/DMS_20111114_164409_04413.tsai_M1.tsai"));

  float res;
  BBox2 image_bbox = camera_bbox( dem, georef, georef, camPtr, 5616, 3744, res, false );
  std::cout << image_bbox << std::endl;
  EXPECT_VECTOR_NEAR( image_bbox.min(), Vector2(87,0), 2 );
}
*/

#endif
