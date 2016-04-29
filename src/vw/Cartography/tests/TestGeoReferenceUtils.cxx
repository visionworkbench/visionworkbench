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

#include <vw/Cartography/GeoReferenceUtils.h>

#include <iomanip>

using namespace vw;
using namespace vw::cartography;
using namespace vw::test;

TEST(GeoReferenceUtils, CropAndResample) {
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

  // Test with a projection that doesn't use degrees.
  input_pp.set_equirectangular( 40, 40, 50, 0, 0 );
  crop_pp =  crop( input_pp, 200, 400 );
  resample_pp = resample( input_pp, 0.5, 2 );
  EXPECT_VECTOR_NEAR( input_pp.pixel_to_lonlat(Vector2(200, 400)),
                      crop_pp.pixel_to_lonlat(Vector2()), 1e-7 );
  EXPECT_VECTOR_NEAR( input_pp.pixel_to_lonlat(Vector2()),
                      crop_pp.pixel_to_lonlat(Vector2(-200,-400)), 1e-7 );
  EXPECT_VECTOR_NEAR( input_pp.pixel_to_lonlat(Vector2()),
                      resample_pp.pixel_to_lonlat(Vector2()), 1e-7 );
  EXPECT_VECTOR_NEAR( input_pp.pixel_to_lonlat(Vector2(100,100)),
                      resample_pp.pixel_to_lonlat(Vector2(50,200)), 1e-7 );
}





