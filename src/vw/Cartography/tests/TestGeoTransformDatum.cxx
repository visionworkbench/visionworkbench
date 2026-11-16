// __BEGIN_LICENSE__
//  Copyright (c) 2006-2026, United States Government as represented by the
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

// GeoTransform passes coordinates between georefs through longitude and
// latitude, with no datum conversion, so it must refuse georefs with different
// datums. Name variants of one non-Earth datum (same shape and prime meridian)
// count as the same datum.
#include <gtest/gtest_VW.h>
#include <test/Helpers.h>

#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoTransform.h>

using namespace vw;

namespace {

std::string wgs84    = "+proj=longlat +datum=WGS84 +no_defs";
std::string wgs84_ab = "+proj=longlat +a=6378137 +rf=298.257223563 +no_defs";
std::string nad83    = "+proj=longlat +datum=NAD83 +no_defs";
std::string wgs72    = "+proj=longlat +ellps=WGS72 +no_defs";
std::string moon     = "+proj=longlat +a=1737400 +b=1737400 +no_defs";
std::string moon2000 = "GEOGCS[\"Moon 2000\",DATUM[\"D_Moon_2000\","
  "SPHEROID[\"Moon_2000_IAU_IAG\",1737400,0]],PRIMEM[\"Greenwich\",0],"
  "UNIT[\"Decimal_Degree\",0.0174532925199433]]";
std::string mars_sph = "+proj=longlat +a=3396190 +b=3396190 +no_defs";
std::string mars_ell = "+proj=longlat +a=3396190 +b=3376200 +no_defs";
std::string utm      = "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs";

GeoReference make_georef(std::string const& srs) {
  OGRSpatialReference sr;
  if (sr.SetFromUserInput(srs.c_str()) != OGRERR_NONE)
    vw_throw(ArgumentErr() << "Cannot parse: " << srs);
  char *wkt = NULL;
  sr.exportToWkt(&wkt);
  GeoReference g;
  g.set_wkt(wkt);
  CPLFree(wkt);
  Matrix3x3 affine;
  affine(0,0) =  0.001;
  affine(1,1) = -0.001;
  affine(2,2) = 1;
  affine(0,2) = 10.0;
  affine(1,2) = 20.0;
  g.set_transform(affine);
  return g;
}

bool same(std::string const& srs1, std::string const& srs2) {
  return same_datum(make_georef(srs1).datum(), make_georef(srs2).datum());
}

Vector2 forward_pixel(std::string const& srs1, std::string const& srs2) {
  GeoReference g1 = make_georef(srs1), g2 = make_georef(srs2);
  GeoTransform t(g1, g2);
  return t.forward(Vector2(100, 200));
}

}

TEST(GeoTransformDatum, SameDatum) {
  EXPECT_TRUE(same(wgs84, wgs84));
  EXPECT_TRUE(same(wgs84, "EPSG:4326"));
  EXPECT_TRUE(same(wgs84, wgs84_ab));
  EXPECT_TRUE(same(wgs84, utm));
  EXPECT_TRUE(same(moon, moon2000));

  EXPECT_FALSE(same(wgs84, nad83));
  EXPECT_FALSE(same(wgs84, wgs72));
  EXPECT_FALSE(same(nad83, wgs72));
  EXPECT_FALSE(same(wgs84, moon));
  EXPECT_FALSE(same(mars_sph, mars_ell));
  EXPECT_FALSE(same(moon, mars_sph));
}

TEST(GeoTransformDatum, ThrowOnDifferentDatums) {
  EXPECT_THROW(forward_pixel(wgs84, nad83),       ArgumentErr);
  EXPECT_THROW(forward_pixel(wgs84, wgs72),       ArgumentErr);
  EXPECT_THROW(forward_pixel(nad83, wgs84),       ArgumentErr);
  EXPECT_THROW(forward_pixel(wgs84, moon),        ArgumentErr);
  EXPECT_THROW(forward_pixel(mars_sph, mars_ell), ArgumentErr);
}

TEST(GeoTransformDatum, PassThroughOnSameDatum) {
  // Both georefs have the same grid, so a pixel maps to itself.
  Vector2 same_pix(100, 199.99999999999636);
  double tol = 1e-8;
  EXPECT_VECTOR_NEAR(forward_pixel(wgs84, wgs84),       same_pix, tol);
  EXPECT_VECTOR_NEAR(forward_pixel(wgs84, "EPSG:4326"), same_pix, tol);
  EXPECT_VECTOR_NEAR(forward_pixel(wgs84, wgs84_ab),    same_pix, tol);
  EXPECT_VECTOR_NEAR(forward_pixel(moon, moon2000),     same_pix, tol);

  // Geographic to projected, same datum.
  EXPECT_VECTOR_NEAR(forward_pixel(wgs84, utm),
                     Vector2(615258603.52521682, -2189649233.8339739), 1e-4);
}
