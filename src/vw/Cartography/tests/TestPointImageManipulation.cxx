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

#include <vw/Cartography/PointImageManipulation.h>
#include <boost/math/special_functions/fpclassify.hpp>

using namespace vw;
using namespace vw::cartography;
using namespace vw::test;

// Spherical transform that is ignorant of flattened datums.
TEST( PointImageManipulation, XYZ_to_LonLat_Estimate ) {
  // Test a full forward and reverse transformation
  Vector3 xyz(-2197110.000000, 1741355.875000, 1898886.875000);
  Vector3 lon_lat_alt = xyz_to_lon_lat_radius_estimate(xyz);
  Vector3 xyz2        = lon_lat_radius_to_xyz_estimate(lon_lat_alt);

  EXPECT_VECTOR_NEAR( xyz, xyz2, 1e-2 );

  // Test to see if things still work for West positive coordinate
  // systems.
  lon_lat_alt = xyz_to_lon_lat_radius_estimate(xyz,false);
  xyz2        = lon_lat_radius_to_xyz_estimate(lon_lat_alt,false);

  EXPECT_VECTOR_NEAR( xyz, xyz2, 1e-2 );

  // See if it still works if using 0-360 range
  xyz[1] = -xyz[1];
  lon_lat_alt = xyz_to_lon_lat_radius_estimate(xyz,true,false);
  xyz2        = lon_lat_radius_to_xyz_estimate(lon_lat_alt);

  EXPECT_VECTOR_NEAR( xyz, xyz2, 1e-2 );
}

// These are the more general operators which can actually handle
// squished datums.
TEST( PointImageManipulation, GeodeticCartesian ) {
  ImageView<Vector3> geodetic(2,2);
  geodetic(0,0) = Vector3( 90, 10, 100 );
  geodetic(0,1) = Vector3(); // An actual valid measurement.
  geodetic(1,0) = Vector3( 10, 89, 9  );
  geodetic(1,1) = Vector3(0, 0, std::numeric_limits<double>::quiet_NaN() ); // invalid measure.

  Datum moon("D_MOON"), earth("WGS84");
  ImageView<Vector3> result_moon  = cartesian_to_geodetic(geodetic_to_cartesian(geodetic, moon),moon);
  ImageView<Vector3> result_earth = cartesian_to_geodetic(geodetic_to_cartesian(geodetic,earth),earth);
  EXPECT_RANGE_NEAR( geodetic.begin(), geodetic.begin()+3, result_moon.begin(), result_moon.begin()+3, 1e-9 );
  EXPECT_TRUE( boost::math::isnan(result_moon(1,1).z()) );
  EXPECT_RANGE_NEAR( geodetic.begin(), geodetic.begin()+3, result_earth.begin(), result_earth.begin()+3, 3e-9 ); // earth is bigger, so bigger error too
  EXPECT_TRUE( boost::math::isnan(result_earth(1,1).z()) );
}

TEST( PointImageManipulation, CartesianGeodetic ) {
  ImageView<Vector3> cartesian(2,2);
  cartesian(0,0) = Vector3(1737892,80,5320);
  cartesian(1,0) = Vector3(50,-190,-80);
  cartesian(0,1) = Vector3(33e8, -1e5, 0);
  cartesian(1,1) = Vector3(); // Invalid in for cartesian

  Datum moon("D_MOON"), earth("WGS84");
  ImageView<Vector3> result_moon  = geodetic_to_cartesian(cartesian_to_geodetic(cartesian, moon),moon);
  ImageView<Vector3> result_earth = geodetic_to_cartesian(cartesian_to_geodetic(cartesian,earth),earth);
  EXPECT_SEQ_NEAR( cartesian, result_moon,  1e-6 );
  EXPECT_SEQ_NEAR( cartesian, result_earth, 1e-6 );
}
