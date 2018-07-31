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


// TestDatum.h
#include <test/Helpers.h>

#include <vw/Cartography/Datum.h>
#include <vw/Core/Stopwatch.h>
#include <boost/assign/std/vector.hpp>

using namespace vw;
using namespace vw::cartography;
using namespace vw::test;

using namespace boost::assign;

vw::Vector3 vermille_2011_cart_to_geodetic( Datum const& d, Vector3 const& cart ) {
  const double a2 = d.semi_major_axis() * d.semi_major_axis();
  const double b2 = d.semi_minor_axis() * d.semi_minor_axis();
  const double e2 = 1 - b2 / a2;
  const double e4 = e2 * e2;

  double xy_dist = sqrt( cart[0] * cart[0] + cart[1] * cart[1] );
  double p = ( cart[0] * cart[0] + cart[1] * cart[1] ) / a2;
  double q = ( 1 - e2 ) * cart[2] * cart[2] / a2;
  double r = ( p + q - e4 ) / 6.0;
  double r3 = r * r * r;

  Vector3 llh;

  double evolute = 8 * r3 + e4 * p * q;
  double u = std::numeric_limits<double>::quiet_NaN();
  if ( evolute > 0 ) {
    // outside the evolute
    double right_inside_pow = sqrt(e4 * p * q);
    double sqrt_evolute = sqrt( evolute );
    u = r + 0.5 * pow(sqrt_evolute + right_inside_pow,2.0/3.0) +
      0.5 * pow(sqrt_evolute - right_inside_pow,2.0/3.0);
  } else if ( fabs(cart[2]) < std::numeric_limits<double>::epsilon() ) {
    // On the equator plane
    llh[1] = 0;
    llh[2] = norm_2( cart ) - d.semi_major_axis();
  } else if ( evolute < 0 and fabs(q) > std::numeric_limits<double>::epsilon() ) {
    // On or inside the evolute
    double atan_result = atan2( sqrt( e4 * p * q ), sqrt( -evolute ) + sqrt(-8 * r3) );
    u = -4 * r * sin( 2.0 / 3.0 * atan_result ) *
      cos( M_PI / 6.0 + 2.0 / 3.0 * atan_result );
  } else if ( fabs(q) < std::numeric_limits<double>::epsilon() and p <= e4 ) {
    // In the singular disc
    llh[2] = -d.semi_major_axis() * sqrt(1 - e2) * sqrt(e2 - p) / sqrt(e2);
    llh[1] = 2 * atan2( sqrt(e4 - p), sqrt(e2*(e2 - p)) + sqrt(1-e2) * sqrt(p) );
  } else {
    // Near the cusps of the evolute
    double inside_pow = sqrt(evolute) + sqrt(e4 * p * q);
    u = r + 0.5 * pow(inside_pow,2.0/3.0) +
      2 * r * r * pow(inside_pow,-2.0/3.0);
  }

  if (!std::isnan(u) ) {
    double v = sqrt( u * u + e4 * q );
    double u_v = u + v;
    double w = e2 * ( u_v - q ) / ( 2 * v );
    double k = u_v / ( w + sqrt( w * w + u_v ) );
    double D = k * xy_dist / ( k + e2 );
    double dist_2 = D * D + cart[2] * cart[2];
    llh[2] = ( k + e2 - 1 ) * sqrt( dist_2 ) / k;
    llh[1] = 2 * atan2( cart[2], sqrt( dist_2 ) + D );
  }

  if ( xy_dist + cart[0] > ( sqrt(2) - 1 ) * cart[1] ) {
    // Longitude is between -135 and 135
    llh[0] = 360.0 * atan2( cart[1], xy_dist + cart[0] ) / M_PI;
  } else if ( xy_dist + cart[1] < ( sqrt(2) + 1 ) * cart[0] ) {
    // Longitude is between -225 and 45
    llh[0] = - 90.0 + 360.0 * atan2( cart[0], xy_dist - cart[1] ) / M_PI;
  } else {
    // Longitude is between -45 and 225
    llh[0] = 90.0 - 360.0 * atan2( cart[0], xy_dist + cart[1] ) / M_PI;
  }
  llh[1] *= 180.0 / M_PI;

  return llh;
}


TEST( Datum, SetFromName) {
  // Verify that all of our standard names work
  Datum d;
  EXPECT_NO_THROW(d.set_well_known_datum("WGS84" ));
  EXPECT_NO_THROW(d.set_well_known_datum("WGS72" ));
  EXPECT_NO_THROW(d.set_well_known_datum("NAD83" ));
  EXPECT_NO_THROW(d.set_well_known_datum("NAD27" ));
  EXPECT_NO_THROW(d.set_well_known_datum("D_MOON"));
  EXPECT_NO_THROW(d.set_well_known_datum("D_MARS"));
  EXPECT_NO_THROW(d.set_well_known_datum("MOLA"  ));
  EXPECT_THROW(d.set_well_known_datum("MOOLA"), vw::InputErr);
}

TEST( Datum, SetFromProjString) {
  // This grid file does not exist but we don't check until a conversion is attempted.
  std::string proj_str = "+proj=longlat +datum=NAD83 +no_defs +nadgrids=fake_grid.ct2";
  Datum d;
  d.set_datum_from_proj_str(proj_str);
  EXPECT_EQ(proj_str, d.proj4_str());
  
  EXPECT_EQ("North_American_Datum_1983", d.name());
  EXPECT_EQ("GRS 1980", d.spheroid_name());
  EXPECT_EQ("Greenwich", d.meridian_name());
}

TEST( Datum, GeodeticConversion ) {
  Datum datum("WGS84");

  std::vector<double> values;
  values += -8000000.,-7000000.,-6500000.,-6100000.,-5788000.,-500000.,-100000.,-50000.,-10000.,-6000.,-1000.,-200.,-1.,0.,8000000.,7000000.,6500000.,6100000.,5788000.,500000.,100000.,50000.,10000.,6000.,1000.,200.,1;

  // The precision seems to be around a center when the vector is near
  // the center of the earth, which is the worse possible
  // situation. Things are better near the surface.
  for ( size_t ix = 0; ix < values.size(); ix++ ) {
    for ( size_t iy = 0; iy < values.size(); iy++ ) {
      for ( size_t iz = 0; iz < values.size(); iz++ ) {
        Vector3 test_xyz(values[ix],values[iy],values[iz]);
        EXPECT_VECTOR_NEAR( datum.geodetic_to_cartesian(vermille_2011_cart_to_geodetic(datum,test_xyz)),
                            test_xyz, std::max(1e-12 * norm_2(test_xyz), 1e-2) );
      }
    }
  }

  EXPECT_VECTOR_NEAR( datum.geodetic_to_cartesian(datum.cartesian_to_geodetic(Vector3(10000,10,-10))),
                      Vector3(10000,10,-10), 1e-6 );
  EXPECT_VECTOR_NEAR( datum.cartesian_to_geodetic(datum.geodetic_to_cartesian(Vector3(30,-10,173740))),
                      Vector3(30,-10,173740), 1e-6 );

  datum.set_well_known_datum("D_MOON"); // This is a spherical datum

  // The precision seems to be around a center when the vector is near
  // the center of the earth, which is the worse possible
  // situation. Things are better near the surface.
  for ( size_t ix = 0; ix < values.size(); ix++ ) {
    for ( size_t iy = 0; iy < values.size(); iy++ ) {
      for ( size_t iz = 0; iz < values.size(); iz++ ) {
        Vector3 test_xyz(values[ix],values[iy],values[iz]);
        EXPECT_VECTOR_NEAR( datum.geodetic_to_cartesian(vermille_2011_cart_to_geodetic(datum,test_xyz)),
                            test_xyz, std::max(1e-12 * norm_2(test_xyz), 1e-2) );
      }
    }
  }

  EXPECT_VECTOR_NEAR( datum.geodetic_to_cartesian(datum.cartesian_to_geodetic(Vector3(10000,10,-10))),
                      Vector3(10000,10,-10), 1e-6 );
  EXPECT_VECTOR_NEAR( datum.cartesian_to_geodetic(datum.geodetic_to_cartesian(Vector3(30,-10,173740))),
                      Vector3(30,-10,173740), 1e-6 );
}
