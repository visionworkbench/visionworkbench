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


#include <gtest/gtest_VW.h>

#include <vw/Camera/CAHVORModel.h>
#include <vw/Math/Vector.h>
#include <test/Helpers.h>

using namespace vw;
using namespace vw::camera;
using namespace vw::test;

TEST( CAHVORModel, PointDerivatives ) {
  CAHVORModel cahvor;

  // A random camera model to avoid symmetries.
  cahvor.C = Vector3(76,-34,20);
  cahvor.A = normalize(Vector3(1,2,3));
  cahvor.H = Vector3(1300,123,456);
  cahvor.V = Vector3(345,900,157);
  cahvor.O = normalize(Vector3(2,3,4));
  cahvor.R = Vector3(1,-1,1);
  Vector2 input_pixel(0,0);
  Vector3 grad_u, grad_v;
  Matrix3x3 hess_u, hess_v;
  Vector3 P(0,0,0);
  cahvor.get_point_derivatives( P, input_pixel[0], input_pixel[1],
                                grad_u, grad_v, hess_u, hess_v );
  Vector2 ref = cahvor.point_to_pixel( P );
  double delta = 0.001;
  // Calculating Derivatives by Hand
  Vector2 refgx = (cahvor.point_to_pixel( P+Vector3(delta,0,0) )-ref)/delta;
  Vector2 refgy = (cahvor.point_to_pixel( P+Vector3(0,delta,0) )-ref)/delta;
  Vector2 refgz = (cahvor.point_to_pixel( P+Vector3(0,0,delta) )-ref)/delta;
  Vector2 refhxx = ( ( cahvor.point_to_pixel( P+Vector3(2*delta,0,0) ) -
                       cahvor.point_to_pixel( P+Vector3(delta,0,0) ) ) -
                     ( cahvor.point_to_pixel( P+Vector3(delta,0,0) ) -
                       cahvor.point_to_pixel( P+Vector3(0,0,0) ) ) ) / (delta*delta);
  Vector2 refhxy = ( ( cahvor.point_to_pixel( P+Vector3(delta,delta,0) ) -
                       cahvor.point_to_pixel( P+Vector3(0,delta,0) ) ) -
                     ( cahvor.point_to_pixel( P+Vector3(delta,0,0) ) -
                       cahvor.point_to_pixel( P+Vector3(0,0,0) ) ) ) / (delta*delta);
  Vector2 refhxz = ( ( cahvor.point_to_pixel( P+Vector3(delta,0,delta) ) -
                       cahvor.point_to_pixel( P+Vector3(0,0,delta) ) ) -
                     ( cahvor.point_to_pixel( P+Vector3(delta,0,0) ) -
                       cahvor.point_to_pixel( P+Vector3(0,0,0) ) ) ) / (delta*delta);
  Vector2 refhyx = ( ( cahvor.point_to_pixel( P+Vector3(delta,delta,0) ) -
                       cahvor.point_to_pixel( P+Vector3(delta,0,0) ) ) -
                     ( cahvor.point_to_pixel( P+Vector3(0,delta,0) ) -
                       cahvor.point_to_pixel( P+Vector3(0,0,0) ) ) ) / (delta*delta);
  Vector2 refhyy = ( ( cahvor.point_to_pixel( P+Vector3(0,2*delta,0) ) -
                       cahvor.point_to_pixel( P+Vector3(0,delta,0) ) ) -
                     ( cahvor.point_to_pixel( P+Vector3(0,delta,0) ) -
                       cahvor.point_to_pixel( P+Vector3(0,0,0) ) ) ) / (delta*delta);
  Vector2 refhyz = ( ( cahvor.point_to_pixel( P+Vector3(0,delta,delta) ) -
                       cahvor.point_to_pixel( P+Vector3(0,0,delta) ) ) -
                     ( cahvor.point_to_pixel( P+Vector3(0,delta,0) ) -
                       cahvor.point_to_pixel( P+Vector3(0,0,0) ) ) ) / (delta*delta);
  Vector2 refhzx = ( ( cahvor.point_to_pixel( P+Vector3(delta,0,delta) ) -
                       cahvor.point_to_pixel( P+Vector3(delta,0,0) ) ) -
                     ( cahvor.point_to_pixel( P+Vector3(0,0,delta) ) -
                       cahvor.point_to_pixel( P+Vector3(0,0,0) ) ) ) / (delta*delta);
  Vector2 refhzy = ( ( cahvor.point_to_pixel( P+Vector3(0,delta,delta) ) -
                       cahvor.point_to_pixel( P+Vector3(0,delta,0) ) ) -
                     ( cahvor.point_to_pixel( P+Vector3(0,0,delta) ) -
                       cahvor.point_to_pixel( P+Vector3(0,0,0) ) ) ) / (delta*delta);
  Vector2 refhzz = ( ( cahvor.point_to_pixel( P+Vector3(0,0,2*delta) ) -
                       cahvor.point_to_pixel( P+Vector3(0,0,delta) ) ) -
                     ( cahvor.point_to_pixel( P+Vector3(0,0,delta) ) -
                       cahvor.point_to_pixel( P+Vector3(0,0,0) ) ) ) / (delta*delta);

  // Compare Results
  EXPECT_VECTOR_NEAR( input_pixel, ref, 1e-5 );

  EXPECT_NEAR( grad_u[0], refgx[0], 0.01 );
  EXPECT_NEAR( grad_u[1], refgy[0], 0.01 );
  EXPECT_NEAR( grad_u[2], refgz[0], 0.01 );
  EXPECT_NEAR( grad_v[0], refgx[1], 0.01 );
  EXPECT_NEAR( grad_v[1], refgy[1], 0.01 );
  EXPECT_NEAR( grad_v[2], refgz[1], 0.01 );

  EXPECT_NEAR( hess_u(0,0), refhxx[0], 0.01 );
  EXPECT_NEAR( hess_u(1,0), refhyx[0], 0.01 );
  EXPECT_NEAR( hess_u(2,0), refhzx[0], 0.01 );
  EXPECT_NEAR( hess_u(0,1), refhxy[0], 0.01 );
  EXPECT_NEAR( hess_u(1,1), refhyy[0], 0.01 );
  EXPECT_NEAR( hess_u(2,1), refhzy[0], 0.01 );
  EXPECT_NEAR( hess_u(0,2), refhxz[0], 0.01 );
  EXPECT_NEAR( hess_u(1,2), refhyz[0], 0.01 );
  EXPECT_NEAR( hess_u(2,2), refhzz[0], 0.01 );
  EXPECT_NEAR( hess_v(0,0), refhxx[1], 0.01 );
  EXPECT_NEAR( hess_v(1,0), refhyx[1], 0.01 );
  EXPECT_NEAR( hess_v(2,0), refhzx[1], 0.01 );
  EXPECT_NEAR( hess_v(0,1), refhxy[1], 0.01 );
  EXPECT_NEAR( hess_v(1,1), refhyy[1], 0.01 );
  EXPECT_NEAR( hess_v(2,1), refhzy[1], 0.01 );
  EXPECT_NEAR( hess_v(0,2), refhxz[1], 0.01 );
  EXPECT_NEAR( hess_v(1,2), refhyz[1], 0.01 );
  EXPECT_NEAR( hess_v(2,2), refhzz[1], 0.01 );
}

TEST( CAHVORModel, StringRead ) {
  CAHVORModel c(TEST_SRCDIR"/cahvor.txt");

  Vector3 C = Vector3(76,-34,20);
  Vector3 A = normalize(Vector3(1,2,3));
  Vector3 H = Vector3(1300,123,456);
  Vector3 V = Vector3(345,900,157);
  Vector3 O = normalize(Vector3(2,3,4));
  Vector3 R = Vector3(1,-1,1);

  EXPECT_VECTOR_NEAR( C, c.C, 1e-5 );

  EXPECT_VECTOR_NEAR( A, c.A, 1e-5 );
  EXPECT_VECTOR_NEAR( H, c.H, 1e-5 );
  EXPECT_VECTOR_NEAR( V, c.V, 1e-5 );
  EXPECT_VECTOR_NEAR( O, c.O, 1e-5 );
  EXPECT_VECTOR_NEAR( R, c.R, 1e-5 );
}

TEST( CAHVORModel, StringWriteRead ) {
  CAHVORModel c;
  UnlinkName file("StringWriteRead.txt");

  c.C = Vector3(76,-34,20);
  c.A = normalize(Vector3(1,2,3));
  c.H = Vector3(1300,123,456);
  c.V = Vector3(345,900,157);
  c.O = normalize(Vector3(2,3,4));
  c.R = Vector3(1,-1,1);
  c.write(file);

  CAHVORModel c2(file);

  EXPECT_VECTOR_NEAR( c2.C, c.C, 1e-5 );
  EXPECT_VECTOR_NEAR( c2.C, c.camera_center(), 1e-5 );
  EXPECT_VECTOR_NEAR( c2.A, c.A, 1e-5 );
  EXPECT_VECTOR_NEAR( c2.H, c.H, 1e-5 );
  EXPECT_VECTOR_NEAR( c2.V, c.V, 1e-5 );
  EXPECT_VECTOR_NEAR( c2.O, c.O, 1e-5 );
  EXPECT_VECTOR_NEAR( c2.R, c.R, 1e-5 );

  EXPECT_STREQ( "CAHVOR", c2.type().c_str() );
}

TEST( CAHVORModel, ForwardReverse) {
  CAHVORModel cahvor(Vector3(0.491222,-0.0717236,-1.24143),
                     Vector3(0.921657,-0.230518,0.312107),
                     Vector3(757.076,1071.6,160.227),
                     Vector3(91.7479,-27.7504,1319.48),
                     Vector3(0.920759,-0.206185,0.331197),
                     Vector3(0.00096,-0.002183,0.018547));

  for ( uint32 i = 100; i < 901; i += 100 ) {
    for ( uint32 j = 100; j < 901; j+= 100 ) {
      Vector2 pixel( i,j );
      Vector3 unit = cahvor.pixel_to_vector( pixel );
      Vector3 point = cahvor.C + 30*unit;
      Vector2 result = cahvor.point_to_pixel( point );
      EXPECT_VECTOR_NEAR( result, pixel, 1e-2 );
    }
  }
}
