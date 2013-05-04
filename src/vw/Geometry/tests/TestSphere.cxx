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
#include <vw/Geometry/Sphere.h>
#include <vw/Math/Vector.h>

using namespace vw;

TEST(TestSphere, Sphere) {
  SphereN b0;
  ASSERT_TRUE( b0.empty() );
  b0.grow(Vector3(0, 0, 0));
  ASSERT_TRUE( !b0.empty() );
  b0.grow(Vector3(1, 0, 0));
  ASSERT_TRUE( !b0.empty() );

  // unit ball centered at (0,0,0)
  Sphere3 b1;
  b1.grow(Vector3(0, 0, 0));
  b1.grow(Vector3(0, 0, 1));
  b1.grow(Vector3(0, 0, -1));
  b1.grow(Vector3(0, 1, 0));
  b1.grow(Vector3(0, -1, 0));
  b1.grow(Vector3(1, 0, 0));
  b1.grow(Vector3(-1, 0, 0));
  EXPECT_EQ( Vector3(0, 0, 0), b1.center() );
  EXPECT_TRUE( b1.contains(Vector3(0.5, 0.5, 0.5)) );
  EXPECT_TRUE( !b1.contains(Vector3(0.9, 0.9, 1.5)) );
  EXPECT_TRUE( b1.intersects(b1) );

  Sphere3 b2( Vector3( 0, 0, 0 ), 0 );
  EXPECT_TRUE( !b2.empty() );
}

