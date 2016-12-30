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

#include <test/Helpers.h>
#include <gtest/gtest_VW.h>
#include <vw/Math/Functions.h>

// This is tested here to consolidate otherwise tiny files
#include <vw/Math/BresenhamLine.h>

using namespace vw;

static const double DELTA = 1e-15;

TEST(Functions, ERF) {
  EXPECT_NEAR( -0.9999779095030014,    vw::math::impl::erf(-3.0)   , DELTA);
  EXPECT_NEAR( -0.8427007929497149,    vw::math::impl::erf(-1e0)   , DELTA);
  EXPECT_NEAR( -0.1124629160182849,    vw::math::impl::erf(-1e-1)  , DELTA);
  EXPECT_NEAR( -0.01128341555584962,   vw::math::impl::erf(-1e-2)  , DELTA);
  EXPECT_NEAR( -0.0001128379163334249, vw::math::impl::erf(-1e-4)  , DELTA);
  EXPECT_NEAR( 0,                      vw::math::impl::erf(0)      , DELTA);
  EXPECT_NEAR( 0.0001128379163334249,  vw::math::impl::erf(1e-4)   , DELTA);
  EXPECT_NEAR( 0.01128341555584962,    vw::math::impl::erf(1e-2)   , DELTA);
  EXPECT_NEAR( 0.1124629160182849,     vw::math::impl::erf(1e-1)   , DELTA);
  EXPECT_NEAR( 0.8427007929497149,     vw::math::impl::erf(1e0)    , DELTA);
  EXPECT_NEAR( 0.9999779095030014,     vw::math::impl::erf(3.0)    , DELTA);
}

TEST(Functions, ERFC) {
  EXPECT_NEAR( 1.999977909503001,      vw::math::impl::erfc(-3.0)  , DELTA);
  EXPECT_NEAR( 1.842700792949715,      vw::math::impl::erfc(-1e0)  , DELTA);
  EXPECT_NEAR( 1.112462916018285,      vw::math::impl::erfc(-1e-1) , DELTA);
  EXPECT_NEAR( 1.011283415555850,      vw::math::impl::erfc(-1e-2) , DELTA);
  EXPECT_NEAR( 1.000112837916333,      vw::math::impl::erfc(-1e-4) , DELTA);
  EXPECT_NEAR( 1,                      vw::math::impl::erfc(0)     , DELTA);
  EXPECT_NEAR( 0.9998871620836666,     vw::math::impl::erfc(1e-4)  , DELTA);
  EXPECT_NEAR( 0.9887165844441504,     vw::math::impl::erfc(1e-2)  , DELTA);
  EXPECT_NEAR( 0.8875370839817151,     vw::math::impl::erfc(1e-1)  , DELTA);
  EXPECT_NEAR( 0.1572992070502851,     vw::math::impl::erfc(1e0)   , DELTA);
  EXPECT_NEAR( 0.00002209049699858544, vw::math::impl::erfc(3.0)   , DELTA);
}

TEST(BresenhamLine, BresenhamLine) {

  vw::math::BresenhamLine lineA(0,0,  5,10);
  vw::math::BresenhamLine lineB(0,0, 10, 5);

  // Check steep line
  EXPECT_VECTOR_EQ(Vector2i(0,0), *lineA); lineA++;
  EXPECT_VECTOR_EQ(Vector2i(0,1), *lineA); lineA++;
  EXPECT_VECTOR_EQ(Vector2i(1,2), *lineA); lineA++;
  EXPECT_VECTOR_EQ(Vector2i(1,3), *lineA); lineA++;

  // Check shallow line
  EXPECT_VECTOR_EQ(Vector2i(0,0), *lineB); lineB++;
  EXPECT_VECTOR_EQ(Vector2i(1,0), *lineB); lineB++;
  EXPECT_VECTOR_EQ(Vector2i(2,1), *lineB); lineB++;
  EXPECT_VECTOR_EQ(Vector2i(3,1), *lineB); lineB++;
  
  // Check that the line ends where we expect
  lineB++;
  lineB++;
  lineB++;
  lineB++;
  lineB++;
  EXPECT_VECTOR_EQ(Vector2i(9,4), *lineB);
  lineB++;
  EXPECT_FALSE(lineB.is_good());
  
}

