// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>
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

