// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestBBall.h

#include <iostream>

#include <cxxtest/TestSuite.h>
#include <vw/Math/Vector.h>
#include <vw/Geometry/Sphere.h>

using namespace vw;

class TestSphere : public CxxTest::TestSuite
{
public:

  void test_bball()
  {
    SphereN b0;
    TS_ASSERT( b0.empty() );
    b0.grow(Vector3(0, 0, 0));
    TS_ASSERT( !b0.empty() );
    b0.grow(Vector3(1, 0, 0));
    TS_ASSERT( !b0.empty() );

    // unit ball centered at (0,0,0)
    Sphere3 b1;
    b1.grow(Vector3(0, 0, 0));
    b1.grow(Vector3(0, 0, 1));
    b1.grow(Vector3(0, 0, -1));
    b1.grow(Vector3(0, 1, 0));
    b1.grow(Vector3(0, -1, 0));
    b1.grow(Vector3(1, 0, 0));
    b1.grow(Vector3(-1, 0, 0));
    TS_ASSERT_EQUALS(b1.center(), Vector3(0, 0, 0));
    TS_ASSERT( b1.contains(Vector3(0.5, 0.5, 0.5)) );
    TS_ASSERT( !b1.contains(Vector3(0.9, 0.9, 1.5)) );
    TS_ASSERT( b1.intersects(b1) );

    Sphere3 b2( Vector3( 0, 0, 0 ), 0 );
    TS_ASSERT( !b2.empty() );
  }

}; // class TestSphere
