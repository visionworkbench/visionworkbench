// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
//
// Copyright 2006 Carnegie Mellon University. All rights reserved.
//
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
//
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
//
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
