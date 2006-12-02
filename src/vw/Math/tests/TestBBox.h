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

// TestVector.h
#include <cxxtest/TestSuite.h>
#include <vw/Math/BBox.h>

using namespace vw;

class TestVector : public CxxTest::TestSuite
{
public:

  void test_bbox()
  {
    // Default constructor (floating-point type)
    // Hard to know how to test this one....
    BBox<float,2> b1;

    // Default constructor (integer type)
    BBox<uint8,1> b2;
    TS_ASSERT_EQUALS( b2.min()[0], 255 );
    TS_ASSERT_EQUALS( b2.max()[0], 0 );
    TS_ASSERT( b2.empty() );

    // Two-vector constructor
    BBox2 b3( Vector2(1,2), Vector2(3,4) );
    TS_ASSERT_EQUALS( b3.min()[0], 1 );
    TS_ASSERT_EQUALS( b3.min()[1], 2 );
    TS_ASSERT_EQUALS( b3.max()[0], 3 );
    TS_ASSERT_EQUALS( b3.max()[1], 4 );

    // Four-number constructor
    BBox2 b4(1,2,1,4);
    TS_ASSERT_EQUALS( b4.min()[0], 1 );
    TS_ASSERT_EQUALS( b4.min()[1], 2 );
    TS_ASSERT_EQUALS( b4.max()[0], 2 );
    TS_ASSERT_EQUALS( b4.max()[1], 6 );

    // Assorted accessors
    TS_ASSERT_EQUALS( b4.size()[0], 1 );
    TS_ASSERT_EQUALS( b4.size()[1], 4 );
    TS_ASSERT_EQUALS( b4.center()[0], 1.5 );
    TS_ASSERT_EQUALS( b4.center()[1], 4 );
    TS_ASSERT_EQUALS( b4.min()[0], 1 );
    TS_ASSERT_EQUALS( b4.min()[1], 2 );
    TS_ASSERT_EQUALS( b4.max()[0], 2 );
    TS_ASSERT_EQUALS( b4.max()[1], 6 );
    TS_ASSERT_EQUALS( b4.width(), 1 );
    TS_ASSERT_EQUALS( b4.height(), 4 );
  }

  void test_bbox_ops()
  {
    BBox2 b4(1,2,1,4);

    // Grow to include a point
    b4.grow( Vector2(0,7) );
    TS_ASSERT_EQUALS( b4.min()[0], 0 );
    TS_ASSERT_EQUALS( b4.min()[1], 2 );
    TS_ASSERT_EQUALS( b4.max()[0], 2 );
    TS_ASSERT_EQUALS( b4.max()[1], 7 );

    // Grow to include a bbox
    b4.grow( BBox2(0,3,3,3) );
    TS_ASSERT_EQUALS( b4.min()[0], 0 );
    TS_ASSERT_EQUALS( b4.min()[1], 2 );
    TS_ASSERT_EQUALS( b4.max()[0], 3 );
    TS_ASSERT_EQUALS( b4.max()[1], 7 );
    
    // Crop to bbox
    b4.crop( BBox2(2,0,1,4) );
    TS_ASSERT_EQUALS( b4.min()[0], 2 );
    TS_ASSERT_EQUALS( b4.min()[1], 2 );
    TS_ASSERT_EQUALS( b4.max()[0], 3 );
    TS_ASSERT_EQUALS( b4.max()[1], 4 );
    
    // Expand
    b4.expand(2);
    TS_ASSERT_EQUALS( b4.min()[0], 0 );
    TS_ASSERT_EQUALS( b4.min()[1], 0 );
    TS_ASSERT_EQUALS( b4.max()[0], 5 );
    TS_ASSERT_EQUALS( b4.max()[1], 6 );
    
    // Contract
    b4.contract(1);
    TS_ASSERT_EQUALS( b4.min()[0], 1 );
    TS_ASSERT_EQUALS( b4.min()[1], 1 );
    TS_ASSERT_EQUALS( b4.max()[0], 4 );
    TS_ASSERT_EQUALS( b4.max()[1], 5 );

    // Contains points
    TS_ASSERT_EQUALS( b4.contains( Vector2(2,3) ), true );
    TS_ASSERT_EQUALS( b4.contains( Vector2(0,3) ), false );
    TS_ASSERT_EQUALS( b4.contains( Vector2(2,6) ), false );
    TS_ASSERT_EQUALS( b4.contains( b4.min() ), true );
    TS_ASSERT_EQUALS( b4.contains( b4.max() ), false );

    // Contains bboxes
    TS_ASSERT_EQUALS( b4.contains( BBox2(2,2,2,2) ), true );
    TS_ASSERT_EQUALS( b4.contains( BBox2(0,2,4,2) ), false );
    TS_ASSERT_EQUALS( b4.contains( BBox2(2,0,2,4) ), false );
    TS_ASSERT_EQUALS( b4.contains( BBox2(2,2,4,2) ), false );
    TS_ASSERT_EQUALS( b4.contains( BBox2(2,2,4,4) ), false );
    TS_ASSERT_EQUALS( b4.contains( b4 ), true );
    
    // Intersects
    TS_ASSERT_EQUALS( b4.intersects( BBox2(2,2,2,2) ), true );
    TS_ASSERT_EQUALS( b4.intersects( BBox2(0,2,4,2) ), true );
    TS_ASSERT_EQUALS( b4.intersects( BBox2(2,2,4,4) ), true );
    TS_ASSERT_EQUALS( b4.intersects( BBox2(0,0,1,1) ), false );
    TS_ASSERT_EQUALS( b4.intersects( BBox2(0,0,2,2) ), true );
    TS_ASSERT_EQUALS( b4.intersects( BBox2(2,5,2,2) ), false );

    // Empty
    TS_ASSERT_EQUALS( b4.empty(), false );
    b4.min()[0] = b4.max()[0];
    TS_ASSERT_EQUALS( b4.empty(), true );
  }
  
  void test_bbox_math()
  {
    BBox2 b1(1,2,2,2), b2;

    // Scale-assignment
    b1 *= 2;
    TS_ASSERT_EQUALS( b1.min()[0], 2 );
    TS_ASSERT_EQUALS( b1.min()[1], 4 );
    TS_ASSERT_EQUALS( b1.max()[0], 6 );
    TS_ASSERT_EQUALS( b1.max()[1], 8 );

    // Sum-assignment
    b1 += Vector2(2,3);
    TS_ASSERT_EQUALS( b1.min()[0], 4 );
    TS_ASSERT_EQUALS( b1.min()[1], 7 );
    TS_ASSERT_EQUALS( b1.max()[0], 8 );
    TS_ASSERT_EQUALS( b1.max()[1], 11 );

    // Difference-assignment
    b1 -= Vector2(4,4);
    TS_ASSERT_EQUALS( b1.min()[0], 0 );
    TS_ASSERT_EQUALS( b1.min()[1], 3 );
    TS_ASSERT_EQUALS( b1.max()[0], 4 );
    TS_ASSERT_EQUALS( b1.max()[1], 7 );
    
    // Left scale
    b2 = 2*b1;
    TS_ASSERT_EQUALS( b2.min()[0], 0 );
    TS_ASSERT_EQUALS( b2.min()[1], 6 );
    TS_ASSERT_EQUALS( b2.max()[0], 8 );
    TS_ASSERT_EQUALS( b2.max()[1], 14 );

    // Right scale
    b2 = b1*3;
    TS_ASSERT_EQUALS( b2.min()[0], 0 );
    TS_ASSERT_EQUALS( b2.min()[1], 9 );
    TS_ASSERT_EQUALS( b2.max()[0], 12 );
    TS_ASSERT_EQUALS( b2.max()[1], 21 );

    // Right quotient
    b2 = b1/2;
    TS_ASSERT_EQUALS( b2.min()[0], 0 );
    TS_ASSERT_EQUALS( b2.min()[1], 1.5 );
    TS_ASSERT_EQUALS( b2.max()[0], 2 );
    TS_ASSERT_EQUALS( b2.max()[1], 3.5 );

    // Right vector sum
    b2 = b1 + Vector2(1,1);
    TS_ASSERT_EQUALS( b2.min()[0], 1 );
    TS_ASSERT_EQUALS( b2.min()[1], 4 );
    TS_ASSERT_EQUALS( b2.max()[0], 5 );
    TS_ASSERT_EQUALS( b2.max()[1], 8 );

    // Left vector sum
    b2 = Vector2(1,1) + b1;
    TS_ASSERT_EQUALS( b2.min()[0], 1 );
    TS_ASSERT_EQUALS( b2.min()[1], 4 );
    TS_ASSERT_EQUALS( b2.max()[0], 5 );
    TS_ASSERT_EQUALS( b2.max()[1], 8 );

    // Vector sum
    b2 = b1 - Vector2(1,1);
    TS_ASSERT_EQUALS( b2.min()[0], -1 );
    TS_ASSERT_EQUALS( b2.min()[1], 2 );
    TS_ASSERT_EQUALS( b2.max()[0], 3 );
    TS_ASSERT_EQUALS( b2.max()[1], 6 );
  }

}; // class TestVector
