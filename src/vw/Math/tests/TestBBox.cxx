// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest_VW.h>
#include <vw/Math/Vector.h>
#include <vw/Math/BBox.h>

using namespace vw;

TEST(BBox, Static) {
  // Default constructor (floating-point type)
  // Hard to know how to test this one....
  BBox<float,2> b1;

  // Default constructor (integer type)
  BBox<uint8,1> b2;
  EXPECT_EQ( 255, b2.min()[0] );
  EXPECT_EQ( 0, b2.max()[0] );
  EXPECT_TRUE( b2.empty() );

  // Two-vector constructor
  BBox2 b3( Vector2(1,2), Vector2(3,4) );
  EXPECT_DOUBLE_EQ( 1, b3.min()[0] );
  EXPECT_DOUBLE_EQ( 2, b3.min()[1] );
  EXPECT_DOUBLE_EQ( 3, b3.max()[0] );
  EXPECT_DOUBLE_EQ( 4, b3.max()[1] );

  // Four-number constructor
  BBox2 b4(1,2,1,4);
  EXPECT_DOUBLE_EQ( 1, b4.min()[0] );
  EXPECT_DOUBLE_EQ( 2, b4.min()[1] );
  EXPECT_DOUBLE_EQ( 2, b4.max()[0] );
  EXPECT_DOUBLE_EQ( 6, b4.max()[1] );

  // Six-number constructor
  BBox3 b5(1,2,3,1,4,6);
  EXPECT_DOUBLE_EQ( 1, b5.min()[0] );
  EXPECT_DOUBLE_EQ( 2, b5.min()[1] );
  EXPECT_DOUBLE_EQ( 3, b5.min()[2] );
  EXPECT_DOUBLE_EQ( 2, b5.max()[0] );
  EXPECT_DOUBLE_EQ( 6, b5.max()[1] );
  EXPECT_DOUBLE_EQ( 9, b5.max()[2] );

  // Copy constructor
  BBox3 b6(b5);
  EXPECT_DOUBLE_EQ( 1, b6.min()[0] );
  EXPECT_DOUBLE_EQ( 2, b6.min()[1] );
  EXPECT_DOUBLE_EQ( 3, b6.min()[2] );
  EXPECT_DOUBLE_EQ( 2, b6.max()[0] );
  EXPECT_DOUBLE_EQ( 6, b6.max()[1] );
  EXPECT_DOUBLE_EQ( 9, b6.max()[2] );

  // Assorted accessors
  EXPECT_EQ( 1, b6.size()[0] );
  EXPECT_EQ( 4, b6.size()[1] );
  EXPECT_EQ( 6, b6.size()[2] );
  EXPECT_DOUBLE_EQ( 1.5, b6.center()[0] );
  EXPECT_DOUBLE_EQ( 4, b6.center()[1] );
  EXPECT_DOUBLE_EQ( 6, b6.center()[2] );
  EXPECT_DOUBLE_EQ( 1, b6.min()[0] );
  EXPECT_DOUBLE_EQ( 2, b6.min()[1] );
  EXPECT_DOUBLE_EQ( 3, b6.min()[2] );
  EXPECT_DOUBLE_EQ( 2, b6.max()[0] );
  EXPECT_DOUBLE_EQ( 6, b6.max()[1] );
  EXPECT_DOUBLE_EQ( 9, b6.max()[2] );
  EXPECT_EQ( 1, b6.width() );
  EXPECT_EQ( 4, b6.height() );
  EXPECT_EQ( 6, b6.depth() );
}

TEST(BBox, Dynamic) {
  // Default constructor (floating-point type)
  // Hard to know how to test this one....
  BBox<float> b1;

  // Default constructor (integer type)
  BBox<uint8> b2;
  EXPECT_TRUE( b2.empty() );

  // Two-vector constructor
  BBoxN b3( Vector2(1,2), Vector2(3,4) );
  EXPECT_DOUBLE_EQ( 1, b3.min()[0] );
  EXPECT_DOUBLE_EQ( 2, b3.min()[1] );
  EXPECT_DOUBLE_EQ( 3, b3.max()[0] );
  EXPECT_DOUBLE_EQ( 4, b3.max()[1] );

  // Four-number constructor
  BBoxN b4(1,2,1,4);
  EXPECT_DOUBLE_EQ( 1, b4.min()[0] );
  EXPECT_DOUBLE_EQ( 2, b4.min()[1] );
  EXPECT_DOUBLE_EQ( 2, b4.max()[0] );
  EXPECT_DOUBLE_EQ( 6, b4.max()[1] );

  // Six-number constructor
  BBoxN b5(1,2,3,1,4,6);
  EXPECT_DOUBLE_EQ( 1, b5.min()[0] );
  EXPECT_DOUBLE_EQ( 2, b5.min()[1] );
  EXPECT_DOUBLE_EQ( 3, b5.min()[2] );
  EXPECT_DOUBLE_EQ( 2, b5.max()[0] );
  EXPECT_DOUBLE_EQ( 6, b5.max()[1] );
  EXPECT_DOUBLE_EQ( 9, b5.max()[2] );

  // Copy constructor
  BBoxN b6(b5);
  EXPECT_DOUBLE_EQ( 1, b6.min()[0] );
  EXPECT_DOUBLE_EQ( 2, b6.min()[1] );
  EXPECT_DOUBLE_EQ( 3, b6.min()[2] );
  EXPECT_DOUBLE_EQ( 2, b6.max()[0] );
  EXPECT_DOUBLE_EQ( 6, b6.max()[1] );
  EXPECT_DOUBLE_EQ( 9, b6.max()[2] );

  // Assorted accessors
  EXPECT_EQ( 1, b6.size()[0] );
  EXPECT_EQ( 4, b6.size()[1] );
  EXPECT_EQ( 6, b6.size()[2] );
  EXPECT_DOUBLE_EQ( 1.5, b6.center()[0] );
  EXPECT_DOUBLE_EQ( 4, b6.center()[1] );
  EXPECT_DOUBLE_EQ( 6, b6.center()[2] );
  EXPECT_DOUBLE_EQ( 1, b6.min()[0] );
  EXPECT_DOUBLE_EQ( 2, b6.min()[1] );
  EXPECT_DOUBLE_EQ( 3, b6.min()[2] );
  EXPECT_DOUBLE_EQ( 2, b6.max()[0] );
  EXPECT_DOUBLE_EQ( 6, b6.max()[1] );
  EXPECT_DOUBLE_EQ( 9, b6.max()[2] );
  EXPECT_EQ( 1, b6.width() );
  EXPECT_EQ( 4, b6.height() );
  EXPECT_EQ( 6, b6.depth() );
}

TEST(BBox, Ops) {
  BBox2 b4(1,2,1,4);
  BBoxN b5(1,2,1,4);

  // Grow to include a point
  b4.grow( Vector2(0,7) );
  EXPECT_DOUBLE_EQ( 0, b4.min()[0] );
  EXPECT_DOUBLE_EQ( 2, b4.min()[1] );
  EXPECT_DOUBLE_EQ( 2, b4.max()[0] );
  EXPECT_DOUBLE_EQ( 7, b4.max()[1] );

  b5.grow( Vector2(0,7) );
  EXPECT_DOUBLE_EQ( 0, b5.min()[0] );
  EXPECT_DOUBLE_EQ( 2, b5.min()[1] );
  EXPECT_DOUBLE_EQ( 2, b5.max()[0] );
  EXPECT_DOUBLE_EQ( 7, b5.max()[1] );

  // Grow to include a bbox
  b4.grow( BBox2(0,3,3,3) );
  EXPECT_DOUBLE_EQ( 0, b4.min()[0] );
  EXPECT_DOUBLE_EQ( 2, b4.min()[1] );
  EXPECT_DOUBLE_EQ( 3, b4.max()[0] );
  EXPECT_DOUBLE_EQ( 7, b4.max()[1] );

  b5.grow( BBox2(0,3,3,3) );
  EXPECT_DOUBLE_EQ( 0, b5.min()[0] );
  EXPECT_DOUBLE_EQ( 2, b5.min()[1] );
  EXPECT_DOUBLE_EQ( 3, b5.max()[0] );
  EXPECT_DOUBLE_EQ( 7, b5.max()[1] );

  // Crop to bbox
  b4.crop( BBox2(2,0,1,4) );
  EXPECT_DOUBLE_EQ( 2, b4.min()[0] );
  EXPECT_DOUBLE_EQ( 2, b4.min()[1] );
  EXPECT_DOUBLE_EQ( 3, b4.max()[0] );
  EXPECT_DOUBLE_EQ( 4, b4.max()[1] );

  b5.crop( BBox2(2,0,1,4) );
  EXPECT_DOUBLE_EQ( 2, b5.min()[0] );
  EXPECT_DOUBLE_EQ( 2, b5.min()[1] );
  EXPECT_DOUBLE_EQ( 3, b5.max()[0] );
  EXPECT_DOUBLE_EQ( 4, b5.max()[1] );

  // Expand
  b4.expand(2);
  EXPECT_DOUBLE_EQ( 0, b4.min()[0] );
  EXPECT_DOUBLE_EQ( 0, b4.min()[1] );
  EXPECT_DOUBLE_EQ( 5, b4.max()[0] );
  EXPECT_DOUBLE_EQ( 6, b4.max()[1] );

  b5.expand(2);
  EXPECT_DOUBLE_EQ( 0, b5.min()[0] );
  EXPECT_DOUBLE_EQ( 0, b5.min()[1] );
  EXPECT_DOUBLE_EQ( 5, b5.max()[0] );
  EXPECT_DOUBLE_EQ( 6, b5.max()[1] );

  // Contract
  b4.contract(1);
  EXPECT_DOUBLE_EQ( 1, b4.min()[0] );
  EXPECT_DOUBLE_EQ( 1, b4.min()[1] );
  EXPECT_DOUBLE_EQ( 4, b4.max()[0] );
  EXPECT_DOUBLE_EQ( 5, b4.max()[1] );

  b5.contract(1);
  EXPECT_DOUBLE_EQ( 1, b5.min()[0] );
  EXPECT_DOUBLE_EQ( 1, b5.min()[1] );
  EXPECT_DOUBLE_EQ( 4, b5.max()[0] );
  EXPECT_DOUBLE_EQ( 5, b5.max()[1] );

  // Contains points
  EXPECT_TRUE  ( b4.contains ( Vector2(2,3) ));
  EXPECT_FALSE ( b4.contains ( Vector2(0,3) ));
  EXPECT_FALSE ( b4.contains ( Vector2(2,6) ));
  EXPECT_TRUE  ( b4.contains ( b4.min() ));
  EXPECT_FALSE ( b4.contains ( b4.max() ));

  EXPECT_TRUE  ( b5.contains ( Vector2(2,3) ));
  EXPECT_FALSE ( b5.contains ( Vector2(0,3) ));
  EXPECT_FALSE ( b5.contains ( Vector2(2,6) ));
  EXPECT_TRUE  ( b5.contains ( b4.min() ));
  EXPECT_FALSE ( b5.contains ( b4.max() ));

  // Contains bboxes
  EXPECT_TRUE  ( b4.contains ( BBox2(2,2,2,2) ));
  EXPECT_TRUE  ( b4.contains ( BBoxN(2,2,2,2) ));
  EXPECT_FALSE ( b4.contains ( BBox2(0,2,4,2) ));
  EXPECT_FALSE ( b4.contains ( BBoxN(2,0,2,4) ));
  EXPECT_FALSE ( b4.contains ( BBox2(2,2,4,2) ));
  EXPECT_FALSE ( b4.contains ( BBoxN(2,2,4,4) ));
  EXPECT_TRUE  ( b4.contains ( b4 ));

  EXPECT_TRUE  ( b5.contains ( BBox2(2,2,2,2) ));
  EXPECT_TRUE  ( b5.contains ( BBoxN(2,2,2,2) ));
  EXPECT_FALSE ( b5.contains ( BBox2(0,2,4,2) ));
  EXPECT_FALSE ( b5.contains ( BBoxN(2,0,2,4) ));
  EXPECT_FALSE ( b5.contains ( BBox2(2,2,4,2) ));
  EXPECT_FALSE ( b5.contains ( BBoxN(2,2,4,4) ));
  EXPECT_TRUE  ( b5.contains ( b4 ));

  // Intersects
  EXPECT_TRUE  ( b4.intersects ( BBox2(2,2,2,2) ));
  EXPECT_TRUE  ( b4.intersects ( BBoxN(0,2,4,2) ));
  EXPECT_TRUE  ( b4.intersects ( BBox2(2,2,4,4) ));
  EXPECT_FALSE ( b4.intersects ( BBoxN(0,0,1,1) ));
  EXPECT_TRUE  ( b4.intersects ( BBox2(0,0,2,2) ));
  EXPECT_FALSE ( b4.intersects ( BBox2(2,5,2,2) ));

  EXPECT_TRUE  ( b5.intersects ( BBox2(2,2,2,2) ));
  EXPECT_TRUE  ( b5.intersects ( BBoxN(0,2,4,2) ));
  EXPECT_TRUE  ( b5.intersects ( BBox2(2,2,4,4) ));
  EXPECT_FALSE ( b5.intersects ( BBoxN(0,0,1,1) ));
  EXPECT_TRUE  ( b5.intersects ( BBox2(0,0,2,2) ));
  EXPECT_FALSE ( b5.intersects ( BBox2(2,5,2,2) ));

  // Empty
  EXPECT_FALSE( b4.empty());
  b4.min()[0] = b4.max()[0];
  EXPECT_TRUE( b4.empty());

  EXPECT_FALSE( b5.empty());
  b5.min()[0] = b5.max()[0];
  EXPECT_TRUE( b5.empty());
}

TEST(BBox, Math) {
  BBox2 b1(1,2,2,2), b2;
  BBoxN b3(1,2,2,2);

  // Scale-assignment
  b1 *= 2;
  EXPECT_DOUBLE_EQ( 2, b1.min()[0] );
  EXPECT_DOUBLE_EQ( 4, b1.min()[1] );
  EXPECT_DOUBLE_EQ( 6, b1.max()[0] );
  EXPECT_DOUBLE_EQ( 8, b1.max()[1] );

  b3 *= 2;
  EXPECT_DOUBLE_EQ( 2, b3.min()[0] );
  EXPECT_DOUBLE_EQ( 4, b3.min()[1] );
  EXPECT_DOUBLE_EQ( 6, b3.max()[0] );
  EXPECT_DOUBLE_EQ( 8, b3.max()[1] );

  // Sum-assignment
  b1 += Vector2(2,3);
  EXPECT_DOUBLE_EQ( 4, b1.min()[0] );
  EXPECT_DOUBLE_EQ( 7, b1.min()[1] );
  EXPECT_DOUBLE_EQ( 8, b1.max()[0] );
  EXPECT_DOUBLE_EQ( 11, b1.max()[1] );

  b3 += Vector2(2,3);
  EXPECT_DOUBLE_EQ( 4, b3.min()[0] );
  EXPECT_DOUBLE_EQ( 7, b3.min()[1] );
  EXPECT_DOUBLE_EQ( 8, b3.max()[0] );
  EXPECT_DOUBLE_EQ( 11, b3.max()[1] );

  // Difference-assignment
  b1 -= Vector2(4,4);
  EXPECT_DOUBLE_EQ( 0, b1.min()[0] );
  EXPECT_DOUBLE_EQ( 3, b1.min()[1] );
  EXPECT_DOUBLE_EQ( 4, b1.max()[0] );
  EXPECT_DOUBLE_EQ( 7, b1.max()[1] );

  b3 -= Vector2(4,4);
  EXPECT_DOUBLE_EQ( 0, b3.min()[0] );
  EXPECT_DOUBLE_EQ( 3, b3.min()[1] );
  EXPECT_DOUBLE_EQ( 4, b3.max()[0] );
  EXPECT_DOUBLE_EQ( 7, b3.max()[1] );

  // Left scale
  b2 = 2*b1;
  EXPECT_DOUBLE_EQ( 0, b2.min()[0] );
  EXPECT_DOUBLE_EQ( 6, b2.min()[1] );
  EXPECT_DOUBLE_EQ( 8, b2.max()[0] );
  EXPECT_DOUBLE_EQ( 14, b2.max()[1] );

  b2 = 2*b3;
  EXPECT_DOUBLE_EQ( 0, b2.min()[0] );
  EXPECT_DOUBLE_EQ( 6, b2.min()[1] );
  EXPECT_DOUBLE_EQ( 8, b2.max()[0] );
  EXPECT_DOUBLE_EQ( 14, b2.max()[1] );

  // Right scale
  b2 = b1*3;
  EXPECT_DOUBLE_EQ( 0, b2.min()[0] );
  EXPECT_DOUBLE_EQ( 9, b2.min()[1] );
  EXPECT_DOUBLE_EQ( 12, b2.max()[0] );
  EXPECT_DOUBLE_EQ( 21, b2.max()[1] );

  b2 = b3*3;
  EXPECT_DOUBLE_EQ( 0, b2.min()[0] );
  EXPECT_DOUBLE_EQ( 9, b2.min()[1] );
  EXPECT_DOUBLE_EQ( 12, b2.max()[0] );
  EXPECT_DOUBLE_EQ( 21, b2.max()[1] );

  // Right quotient
  b2 = b1/2;
  EXPECT_DOUBLE_EQ( 0, b2.min()[0] );
  EXPECT_DOUBLE_EQ( 1.5, b2.min()[1] );
  EXPECT_DOUBLE_EQ( 2, b2.max()[0] );
  EXPECT_DOUBLE_EQ( 3.5, b2.max()[1] );

  b2 = b3/2;
  EXPECT_DOUBLE_EQ( 0, b2.min()[0] );
  EXPECT_DOUBLE_EQ( 1.5, b2.min()[1] );
  EXPECT_DOUBLE_EQ( 2, b2.max()[0] );
  EXPECT_DOUBLE_EQ( 3.5, b2.max()[1] );

  // Right vector sum
  b2 = b1 + Vector2(1,1);
  EXPECT_DOUBLE_EQ( 1, b2.min()[0] );
  EXPECT_DOUBLE_EQ( 4, b2.min()[1] );
  EXPECT_DOUBLE_EQ( 5, b2.max()[0] );
  EXPECT_DOUBLE_EQ( 8, b2.max()[1] );

  b2 = b3 + Vector2(1,1);
  EXPECT_DOUBLE_EQ( 1, b2.min()[0] );
  EXPECT_DOUBLE_EQ( 4, b2.min()[1] );
  EXPECT_DOUBLE_EQ( 5, b2.max()[0] );
  EXPECT_DOUBLE_EQ( 8, b2.max()[1] );

  // Left vector sum
  b2 = Vector2(1,1) + b1;
  EXPECT_DOUBLE_EQ( 1, b2.min()[0] );
  EXPECT_DOUBLE_EQ( 4, b2.min()[1] );
  EXPECT_DOUBLE_EQ( 5, b2.max()[0] );
  EXPECT_DOUBLE_EQ( 8, b2.max()[1] );

  b2 = Vector2(1,1) + b3;
  EXPECT_DOUBLE_EQ( 1, b2.min()[0] );
  EXPECT_DOUBLE_EQ( 4, b2.min()[1] );
  EXPECT_DOUBLE_EQ( 5, b2.max()[0] );
  EXPECT_DOUBLE_EQ( 8, b2.max()[1] );

  // Vector sum
  b2 = b1 - Vector2(1,1);
  EXPECT_DOUBLE_EQ( -1, b2.min()[0] );
  EXPECT_DOUBLE_EQ( 2, b2.min()[1] );
  EXPECT_DOUBLE_EQ( 3, b2.max()[0] );
  EXPECT_DOUBLE_EQ( 6, b2.max()[1] );

  b2 = b3 - Vector2(1,1);
  EXPECT_DOUBLE_EQ( -1, b2.min()[0] );
  EXPECT_DOUBLE_EQ( 2, b2.min()[1] );
  EXPECT_DOUBLE_EQ( 3, b2.max()[0] );
  EXPECT_DOUBLE_EQ( 6, b2.max()[1] );
}
