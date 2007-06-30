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

#include <sstream>

#include <cxxtest/TestSuite.h>
#include <vw/Math/Vector.h>
#include <vw/Math/BBox.h>
#include <vw/Math/SpatialTree.h>

using namespace std;
using namespace vw;
using namespace vw::math;

#define TEST_SPATIAL_TREE_PRINT_RESULT "+ Min[Vector2(0,0)] Max[Vector2(16,16)]\n  Min[Vector2(9,9)] Max[Vector2(9,9)]\n  + Min[Vector2(8,8)] Max[Vector2(16,8)]\n  + Min[Vector2(8,8)] Max[Vector2(16,8)]\n  + Min[Vector2(0,8)] Max[Vector2(8,16)]\n  + Min[Vector2(0,0)] Max[Vector2(8,8)]\n    Min[Vector2(1.5,3)] Max[Vector2(2,5)]\n    + Min[Vector2(4,4)] Max[Vector2(8,4)]\n    + Min[Vector2(4,4)] Max[Vector2(8,4)]\n    + Min[Vector2(0,4)] Max[Vector2(4,8)]\n    + Min[Vector2(0,0)] Max[Vector2(4,4)]\n      + Min[Vector2(2,2)] Max[Vector2(4,2)]\n      + Min[Vector2(2,2)] Max[Vector2(4,2)]\n      + Min[Vector2(0,2)] Max[Vector2(2,4)]\n        Min[Vector2(1,2)] Max[Vector2(1.75,4)]\n      + Min[Vector2(0,0)] Max[Vector2(2,2)]\n        + Min[Vector2(1,1)] Max[Vector2(2,1)]\n        + Min[Vector2(1,1)] Max[Vector2(2,1)]\n        + Min[Vector2(0,1)] Max[Vector2(1,2)]\n        + Min[Vector2(0,0)] Max[Vector2(1,1)]\n          + Min[Vector2(0.5,0.5)] Max[Vector2(1,0.5)]\n          + Min[Vector2(0.5,0.5)] Max[Vector2(1,0.5)]\n          + Min[Vector2(0,0.5)] Max[Vector2(0.5,1)]\n          + Min[Vector2(0,0)] Max[Vector2(0.5,0.5)]\n            + Min[Vector2(0.25,0.25)] Max[Vector2(0.5,0.25)]\n            + Min[Vector2(0.25,0.25)] Max[Vector2(0.5,0.25)]\n            + Min[Vector2(0,0.25)] Max[Vector2(0.25,0.5)]\n            + Min[Vector2(0,0)] Max[Vector2(0.25,0.25)]\n              + Min[Vector2(0.125,0.125)] Max[Vector2(0.25,0.125)]\n              + Min[Vector2(0.125,0.125)] Max[Vector2(0.25,0.125)]\n              + Min[Vector2(0,0.125)] Max[Vector2(0.125,0.25)]\n              + Min[Vector2(0,0)] Max[Vector2(0.125,0.125)]\n                Min[Vector2(0.1,0.1)] Max[Vector2(0.1,0.1)]\n"

#define TEST_SPATIAL_TREE_VRML_RESULT "#VRML V1.0 ascii\n#\nSeparator {\n  DrawStyle { style   LINES lineWidth 1 }\n  BaseColor { rgb [ 1 0 0 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         0 0 0,\n         16 0 0,\n         16 16 0,\n         0 16 0,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 1 0 0 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         9 9 0,\n         9 9 0,\n         9 9 0,\n         9 9 0,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 0 1 0 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         8 8 -0.5,\n         16 8 -0.5,\n         16 8 -0.5,\n         8 8 -0.5,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 0 1 0 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         8 8 -0.5,\n         16 8 -0.5,\n         16 8 -0.5,\n         8 8 -0.5,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 0 1 0 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         0 8 -0.5,\n         8 8 -0.5,\n         8 16 -0.5,\n         0 16 -0.5,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 0 1 0 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         0 0 -0.5,\n         8 0 -0.5,\n         8 8 -0.5,\n         0 8 -0.5,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 0 1 0 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         1.5 3 -0.5,\n         2 3 -0.5,\n         2 5 -0.5,\n         1.5 5 -0.5,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 0 0 1 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         4 4 -1,\n         8 4 -1,\n         8 4 -1,\n         4 4 -1,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 0 0 1 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         4 4 -1,\n         8 4 -1,\n         8 4 -1,\n         4 4 -1,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 0 0 1 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         0 4 -1,\n         4 4 -1,\n         4 8 -1,\n         0 8 -1,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 0 0 1 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         0 0 -1,\n         4 0 -1,\n         4 4 -1,\n         0 4 -1,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 1 0 1 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         2 2 -1.5,\n         4 2 -1.5,\n         4 2 -1.5,\n         2 2 -1.5,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 1 0 1 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         2 2 -1.5,\n         4 2 -1.5,\n         4 2 -1.5,\n         2 2 -1.5,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 1 0 1 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         0 2 -1.5,\n         2 2 -1.5,\n         2 4 -1.5,\n         0 4 -1.5,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 1 0 1 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         1 2 -1.5,\n         1.75 2 -1.5,\n         1.75 4 -1.5,\n         1 4 -1.5,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 1 0 1 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         0 0 -1.5,\n         2 0 -1.5,\n         2 2 -1.5,\n         0 2 -1.5,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 0 1 1 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         1 1 -2,\n         2 1 -2,\n         2 1 -2,\n         1 1 -2,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 0 1 1 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         1 1 -2,\n         2 1 -2,\n         2 1 -2,\n         1 1 -2,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 0 1 1 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         0 1 -2,\n         1 1 -2,\n         1 2 -2,\n         0 2 -2,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 0 1 1 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         0 0 -2,\n         1 0 -2,\n         1 1 -2,\n         0 1 -2,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 1 1 0 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         0.5 0.5 -2.5,\n         1 0.5 -2.5,\n         1 0.5 -2.5,\n         0.5 0.5 -2.5,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 1 1 0 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         0.5 0.5 -2.5,\n         1 0.5 -2.5,\n         1 0.5 -2.5,\n         0.5 0.5 -2.5,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 1 1 0 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         0 0.5 -2.5,\n         0.5 0.5 -2.5,\n         0.5 1 -2.5,\n         0 1 -2.5,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 1 1 0 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         0 0 -2.5,\n         0.5 0 -2.5,\n         0.5 0.5 -2.5,\n         0 0.5 -2.5,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 0 0 0 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         0.25 0.25 -3,\n         0.5 0.25 -3,\n         0.5 0.25 -3,\n         0.25 0.25 -3,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 0 0 0 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         0.25 0.25 -3,\n         0.5 0.25 -3,\n         0.5 0.25 -3,\n         0.25 0.25 -3,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 0 0 0 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         0 0.25 -3,\n         0.25 0.25 -3,\n         0.25 0.5 -3,\n         0 0.5 -3,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 0 0 0 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         0 0 -3,\n         0.25 0 -3,\n         0.25 0.25 -3,\n         0 0.25 -3,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 1 1 1 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         0.125 0.125 -3.5,\n         0.25 0.125 -3.5,\n         0.25 0.125 -3.5,\n         0.125 0.125 -3.5,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 1 1 1 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         0.125 0.125 -3.5,\n         0.25 0.125 -3.5,\n         0.25 0.125 -3.5,\n         0.125 0.125 -3.5,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 1 1 1 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         0 0.125 -3.5,\n         0.125 0.125 -3.5,\n         0.125 0.25 -3.5,\n         0 0.25 -3.5,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 1 1 1 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         0 0 -3.5,\n         0.125 0 -3.5,\n         0.125 0.125 -3.5,\n         0 0.125 -3.5,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n  BaseColor { rgb [ 1 1 1 ] }\n  FaceSet {\n    vertexProperty VertexProperty {\n      vertex [\n         0.1 0.1 -3.5,\n         0.1 0.1 -3.5,\n         0.1 0.1 -3.5,\n         0.1 0.1 -3.5,\n      ]\n    }\n    numVertices [ 4 ]\n  }\n}\n"


class TestGeomPrimitive : public BBoxN, public GeomPrimitive
{
  bool contains(const Vector<double> &point) const {return BBoxN::contains(point);}
  const BBox<double> &bounding_box() const {return *this;}
};

class TestSpatialTree : public CxxTest::TestSuite
{
public:

  void test_spatial_tree() 
  {
    BBoxN b(0, 0, 1, 1);
    SpatialTree t(b);
    std::list<GeomPrimitive*> l;
    std::list<GeomPrimitive*>::iterator i;

    Vector<double,2> p0(0.1, 0.1);
    TestGeomPrimitive g0;
    g0.grow(p0);
    t.add(&g0);

    Vector<double,2> p1(1, 2), p2(1.75, 4);
    TestGeomPrimitive g1;
    g1.grow(p1);
    g1.grow(p2);
    t.add(&g1);

    Vector<double,2> p3(1.5, 3), p4(2, 5);
    TS_ASSERT_EQUALS( t.contains(p3), &g1 );
    l.clear();
    t.contains(p3, l);
    TS_ASSERT_EQUALS( l.size(), 1 );
    TS_ASSERT_EQUALS( *(l.begin()), &g1 );
    TS_ASSERT_EQUALS( t.contains(p4), (GeomPrimitive*)0 );
    l.clear();
    t.contains(p4, l);
    TS_ASSERT_EQUALS( l.size(), 0 );

    TestGeomPrimitive g2;
    g2.grow(p3);
    g2.grow(p4);
    t.add(&g2);

    Vector<double, 2> p5(1.25, 3.5), p6(1.6, 3.5), p7(1.75, 4.5), p8(1.25, 4.5), p9(8, 8);
    GeomPrimitive *gt1;
    GeomPrimitive *gt2;
    TS_ASSERT_EQUALS( t.contains(p5), &g1 );
    l.clear();
    t.contains(p5, l);
    TS_ASSERT_EQUALS( l.size(), 1 );
    TS_ASSERT_EQUALS( *(l.begin()), &g1 );
    gt1 = t.contains(p6);
    TS_ASSERT( gt1 == &g1 || gt1 == &g2 );
    l.clear();
    t.contains(p6, l);
    TS_ASSERT_EQUALS( l.size(), 2 );
    i = l.begin(); gt1 = *i; i++; gt2 = *i;
    TS_ASSERT( gt1 == &g1 || gt1 == &g2 );
    TS_ASSERT( gt2 == &g1 || gt2 == &g2 );
    TS_ASSERT( gt1 != gt2 );
    TS_ASSERT_EQUALS( t.contains(p7), &g2 );
    l.clear();
    t.contains(p7, l);
    TS_ASSERT_EQUALS( l.size(), 1 );
    TS_ASSERT_EQUALS( *(l.begin()), &g2 );
    TS_ASSERT_EQUALS( t.contains(p8), (GeomPrimitive*)0 );
    l.clear();
    t.contains(p8, l);
    TS_ASSERT_EQUALS( l.size(), 0 );
    TS_ASSERT_EQUALS( t.contains(p9), (GeomPrimitive*)0 );
    l.clear();
    t.contains(p9, l);
    TS_ASSERT_EQUALS( l.size(), 0 );

    Vector<double,2> p10(9, 9);
    TestGeomPrimitive g3;
    g3.grow(p10);
    t.add(&g3);

    ostringstream os;
    std::string printstr(TEST_SPATIAL_TREE_PRINT_RESULT);
    t.print(os);
    TS_ASSERT_EQUALS( os.str(), printstr );

    ostringstream os2;
    std::string vrmlstr(TEST_SPATIAL_TREE_VRML_RESULT);
    t.write_vrml(os2);
    TS_ASSERT_EQUALS( os2.str(), vrmlstr );
  }

}; // class TestSpatialTree
