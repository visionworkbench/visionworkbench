#ifndef __VW_MATH_SPATIAL_TREE_H__
#define __VW_MATH_SPATIAL_TREE_H__

// NOTE: the SpatialTree class is intended to be an N-Dimensional
// heirarchical geometry representation. Right now it only deals with
// 2-D geometry!

// Also note that this class is not intended to implement a Binary
// Space Partition (BSP)...

#ifdef __APPLE__
#include <float.h>			   // for DBL_MAX
#else
#include <values.h>			   // for DBL_MAX
#endif

#include <assert.h>

// STL includes:
#include <iostream>			   // debugging
using namespace std;

namespace vw
{
  namespace math
  {
    class BBox2D;

    class GeomPrimitive
    {
    public:
      virtual bool Contains(const class Point2D &point) const = 0;
      virtual const BBox2D &BoundingBox() const = 0;
    };

    class Point2D
    {
    public:
      Point2D() { m_x = m_y = 0.0; }
      Point2D(const double x, const double y = 0.0) { m_x = x; m_y = y; }
      Point2D(const Point2D &point) { m_x = point.m_x; m_y = point.m_y; }
      void X(const double x) { m_x = x; }
      void Y(const double y) { m_y = y; }
      double X() const { return m_x; }
      double Y() const { return m_y; }
      void AssignMax(const Point2D &point1, const Point2D &point2)
      {
	m_x = (point1.m_x > point2.m_x) ? point1.m_x : point2.m_x;
	m_y = (point1.m_y > point2.m_y) ? point1.m_y : point2.m_y;
      }
      void AssignMin(const Point2D &point1, const Point2D &point2)
      {
	m_x = (point1.m_x < point2.m_x) ? point1.m_x : point2.m_x;
	m_y = (point1.m_y < point2.m_y) ? point1.m_y : point2.m_y;
      }

      Point2D operator+(const Point2D &point) const
      {
	return Point2D(m_x + point.m_x, m_y + point.m_y);
      }

      Point2D operator-(const Point2D &point) const
      {
	return Point2D(m_x - point.m_x, m_y - point.m_y);
      }

      Point2D operator*(const double scalar) const
      {
	return Point2D(m_x*scalar, m_y * scalar);
      }

      bool operator==(const Point2D &point) const
      {
	return ((m_x == point.m_x) && (m_y == point.m_y));
      }

      bool operator<(const Point2D &point) const
      {
	return ((m_x < point.m_x) && (m_y < point.m_y));
      }

      bool operator<=(const Point2D &point) const
      {
	return ((m_x <= point.m_x) && (m_y <= point.m_y));
      }

      bool operator>(const Point2D &point) const
      {
	return ((m_x > point.m_x) && (m_y > point.m_y));
      }

      bool operator>=(const Point2D &point) const
      {
	return ((m_x >= point.m_x) && (m_y >= point.m_y));
      }

    private:
      double m_x, m_y;
    };

    class BBox2D
    {
    public:
      BBox2D() : m_min(DBL_MAX, DBL_MAX), m_max(-DBL_MAX, -DBL_MAX) {}
      BBox2D(const BBox2D &bbox) { m_min = bbox.Min(); m_max = bbox.Max(); }
      BBox2D(Point2D &min, Point2D &max) { m_min = min, m_max = max; }
      void Grow(const Point2D &point)
      {
	if (point.X() > m_max.X())
	  m_max.X(point.X());
	if (point.X() < m_min.X())
	  m_min.X(point.X());

	if (point.Y() > m_max.Y())
	  m_max.Y(point.Y());
	if (point.Y() < m_min.Y())
	  m_min.Y(point.Y());
      }
      void Grow(const BBox2D bbox)
      {
	Grow(bbox.Min());
 	Grow(bbox.Max());
      }
      bool Contains(const Point2D &point) const
      {
	return ((point >= m_min) && (point <= m_max));
      }
      bool Contains(const BBox2D &bbox) const
      {
	return ((bbox.m_min >= m_min) && (bbox.m_max <= m_max));
      }
      const Point2D &Min() const { return m_min; }
      const Point2D &Max() const { return m_max; }
      void Min(Point2D min) { m_min = min; }
      void Max(Point2D max) { m_max = max; }
      void QuadSplit(BBox2D quadrantBBoxes[])
      {
	Point2D center = (m_min + m_max) * 0.5;
	Point2D diagonalVec = (m_max - m_min) * 0.5;
	Point2D xVec(diagonalVec.X(), 0.0);
	Point2D yVec(0.0, diagonalVec.Y());

	// We arbitrarily start in lower left and go clockwise
	quadrantBBoxes[0].Min(m_min);
	quadrantBBoxes[0].Max(center);

	quadrantBBoxes[1].Min(m_min + yVec);
	quadrantBBoxes[1].Max(center + yVec);

	quadrantBBoxes[2].Min(center);
	quadrantBBoxes[2].Max(m_max);

	quadrantBBoxes[3].Min(m_min + xVec);
	quadrantBBoxes[3].Max(center + xVec);

// Debugging:
	assert(quadrantBBoxes[0].Max() >= quadrantBBoxes[0].Min());
	assert(quadrantBBoxes[1].Max() >= quadrantBBoxes[1].Min());
	assert(quadrantBBoxes[2].Max() >= quadrantBBoxes[2].Min());
	assert(quadrantBBoxes[3].Max() >= quadrantBBoxes[3].Min());
	assert(quadrantBBoxes[2].Min() >= quadrantBBoxes[0].Min());
	assert(quadrantBBoxes[2].Max() >= quadrantBBoxes[0].Max());
	assert(quadrantBBoxes[0].Max() == quadrantBBoxes[2].Min());
      }
    private:
      Point2D m_min, m_max;
    };

    struct PrimitiveListElem
    {
      PrimitiveListElem()
      {
	next = 0;
	geomPrimitive = 0;
      }
      PrimitiveListElem *next;
      GeomPrimitive *geomPrimitive;
    };

    class BBoxFunctor
    {
    public:
      virtual void operator()(const BBox2D &bbox, int level = -1) const {}
    };

    class SpatialTree
    {
    public:
      SpatialTree()
      {
	for (int i = 0; i < 4; i++)
	  m_quadrant[i] = 0;
	m_primitiveList = 0;
	m_numPrimitives = 0;		   // mostly for debugging
	m_isSplit = false;
      }
      SpatialTree(int numPrimitives, GeomPrimitive **geomPrimitives)
      {
	for (int i = 0; i < 4; i++)
	  m_quadrant[i] = 0;
	m_primitiveList = 0;
	m_numPrimitives = 0;		   // mostly for debugging
	m_isSplit = false;
	GrowBBox(numPrimitives, geomPrimitives);

	cout << "Bounding Box of Total mesh ="
	     << " Min[" << m_BBox.Min().X() << "," << m_BBox.Min().Y() << "]"
	     << " Max[" << m_BBox.Max().X() << "," << m_BBox.Max().Y() << "]"
	     << endl;

	for (int i = 0; i < numPrimitives; i++)
	  AddGeomPrimitive(geomPrimitives[i]);
      }
      ~SpatialTree()
      {
	for (int i = 0; i < 4; i++)
	{
	  delete m_quadrant[i];
	  m_quadrant[i] = 0;
	}
	while (m_primitiveList != 0)
	{
	  PrimitiveListElem *next = m_primitiveList->next;
	  delete m_primitiveList;
	  m_primitiveList = next;
	}
	m_primitiveList = 0;
	m_numPrimitives = 0;		   // mostly for debugging
	m_isSplit = false;
      }
      BBox2D &BoundingBox() { return m_BBox; }
      GeomPrimitive *Contains(const Point2D &point);
      void Print();
      void WriteVRML(char *fileName, int level = -1);
    private:
      void Apply(const BBoxFunctor &bboxFunctor);
      void WriteVRMLHead(FILE *outFP);
      void WriteVRMLTail(FILE *outFP);
      void WriteVRMLBox(FILE *outFP, const BBox2D &bbox);

      bool AddGeomPrimitive(GeomPrimitive *newPrim);
      bool AddGeomPrimitive(GeomPrimitive *newPrim, BBox2D &newBBox);
      bool IsSplit() { return m_isSplit; }
      void Split(BBox2D quadrantBBoxes[])
      {
	m_isSplit = true;
	for (int i = 0; i < 4; i++)
	{
	  m_quadrant[i] = new SpatialTree;
	  m_quadrant[i]->m_BBox = quadrantBBoxes[i];
	}
      }
      void GrowBBox(int numPrimitives, GeomPrimitive *geomPrimitives[])
      {
	cout << "SpatialTree: growing BBox..." << flush;
	for (int i = 0; i < numPrimitives; i++)
	  m_BBox.Grow(geomPrimitives[i]->BoundingBox());
	cout << "done." << endl;
      }
      BBox2D m_BBox;
      PrimitiveListElem *m_primitiveList;
      int m_numPrimitives;		   // mostly for debugging
      SpatialTree *m_quadrant[4];
      bool m_isSplit;
    };
  }
}

#endif	// _SPATIAL_TREE_H_
