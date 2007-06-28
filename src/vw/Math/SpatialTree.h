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
#include <list>
using namespace std;

#include <vw/Math/Vector.h>
#include <vw/Math/BBox.h>

namespace vw
{
  namespace math
  {
    class GeomPrimitive
    {
    public:
      virtual bool contains(const Vector<double> &point) const = 0;
      virtual const BBox<double> &bounding_box() const = 0;
    };

    class ApplyFunctor
    {
    public:
      virtual ~ApplyFunctor() {}
      virtual bool process_children_first() const {return false; }
      virtual bool should_process(const BBox<double> &bbox, int level = -1) const {return true;}
      virtual bool operator()(GeomPrimitive *prim, int level = -1) {return true;}
    };

    class SpatialTree
    {
    public:
      typedef BBox<double> BBoxT;
      typedef Vector<double> VectorT;

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

      enum SizeType
      {
	SIZE_DIMENSION,
	SIZE_DIMENSION_SQUARED
      };

      SpatialTree(int size, SizeType size_type = SIZE_DIMENSION)
      {
	if (size_type == SIZE_DIMENSION_SQUARED)
	  m_numQuadrants = size;
	else
	  m_numQuadrants = (int)((unsigned)1 << (unsigned)size);
	m_quadrant = new SpatialTree*[m_numQuadrants];
	for (int i = 0; i < m_numQuadrants; i++)
	  m_quadrant[i] = 0;
	m_primitiveList = 0;
	m_numPrimitives = 0;		   // mostly for debugging
	m_isSplit = false;
      }
      SpatialTree(int numPrimitives, GeomPrimitive **geomPrimitives)
      {
	VW_ASSERT( geomPrimitives != 0 && geomPrimitives[0] != 0, ArgumentErr() << "No GeomPrimitives provided." );
	unsigned size = geomPrimitives[0]->bounding_box().min().size();
	m_numQuadrants = (int)((unsigned)1 << size);
	m_quadrant = new SpatialTree*[m_numQuadrants];
	for (int i = 0; i < m_numQuadrants; i++)
	  m_quadrant[i] = 0;
	m_primitiveList = 0;
	m_numPrimitives = 0;		   // mostly for debugging
	m_isSplit = false;
	GrowBBox(numPrimitives, geomPrimitives);

	cout << "Bounding Box of Total mesh ="
	     << " Min[" << m_BBox.min() << "]"
	     << " Max[" << m_BBox.max() << "]"
	     << endl;

	for (int i = 0; i < numPrimitives; i++)
	  AddGeomPrimitive(geomPrimitives[i]);
      }
      ~SpatialTree()
      {
	for (int i = 0; i < m_numQuadrants; i++)
	{
	  delete m_quadrant[i];
	  m_quadrant[i] = 0;
	}
	delete[] m_quadrant;
	m_quadrant = 0;
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
      BBoxT &BoundingBox() { return m_BBox; }
      GeomPrimitive *Contains(const VectorT &point);
      void ContainsAll(const VectorT &point, list<GeomPrimitive*> &prims);
      void Print();
      //NOTE: this can only write a 2D projection (because VRML is 3D)
      void WriteVRML(char *fileName, int level = -1);
    private:
      bool Apply(ApplyFunctor &func, unsigned int level = 0);
      void WriteVRMLHead(FILE *outFP);
      void WriteVRMLTail(FILE *outFP);
      void WriteVRMLBox(FILE *outFP, const BBoxT &bbox);

      bool AddGeomPrimitive(GeomPrimitive *newPrim);
      bool AddGeomPrimitive(GeomPrimitive *newPrim, BBoxT &newBBox);
      bool IsSplit() { return m_isSplit; }
      void Split(BBoxT quadrantBBoxes[])
      {
	m_isSplit = true;
	for (int i = 0; i < m_numQuadrants; i++)
	{
	  m_quadrant[i] = new SpatialTree(m_numQuadrants, SIZE_DIMENSION_SQUARED);
	  m_quadrant[i]->m_BBox = quadrantBBoxes[i];
	}
      }
      void QuadSplit(BBoxT quadrantBBoxes[])
      {
	using namespace vw::math::vector_containment_comparison;
	VectorT center = m_BBox.center();
	VectorT diagonalVec = m_BBox.size() * 0.5;
	VectorT axisVec;
	int size = center.size();

	for (int i = 0; i < m_numQuadrants; i++)
	  quadrantBBoxes[i].min() = center;
	fill(axisVec, 0.0);
	for (int d = 0, width = m_numQuadrants, halfWidth = m_numQuadrants / 2; d < size; d++, width = halfWidth, halfWidth /= 2)
	{
	  assert(d < size - 1 || width > 1);
	  axisVec[d] = diagonalVec[d];
	  for (int i = 0; i < size;)
	  {
	    for (; i < halfWidth; i++)
	      quadrantBBoxes[i].min() -= axisVec;
	    for (; i < width; i++)
	      quadrantBBoxes[i].min() += axisVec;
	  }
	  axisVec[d] = 0.0;
	}
	for (int i = 0; i < m_numQuadrants; i++)
	{
	  quadrantBBoxes[i].max() = max(center, quadrantBBoxes[i].min());
	  quadrantBBoxes[i].min() = min(center, quadrantBBoxes[i].min()); //NOTE: make sure that this is ok
	}

// Debugging:
	if (m_numQuadrants == 4)
	{
	  assert(quadrantBBoxes[0].max() >= quadrantBBoxes[0].min());
	  assert(quadrantBBoxes[1].max() >= quadrantBBoxes[1].min());
	  assert(quadrantBBoxes[2].max() >= quadrantBBoxes[2].min());
	  assert(quadrantBBoxes[3].max() >= quadrantBBoxes[3].min());
	  assert(quadrantBBoxes[3].min() >= quadrantBBoxes[0].min());
	  assert(quadrantBBoxes[3].max() >= quadrantBBoxes[0].max());
	  assert(quadrantBBoxes[0].max() == quadrantBBoxes[3].min());
	}
      }
      void GrowBBox(int numPrimitives, GeomPrimitive *geomPrimitives[])
      {
	cout << "SpatialTree: growing BBox..." << flush;
	for (int i = 0; i < numPrimitives; i++)
	  m_BBox.grow(geomPrimitives[i]->bounding_box());
	cout << "done." << endl;
      }
      int m_numQuadrants;
      BBoxT m_BBox;
      PrimitiveListElem *m_primitiveList;
      int m_numPrimitives;		   // mostly for debugging
      SpatialTree **m_quadrant;
      bool m_isSplit;
    };
  }
}

#endif	// _SPATIAL_TREE_H_
