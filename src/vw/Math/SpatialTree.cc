#include <vw/Math/SpatialTree.h>

#include <iostream>			   // debugging
#include <stdexcept>
using namespace std;

using namespace vw;
using namespace math;

// A spatial tree is a hierarchical subdivision of n-dimensional space

// A SpatialTree node can be in 3 states:
//
// 1) leaf (bbox set, no quads, possibly polygons)
// 2) split (bbox set, quads (some may be leaves), possibly polygons)
// 3) unset (no bbox, no quads, no polygons)

// FIXME!!! To do this incrementally, if the polygon we're adding is
// too big for the top level square, we should add three sibling
// squares to the level and create a parent square that's as big as
// that level (need to pick a direction to grow that best matches new
// polygon). Need to do this until the top level square contains the
// polygon. Should add a parent pointer to facilitate things.

class PrintFunctor : public ApplyFunctor
{
 public:
  PrintFunctor() {}
  virtual bool operator()(GeomPrimitive *prim, int level = -1)
  {
    cout << " Min[" << prim->bounding_box().min() << "]"
	 << " Max[" << prim->bounding_box().max() << "]"
	 << endl;
    return true; // continue processing         
  }
};

//NOTE: this can only write a 2D projection (because VRML is 3D)
class VRMLFunctor : public ApplyFunctor
{
 public:
  VRMLFunctor()
  {
    m_outFP = 0;
    m_selectedLevel = -1;
    m_zSpacing = 0.5;
    m_color[0] = 1; m_color[1] = m_color[2] = 0;
  }
  VRMLFunctor(FILE *outFP)
  {
    m_outFP = outFP;
    m_selectedLevel = -1;
    m_zSpacing = 0.5;
    m_color[0] = 1; m_color[1] = m_color[2] = 0;
  }
  void SelectLevel(int level) { m_selectedLevel = level; }
  void ZSpacing(float spacing) { m_zSpacing = spacing; }
  void Color(float red, float green, float blue)
  {
    m_color[0] = red; m_color[1] = green; m_color[2] = blue;
  }
  virtual bool should_process(const SpatialTree::BBoxT &bbox, int level = -1) const
  {
    ProcessBBox(bbox, level);
    return true; // process box
  }
  virtual bool operator()(GeomPrimitive *prim, int level = -1)
  {
    ProcessBBox(prim->bounding_box(), level);
    return true; // continue processing
  }
 private:
  void ProcessBBox(const SpatialTree::BBoxT &bbox, int level = -1) const
  {
    if (m_outFP != 0)
    {
      if ((m_selectedLevel != -1) && (level != m_selectedLevel))
	return;
    
      if (level > -1)
      {
	float colors[8][3] = {{ 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 },
			      { 0.0, 0.0, 1.0 }, { 1.0, 0.0, 1.0 },
			      { 0.0, 1.0, 1.0 }, { 1.0, 1.0, 0.0 },
			      { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }};
	static const int numColors = 8;	   // r,g,b, m,c,y, b,w
	int index = (level % numColors);
	fprintf(m_outFP, "  BaseColor { rgb [ %f %f %f ] }\n",
		colors[index][0], colors[index][1], colors[index][2]);
      }
      else
      {
	fprintf(m_outFP, "  BaseColor { rgb [ %f %f %f ] }\n",
		m_color[0], m_color[1], m_color[2]);
      }
      float z;
      if (level > -1)
	z = -level * m_zSpacing;
      else
	z = 0.0;
      fprintf(m_outFP, "  FaceSet {\n");
      fprintf(m_outFP, "     vertexProperty VertexProperty {\n");
      fprintf(m_outFP, "	vertex [\n");
      fprintf(m_outFP, "		%f %f %f,\n",
	      bbox.min()[0], bbox.min()[1], z);
      fprintf(m_outFP, "		%f %f %f,\n",
	      bbox.max()[0], bbox.min()[1], z);
      fprintf(m_outFP, "		%f %f %f,\n",
	      bbox.max()[0], bbox.max()[1], z);
      fprintf(m_outFP, "		%f %f %f,\n",
	      bbox.min()[0], bbox.max()[1], z);
      fprintf(m_outFP, "	       ]\n");
      fprintf(m_outFP, "     }\n");
      fprintf(m_outFP, "     numVertices [ 4 ]\n");
      fprintf(m_outFP, "  }\n");
    }
    else
    {
      printf("ERROR in VRMLFunctor: no output file pointer set!\n");
    }
  }
 private:
  FILE *m_outFP;
  int m_selectedLevel;
  float m_zSpacing;
  float m_color[3];
};

class ContainsOneFunctor : public ApplyFunctor
{
 public:
  ContainsOneFunctor(const SpatialTree::VectorT *point) : m_point(point), m_prim(0) {}
  virtual bool process_children_first() const {return true;}
  virtual bool should_process(const BBox<double> &bbox, int level = -1) const
  {
    return bbox.contains(*m_point);
  }
  virtual bool operator()(GeomPrimitive *prim, int level = -1)
  {
    if (prim->bounding_box().contains(*m_point) && prim->contains(*m_point))
    {
      m_prim = prim;
      return false; // stop processing
    }
    return true; // continue processing
  }
  GeomPrimitive *get_primitive() {return m_prim;}
 private:
  const SpatialTree::VectorT *m_point;
  GeomPrimitive *m_prim;
};

GeomPrimitive *
SpatialTree::Contains(const VectorT &point) {
  ContainsOneFunctor func(&point);
  Apply(func);
  return func.get_primitive();
}

class ContainsAllFunctor : public ApplyFunctor
{
 public:
  ContainsAllFunctor(const SpatialTree::VectorT *point) : m_point(point), m_alloc(true)
  {
    m_prims = new list<GeomPrimitive*>;
  }
  ContainsAllFunctor(const SpatialTree::VectorT *point, list<GeomPrimitive*> *prims)
    : m_point(point), m_alloc(false), m_prims(prims) {}
  virtual ~ContainsAllFunctor()
  {
    if (m_alloc)
      delete m_prims;
  }
  virtual bool process_children_first() const {return true;}
  virtual bool should_process(const BBox<double> &bbox, int level = -1) const
  {
    return bbox.contains(*m_point);
  }
  virtual bool operator()(GeomPrimitive *prim, int level = -1)
  {
    if (prim->bounding_box().contains(*m_point) && prim->contains(*m_point))
      m_prims->push_back(prim);
    return true; // continue processing
  }
  list<GeomPrimitive*> *get_primitives() {return m_prims;}
 private:
  const SpatialTree::VectorT *m_point;
  bool m_alloc;
  list<GeomPrimitive*> *m_prims;
};

void
SpatialTree::ContainsAll(const VectorT &point, list<GeomPrimitive*> &prims) {
  ContainsAllFunctor func(&point, &prims);
  Apply(func);
}

#if 0
class AllOverlapsFunctor : public ApplyFunctor
{
 public:
  AllOverlapsFunctor(GeomPrimitive *overlapPrim) : m_overlapPrim(overlapPrim), m_alloc(true)
  {
    m_overlaps = new list<pair<GeomPrimitive*, GeomPrimitive*> >;
  }
  AllOverlapsFunctor(GeomPrimitive *overlapPrim,
                     list<pair<GeomPrimitive*, GeomPrimitive*> > *overlaps)
                       : m_overlapPrim(overlapPrim), m_alloc(false), m_overlaps(overlaps) {}                    
  virtual ~AllOverlapsFunctor()
  {
    if (m_alloc)
      delete m_overlaps;
  }
  virtual bool process_children_first() const {return true;}
  virtual bool should_process(const BBox<double> &bbox, int level = -1) const
  {
    return bbox.contains(*m_point);
  }
  virtual bool operator()(GeomPrimitive *prim, int level = -1)
  {
    if (prim != m_overlapPrim && prim->bounding_box().intersects(m_overlapPrim->bounding_box()))
      overlaps->push_back(make_pair(m_overlapPrim, prim));
    return true; // continue processing
  }
  list<pair<GeomPrimitive*, GeomPrimitive*> > *get_overlaps() {returm m_overlaps;}
 private:
  GeomPrimitive *m_overlapPrim;
  bool m_alloc;
  list<pair<GeomPrimitive*, GeomPrimitive*> > *m_overlaps;
};

void
SpatialTree::AllOverlaps(const VectorT &point, list<pair<GeomPrimitive*, GeomPrimitive*> > &overlaps) {
  ContainsAllFunctor func(&point, &overlaps);
  ContainsHelper(func);
}
#endif

#if 0
GeomPrimitive *
SpatialTree::Contains(const VectorT &point)
{
  // Check quadrant BBox to see if it contains the point... 
  if (m_BBox.contains(point))
  {
    if (IsSplit())
    {
      // No primitives at this level contained the point so recurse down
      // into the quadrant that contains the point
      // This is inneficient... we should do a breadth-first search
      for (int i = 0; i < m_numQuadrants; i++)
      {
	GeomPrimitive *geomPrimitive = m_quadrant[i]->Contains(point);

	if (geomPrimitive != 0)
	  return geomPrimitive;
      }
    }

    // If this quadrant has any associated primitives, see if one of
    // them contains the point
    if (m_primitiveList != 0)
    {
      PrimitiveListElem *listElem = m_primitiveList;
      while (listElem != 0)
      {
	if (listElem->geomPrimitive->bounding_box().contains(point))
	  if (listElem->geomPrimitive->contains(point))
	    return listElem->geomPrimitive;
	listElem = listElem->next;
      }
    }
  }

  return 0;
}
#endif

void
SpatialTree::Print()
{
  PrintFunctor func;
  Apply(func);
}

void
SpatialTree::WriteVRML(char *fileName, int level)
{
  FILE *outFP = 0;

  assert(m_numQuadrants >= 4);

  if ((outFP = fopen(fileName, "w" )) == 0)
  {
    fprintf(stderr, "WriteVRML(): cannot open output file: %s\n", fileName);
    exit(-1);
  }

  VRMLFunctor vrmlFunctor(outFP);
  vrmlFunctor.SelectLevel(level);

  WriteVRMLHead(outFP);
  Apply(vrmlFunctor);
  WriteVRMLTail(outFP);

  fclose(outFP);
}

bool
SpatialTree::Apply(ApplyFunctor &func, unsigned int level /* = 0 */)
{
  if (!func.should_process(m_BBox, level))
    return true; // continue processing one level above
    
  if (func.process_children_first() && IsSplit())
  {
    unsigned int next_level = level + 1;
    for (int i = 0; i < m_numQuadrants; i++)
    {
      if (!m_quadrant[i]->Apply(func, next_level))
        return false; // do not continue processing one level above
    }
  }
  
  if (m_primitiveList != 0)
  {
    PrimitiveListElem *listElem = m_primitiveList;
    while (listElem != 0)
    {
      if (!func(listElem->geomPrimitive, level))
        return false; // do not continue processing one level above
      listElem = listElem->next;
    }
  }

  if (!func.process_children_first() && IsSplit())
  {
    unsigned int next_level = level + 1;
    for (int i = 0; i < m_numQuadrants; i++)
    {
      if (!m_quadrant[i]->Apply(func, next_level))
        return false; // do not continue processing one level above
    }
  }
  
  return true; // continue processing one level above
}

void
SpatialTree::WriteVRMLHead(FILE *outFP)
{
  fprintf(outFP, "#VRML V1.0 ascii\n#\n");
  fprintf(outFP, "Separator\n");
  fprintf(outFP, "{\n");
  fprintf(outFP, "  DrawStyle { style   LINES lineWidth 1 }\n");
}

void
SpatialTree::WriteVRMLTail(FILE *outFP)
{
  fprintf(outFP, "}\n");
}


bool
SpatialTree::AddGeomPrimitive(GeomPrimitive *newPrim, BBoxT &newBBox)
{
  // Only add if it fits inside this levels bounding box...
  if (m_BBox.contains(newBBox))
  {
    if (IsSplit())
    {
      // Try to add the GeomPrimitive to a child quad:
      for (int i = 0; i < m_numQuadrants; i++)
	if (m_quadrant[i]->AddGeomPrimitive(newPrim, newBBox))
	  return true;
    }
    else					   // no quads
    {
      BBoxT quadrantBBoxes[m_numQuadrants];

      QuadSplit(quadrantBBoxes);
      // Check to see if new bbox would fit in a child quad... if so create
      // the quadrants
      for (int i = 0; i < m_numQuadrants; i++)
      {
	if (quadrantBBoxes[i].contains(newBBox)) // we need to make quads...
	{
	  Split(quadrantBBoxes);
	  // Add the GeomPrimitive to this child quad (we know it fits):
	  return m_quadrant[i]->AddGeomPrimitive(newPrim, newBBox);
	}
      }
    }

    // If the newPrim bbox wouldn't fit in any child quad we add it to
    // this Quad
    PrimitiveListElem *newListElem = new PrimitiveListElem;

    newListElem->geomPrimitive = newPrim;
    newListElem->next = m_primitiveList;
    m_primitiveList = newListElem;
    m_numPrimitives++;

    return true;
  }

  return false;
}

bool
SpatialTree::AddGeomPrimitive(GeomPrimitive *newPrim)
{
  BBoxT newBBox(newPrim->bounding_box()); // make copy of bbox for speed...

  return AddGeomPrimitive(newPrim, newBBox);
}
