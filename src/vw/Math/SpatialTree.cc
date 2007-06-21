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

class PrintFunctor : public BBoxFunctor
{
 public:
  PrintFunctor() {}
  virtual void operator()(const BBox2D &bbox, int level = -1) const
  {
    cout << " Min[" << bbox.Min().X() << "," << bbox.Min().Y() << "]"
	 << " Max[" << bbox.Max().X() << "," << bbox.Max().Y() << "]"
	 << endl;
  }
};

class VRMLFunctor : public BBoxFunctor
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
  virtual void operator()(const BBox2D &bbox, int level = -1) const
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
	      bbox.Min().X(), bbox.Min().Y(), z);
      fprintf(m_outFP, "		%f %f %f,\n",
	      bbox.Max().X(), bbox.Min().Y(), z);
      fprintf(m_outFP, "		%f %f %f,\n",
	      bbox.Max().X(), bbox.Max().Y(), z);
      fprintf(m_outFP, "		%f %f %f,\n",
	      bbox.Min().X(), bbox.Max().Y(), z);
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

GeomPrimitive *
SpatialTree::Contains(const Point2D &point)
{
  // Check quadrant BBox to see if it contains the point... 
  if (m_BBox.Contains(point))
  {
    if (IsSplit())
    {
      // No primitives at this level contained the point so recurse down
      // into the quadrant that contains the point
      // This is inneficient... we should do a breadth-first search
      for (int i = 0; i < 4; i++)
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
	if (listElem->geomPrimitive->BoundingBox().Contains(point))
	  if (listElem->geomPrimitive->Contains(point))
	    return listElem->geomPrimitive;
	listElem = listElem->next;
      }
    }
  }

  return 0;
}

void
SpatialTree::Print()
{
  Apply(PrintFunctor());
}

void
SpatialTree::WriteVRML(char *fileName, int level)
{
  FILE *outFP = 0;

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

void
SpatialTree::Apply(const BBoxFunctor &bboxFunctor)
{
  static unsigned int level = 0;

  bboxFunctor(m_BBox, level);

  if (m_primitiveList != 0)
  {
    PrimitiveListElem *listElem = m_primitiveList;
    while (listElem != 0)
    {
      BBox2D bbox = listElem->geomPrimitive->BoundingBox();
      bboxFunctor(bbox, level);
      listElem = listElem->next;
    }
  }

  if (IsSplit())
  {
    level++;
    for (int i = 0; i < 4; i++)
      m_quadrant[i]->Apply(bboxFunctor);
    level--;
  }
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
SpatialTree::AddGeomPrimitive(GeomPrimitive *newPrim, BBox2D &newBBox)
{
  // Only add if it fits inside this levels bounding box...
  if (m_BBox.Contains(newBBox))
  {
    if (IsSplit())
    {
      // Try to add the GeomPrimitive to a child quad:
      for (int i = 0; i < 4; i++)
	if (m_quadrant[i]->AddGeomPrimitive(newPrim, newBBox))
	  return true;
    }
    else					   // no quads
    {
      BBox2D quadrantBBoxes[4];

      m_BBox.QuadSplit(quadrantBBoxes);
      // Check to see if new bbox would fit in a child quad... if so create
      // the quadrants
      for (int i = 0; i < 4; i++)
      {
	if (quadrantBBoxes[i].Contains(newBBox)) // we need to make quads...
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
  BBox2D newBBox(newPrim->BoundingBox()); // make copy of bbox for speed...

  return AddGeomPrimitive(newPrim, newBBox);
}
