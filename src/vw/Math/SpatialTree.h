#ifndef __VW_MATH_SPATIAL_TREE_H__
#define __VW_MATH_SPATIAL_TREE_H__

// Note that this class is not intended to implement a Binary
// Space Partition (BSP)...

// STL includes:
#include <iostream>
#include <list>

#include <vw/Math/Vector.h>
#include <vw/Math/BBox.h>

namespace vw {
namespace math {

  class GeomPrimitive {
  public:
    virtual bool contains(const Vector<double> &point) const = 0;
    virtual const BBox<double> &bounding_box() const = 0;
  };

  class SpatialTree {
  public:
    typedef BBox<double> BBoxT;
    typedef Vector<double> VectorT;

    struct PrimitiveListElem {
      PrimitiveListElem() {
        next = 0;
        prim = 0;
      }
      PrimitiveListElem *next;
      GeomPrimitive *prim;
    };
    
    struct SpatialTreeNode {
      SpatialTreeNode(int num_quadrants) {
        m_quadrant = new SpatialTreeNode*[num_quadrants];
        for (int i = 0; i < num_quadrants; i++)
          m_quadrant[i] = 0;
        m_primitive_list = 0;
        m_num_primitives = 0; // mostly for debugging
        m_is_split = false;
      }
      ~SpatialTreeNode() {
        delete[] m_quadrant;
        m_quadrant = 0;   
      }
      bool is_split() { return m_is_split; }
      BBoxT &bounding_box() { return m_bbox; }
      BBoxT m_bbox;
      PrimitiveListElem *m_primitive_list;
      int m_num_primitives; // mostly for debugging
      SpatialTreeNode **m_quadrant;
      bool m_is_split;
    };

    SpatialTree(BBoxT bbox);
    SpatialTree(int num_primitives, GeomPrimitive **prims);
    ~SpatialTree();
    void add(GeomPrimitive *prim);
    BBoxT &bounding_box() { return m_root_node->m_bbox; }
    GeomPrimitive *contains(const VectorT &point);
    void contains(const VectorT &point, std::list<GeomPrimitive*> &prims);
    void overlap_pairs(std::list<std::pair<GeomPrimitive*, GeomPrimitive*> > &overlaps);
    void print(std::ostream &os = std::cout);
    //NOTE: this can only write a 2D projection (because VRML is 3D)
    void write_vrml(char *fn, int level = -1);
    void write_vrml(std::ostream &os = std::cout, int level = -1);
  private:
    int m_dim;
    int m_num_quadrants;
    SpatialTreeNode *m_root_node;
  };

}} // namespace vw::math

#endif // _SPATIAL_TREE_H_
