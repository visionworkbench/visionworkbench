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

/// \file DisjointSet.h
///
/// Provides a disjoint set data structure.
#ifndef __VW_MATH__DISJOINT_SET_H__
#define __VW_MATH__DISJOINT_SET_H__

// This is an efficient disjoint set data structure. It allows the user
// to insert new elements as single-element sets, insert new elements into
// existing sets, combine sets, and find the set that contains a given
// element.

// Note that find(Elem) is much more efficient than find(ElemT), so it is
// strongly recommended that you store the Elem returned by insert().


#include <list>

#include <vw/Core/Exception.h>


namespace vw {
namespace math {

  template <class ElemT>
  class DisjointSet {
  private:
    struct ElemNode {
      ElemNode(const ElemT &elem_) : elem(elem_), parent(0), rank(0) {}
      ElemNode(const ElemT &elem_, ElemNode *parent_) : elem(elem_), parent(parent_), rank(0) {}
      ~ElemNode() {}
      ElemT elem;
      ElemNode *parent;
      unsigned rank;
    };
    
  public:
    typedef ElemNode* Elem;
    typedef ElemNode* Set;
    
    DisjointSet() {}
    ~DisjointSet() {
      typename std::list<ElemNode*>::iterator i;
      for (i = elems.begin(); i != elems.end(); i++)
        delete *i;
      num_elems = 0;
    }
  
    Elem insert(const ElemT &e) {
      ElemNode *n = new ElemNode(e);
      elems.push_back(n);
      num_elems++;
      return n;
    }
    Elem insert(const ElemT &e, Set s) {
      VW_ASSERT(s && !s->parent, ArgumentErr() << "Not a valid set!");
      ElemNode *n = new ElemNode(e, s);
      elems.push_back(n);
      num_elems++;
      return n;
    }
    
    Set combine(Set s1, Set s2) {
      VW_ASSERT(s1 && !s1->parent, ArgumentErr() << "Not a valid set!");
      VW_ASSERT(s2 && !s2->parent, ArgumentErr() << "Not a valid set!");
      if (s1 == s2)
        return s1;
      if (s1->rank > s2->rank) {
        s2->parent = s1;
        return s1;
      }
      if (s1->rank < s2->rank) {
        s1->parent = s2;
        return s2;
      }
      s2->parent = s1;
      s1->rank++;
      return s1;
    }
    
    Set find(Elem e) {
      VW_ASSERT(e, ArgumentErr() << "Not a valid element!");
      ElemNode *n;
      ElemNode *s;
      ElemNode *p;
      for (n = e; n->parent; n = n->parent);
      s = n;
      for (n = e, p = e->parent; p && p != s; n = p) {
        p = n->parent;
        n->parent = s;
      }
      return s;
    }
    Set find(const ElemT &e) {
      ElemNode *n = 0;
      typename std::list<ElemNode*>::iterator i;
      for (i = elems.begin(); i != elems.end(); i++) {
        if ((*i)->elem == e) {
          n = *i;
          break;
        }
      }
      if (!n)
        return 0;
      return find(n);
    }
    
  private:
    std::list<ElemNode*> elems;
    unsigned num_elems;
  };

}} // namespace vw::math

#endif // __VW_MATH__DISJOINT_SET_H__
