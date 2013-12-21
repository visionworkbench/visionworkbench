// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


/// \file DisjointSet.h
///
/// Provides a disjoint set data structure.
#ifndef __VW_MATH_DISJOINTSET_H__
#define __VW_MATH_DISJOINTSET_H__

// This is an efficient disjoint set data structure, i.e. a set of
// values is partitioned into multiple disjoint subsets. It allows the
// user to insert new elements as single-element sets, insert new
// elements into existing sets, combine sets, and find the set that
// contains a given element.

// Note that find(Elem) is much more efficient than find(ElemT), so it
// is strongly recommended that you store the Elem returned by
// insert() if you intend to find() it again later.  Also, it is possible
// to store multiple elements with the same value in a DisjointSet, but
// if you do so then find(ElemT) becomes ambiguous.

#include <vw/Core/Exception.h>

#include <list>

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
      else if ((s1->rank) < (s2->rank)) {
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
      for (n = e; n->parent; n = n->parent) ;
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

#endif // __VW_MATH_DISJOINTSET_H__
