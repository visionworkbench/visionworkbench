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

#include <vw/Core/Exception.h>
#include <vw/Math/DisjointSet.h>
#include <vw/Math/MinimumSpanningTree.h>

#include <cstdlib>
#include <algorithm>

namespace {

  int edge_primitive_compare(const void *a_, const void *b_)
  {
    const vw::math::EdgePrimitive **a = (const vw::math::EdgePrimitive**)a_;
    const vw::math::EdgePrimitive **b = (const vw::math::EdgePrimitive**)b_;
    if ((*a)->cost() == (*b)->cost())
      return 0;
    else if ((*a)->cost() < (*b)->cost())
      return -1;
    return 1;
  }

} // namespace


namespace vw {
namespace math {

  MinimumSpanningTree::MinimumSpanningTree(int num_primitives_, EdgePrimitive **prims_) {
    VW_ASSERT(num_primitives > 0, ArgumentErr() << "No primitives provided!");
    int max_node;
    DisjointSet<int> ds;
    DisjointSet<int>::Elem *ds_elem;
    DisjointSet<int>::Set set1, set2;
    int i, j1, j2, k;

    num_primitives = num_primitives_;
    prims = new EdgePrimitive*[num_primitives];
    for (i = 0; i < num_primitives; i++)
      prims[i] = prims_[i];
    qsort(prims, num_primitives, sizeof(EdgePrimitive*), edge_primitive_compare);

    min_node = std::min(prims[0]->node1(), prims[0]->node2());
    max_node = std::max(prims[0]->node1(), prims[0]->node2());
    for (i = 1; i < num_primitives; i++) {
      min_node = std::min(min_node, prims[i]->node1());
      min_node = std::min(min_node, prims[i]->node2());
      max_node = std::max(max_node, prims[i]->node1());
      max_node = std::max(max_node, prims[i]->node2());
    }
    num_nodes = max_node - min_node + 1;

    node_used = new bool[num_nodes];
    for (i = 0; i < num_nodes; i++)
      node_used[i] = false;
    for (i = 0; i < num_primitives; i++) {
      j1 = prims[i]->node1() - min_node;
      j2 = prims[i]->node2() - min_node;
      node_used[j1] = true;
      node_used[j2] = true;
    }

    prim_used = new bool[num_primitives];
    for (i = 0; i < num_primitives; i++)
      prim_used[i] = false;
    num_edges = new int[num_nodes];
    ds_elem = new DisjointSet<int>::Elem[num_nodes];
    for (i = 0, k = min_node; i < num_nodes; i++, k++) {
      num_edges[i] = 0;
      if (node_used[i])
        ds_elem[i] = ds.insert(k);
      else
        ds_elem[i] = 0;
    }

    for (i = 0; i < num_primitives; i++) {
      j1 = prims[i]->node1() - min_node;
      j2 = prims[i]->node2() - min_node;
      set1 = ds.find(ds_elem[j1]);
      set2 = ds.find(ds_elem[j2]);
      if (set1 != set2) {
        ds.combine(set1, set2);
        prim_used[i] = true;
        num_edges[j1]++;
        num_edges[j2]++;
      }
    }

    edges = new EdgePrimitive**[num_nodes];
    for (i = 0; i < num_nodes; i++) {
      edges[i] = new EdgePrimitive*[num_edges[i]];
      num_edges[i] = 0;
    }
    for (i = 0; i < num_primitives; i++) {
      if (prim_used[i]) {
        j1 = prims[i]->node1() - min_node;
        j2 = prims[i]->node2() - min_node;
        edges[j1][num_edges[j1]] = prims[i];
        num_edges[j1]++;
        edges[j2][num_edges[j2]] = prims[i];
        num_edges[j2]++;
      }
    }

    delete[] ds_elem;
  }

  MinimumSpanningTree::~MinimumSpanningTree() {
    int i;
    delete[] prims;
    delete[] node_used;
    delete[] prim_used;
    delete[] num_edges;
    for (i = 0; i < num_nodes; i++)
      delete[] edges[i];
    delete[] edges;
  }

  void MinimumSpanningTree::apply(EdgePrimitiveFunctor &func, int start_node) {
    bool *visited;
    int i, j;
    if (num_nodes <= 0)
      return;
    j = start_node - min_node;
    if (j < 0 || j >= num_nodes)
      j = 0;
    for (; j < num_nodes && !node_used[j]; j++) /* seek */;
    if (j >= num_nodes)
      for (j = 0; j < num_nodes && !node_used[j]; j++) /* seek */;
    VW_ASSERT(j >= 0 && j < num_nodes && node_used[j], LogicErr() << "Unable to find a used node!");
    visited = new bool[num_nodes];
    for (i = 0; i < num_nodes; i++)
      visited[i] = false;
    for (i = j; i < num_nodes; i++) {
      if (node_used[i] && !visited[i])
        apply_(func, i + min_node, visited);
    }
    for (i = 0; i < j; i++) {
      if (node_used[i] && !visited[i])
        apply_(func, i + min_node, visited);
    }
    delete[] visited;
  }

  //FIXME: ack! recursion!
  void MinimumSpanningTree::apply_(EdgePrimitiveFunctor &func, int node, bool *visited) {
    int next_node;
    int i, j;
    j = node - min_node;
    VW_ASSERT(j >= 0 && j < num_nodes && node_used[j], LogicErr() << "Node is invalid!");
    VW_ASSERT(!visited[j], LogicErr() << "Node has already been visited!");
    visited[j] = true;
    for (i = 0; i < num_edges[j]; i++) {
      next_node = (node == edges[j][i]->node1()) ? edges[j][i]->node2() : edges[j][i]->node1();
      if (!visited[next_node - min_node]) {
        func(edges[j][i], node);
        apply_(func, next_node, visited);
      }
    }
  }

}} // namespace vw::math
