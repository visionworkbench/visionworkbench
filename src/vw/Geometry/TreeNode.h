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


// -*- C++ -*-
#ifndef vw_geometry_TreeNode_h
#define vw_geometry_TreeNode_h

#include <algorithm>
#include <deque>
#include <vector>
#include <iterator>

#include <iostream>

namespace vw
{
  namespace geometry {
    /**
     * Generic Tree container.
     *
     * The class defines a tree node, which aggregates a templetaized
     * data type as payload.  Each tree node holds a single-linked list
     * of children and a pointer to the parent node.  A root nodes
     * parent pointer is NULL.
     *
     * Note that this class is not thread-safe. As an example of how to
     * use TreeNode<> in a multi-threaded setup, see the FrameStore
     * class.
     *
     * TreeNodes can be heap, stack or statically allocated. In
     * conequence, TreeNodes don't ownership of their
     * children. Destruction of a TreeNode is therefore not recursive by
     * default. For heap-allocated trees, an explicit recursive_delete
     * method is available.
     */
    template<class DATA>
    class TreeNode
    {
    public:
      /** Vector of tree node pointers. */
      typedef std::vector<TreeNode *> NodeVector;

    protected:
      DATA m_data; /**< Payload from template */
      TreeNode * m_parent; /**< Pointer towards root of the tree. Root node m_parent holds NULL */
      TreeNode * m_first_child; /**< Head of the sibling list. NULL for leave nodes. */
      /**
       * @brief Single-linked list of siblings.
       * The list is NULL terminated. List-head is m_first_child of the
       * parent node.
       */
      TreeNode * m_next_sibling;

    private:
      /**
       * Per class data for the @ancestry() method.
       *
       * It represents a cached vector for ancestry.  This helps to
       * minimize memory allocations.  But makes the ancestry() method
       * not re-entrant.
       */
      mutable NodeVector ancestry_;

    public:
      /**
       * Default constructor
       */
      TreeNode() :
          m_parent(NULL),
          m_first_child(NULL),
          m_next_sibling(NULL) {
        ancestry_.reserve(16);
      }
      /**
       * Initializing constructor.
       */
      TreeNode(TreeNode * parent, DATA const& data) :
          m_data(data),
          m_parent(parent),
          m_first_child(NULL),
          m_next_sibling(NULL) {
        ancestry_.reserve(16);

        if (m_parent != NULL)
          m_parent->register_child(this);
      }

      /**
       * @brief Copy constructor, to preserving the original tree.
       * In order not to mess with the tree structure, only the data payload is copied.
       * So this creates a new root node, with no children.
       */
      TreeNode(TreeNode const& node) :
          m_data(node.m_data),
          m_parent(NULL),
          m_first_child(NULL),
          m_next_sibling(NULL) {
        ancestry_.reserve(16);
      }

      /**
       * Destructor.
       * The destructor only deletes the node itself. All children are
       * marked as root nodes. For recursive deletion of the tree, see
       * recursive_delete().
       */
      ~TreeNode() throw() {
        if (m_parent != NULL) {
          m_parent->deregister_child(this);
        }
        while (m_first_child != NULL) {
          m_first_child->set_parent(NULL); // this alters m_first_child !
        }
      }

      /**
       * @brief Assignement operator, preserving the original tree.
       * In order not to mess with the tree structure, only the data payload is copied.
       * So this creates a new root node, with no children.
       */
      TreeNode& operator= (TreeNode const& rhs) {
        if (&rhs != this) {
          m_data = rhs.m_data;
        }
        return *this;
      }

      /**
       * Recursive destructor of tree.
       * NOTE requires heap-allocated tree!
       */
      void recursive_delete() throw() {
        while (m_first_child) {
          m_first_child->recursive_delete();
        }
        delete this;
      }

      /**
       * @brief Clone node and its offsprings.
       * The new node is becoming a child of param parent.
       * Default is NULL, which creates a new root node.
       */
      TreeNode * clone(TreeNode * parent = NULL) const {
        TreeNode * node = new TreeNode(parent, this->data());
        TreeNode * child = this->m_first_child;
        while (child != NULL) {
          child->clone(node);
          child = child->m_next_sibling;
        }
        return node;
      }

      /**
       * @brief Clone node and its offsprings.
       * The new node is becoming a child of param parent.
       * Default is NULL, which creates a new root node.
       * The front() element of the vector is the root node of the tree.
       *
       * In this version, all nodes are appende to the provided provided
       * param nodes argument. You can not call delete or
       * recursive_delete() on any node of this vector.  Copying and
       * appending to the vector is allowed, as long as the tree stays
       * in weak order.
       *
       * \note This only works, if:
       * assert(nodes.capacity() - nodes.size() > num_offsprings())
       * also, you can not assign/copy-construct the returned vector
       * or increase it's size beyond it's capacity.
       */
      void clone_vec(TreeNode * parent, std::vector<TreeNode>& nodes) const {

        // protect against tree size overrun
        if (nodes.capacity() <= nodes.size()) {
          nodes.clear();
          return;
        }

        nodes.push_back(TreeNode(parent, this->data()));

        TreeNode * node = &nodes.back();
        if (parent) {
          node->set_parent(parent);
        }

        TreeNode * child = this->m_first_child;
        while (child != NULL) {
          child->clone_vec(node, nodes);
          child = child->m_next_sibling;
        }
      }

      /**
       * @brief Calculate number of nodes in the tree.
       * Tree size is not cached anywhere so this method has linear run-time
       * in the number of tree nodes.
       */
      unsigned int num_offsprings() const throw() {
        unsigned int offsprings = 0;
        TreeNode * child = this->m_first_child;
        while (child != NULL) {
          offsprings += 1 + child->num_offsprings();
          child = child->m_next_sibling;
        }
        return offsprings;
      }

      /**
       * Const accessor method to the data.
       */
      DATA const& data() const throw() {
        return m_data;
      }
      /**
       * Accessor method to the data.
       */
      DATA & data() throw() {
        return m_data;
      }

      /**
       * @brief Set the parent of the node.
       * This actually re-parents the sub-tree. So if the node had a parent before,
       * it gets removed from the parent's list ob childs, too.
       *
       * Setting the parent to NULL makes the node a root-node.
       */
      void set_parent(TreeNode * parent) throw() {
        if (m_parent != NULL)
          m_parent->deregister_child(this);
        m_parent = parent;
        if (m_parent != NULL)
          m_parent->register_child(this);
      }

      /**
       * @brief Get the parent of the node.
       * If the node is a root-node, NULL is returned.
       */
      TreeNode * parent() const throw() {
        return m_parent;
      }

      /**
       * Predicate: returns true if node is root-node
       */
      bool is_root() const throw() {
        return m_parent == NULL;
      };

      /**
       * Predicate: returns true if node is a leaf-node
       */
      bool is_leaf() const throw() {
        return m_first_child == NULL;
      };

      /**
       * Return the root node of the tree.
       * This method is linear in the depth of the tree.
       */
      TreeNode * root() const throw() {
        TreeNode * node = const_cast<TreeNode *>(this);
        while (node->m_parent)
          node = node->m_parent;
        return node;
      }

      /**
       * Predicate: returns true if this node is ancestor of the param node.
       */
      bool is_ancestor_of(TreeNode * node) const throw() {
        while (node) {
          node = node->m_parent;
          if (node == this)
            return true;
        }

        return false;
      }

      /**
       * @brief Calculate the ancestry.
       * This method returns a vector containg pointers to nodes
       * from the root node down to and including this node.
       *
       * This method is linear in run-time with the depth of the tree.
       * If param dontUpdate is given, the ancestry is not re-calculated,
       * but a reference to the previously calculated ancestry is returned.
       */
      NodeVector const& ancestry(bool dontUpdate = false) const {
        if (!dontUpdate) {
          TreeNode * iter = const_cast<TreeNode *>(this);
          ancestry_.clear();
          while (iter) {
            ancestry_.push_back(iter);
            iter = iter->m_parent;
          }
          std::reverse(ancestry_.begin(), ancestry_.end());
        }
        return ancestry_;
      }

      /**
       * @brief Calculate the last common ancestor.
       * Returns -1 if the nodes are not connected.  If either node is
       * a direct ancestor of the other, that index in ancestry() is returned.
       *
       * The algorithm is of complexity O(n * log(n)) with respect to the
       * depth of the tree.
       */
      int last_common_ancestor_index(TreeNode const * node) const {
        NodeVector const& a1 = node->ancestry(); // min length is 1
        NodeVector const& a2 = this->ancestry(); // min length is 1

        // unique mismatch algorithm with O(log(n))
        // it's pretty much a lower_bound algorithm,
        // only that we don't compare a single element
        // operator < in this case is defined as ==
        // and we return lower_bound() - 1

        typedef typename NodeVector::const_iterator::difference_type DistanceType;

        DistanceType len = std::min(a1.size(), a2.size());
        DistanceType base = 0;
        DistanceType index;
        DistanceType half;

        while (len > 0) {
          half = len >> 1;
          index = base;
          index += half;

          if (a1[index] == a2[index]) {
            base = index;
            ++base;
            len = len - half - 1;
          }
          else {
            len = half;
          }
        }

        return ((int)base) - 1;
      }
      /**
       * @brief Calculate the last common ancestor.
       * Returns NULL if the nodes are not connected.  If either node is
       * a direct ancestor of the other, that node is returned.
       *
       * The algorithm is of complexity O(n * log(n)) with respect to the
       * depth of the tree.
       */
      TreeNode * last_common_ancestor(TreeNode const * node) const {
        int index = last_common_ancestor_index(node);

        return (index >= 0)? ancestry_[index] : NULL;
      }

      /**
       * Return a vector of pointers to the direct childs of the node.
       */
      NodeVector children() const {
        NodeVector c;

        TreeNode * node = this->m_first_child;
        while (node) {
          c.push_back(node);
          node = node->m_next_sibling;
        }

        return c;
      }

      /**
       * Return a vector of pointers to the direct childs of the node.
       */
      template<typename OutputIterator>
      void copy_children(OutputIterator& target) const {
        TreeNode * node = this->m_first_child;
        while (node) {
          (*target) = node;
          ++target;
          node = node->m_next_sibling;
        }
      }

      /**
       * @brief Tree travesal template.
       * Visits all nodes in pre-order.
       * The method accepts 3 functors as arguments.
       * For each node func(node) is called.
       * Before each descend down(node) is called.
       * After each ascend up(node) is called.
       */
      template<typename F, typename D, typename U>
      void
      pre_order_traverse(F& func, D& down, U& up) const {
        func(this);
        if (m_first_child) {
          down(this);
          m_first_child->pre_order_traverse(func, down, up);
          up(this);
        }
        if (m_next_sibling) {
          m_next_sibling->pre_order_traverse(func, down, up);
        }
      }


      /**
       * @brief Tree travresal template
       * Visits all nodes in post-order.
       * The method accepts 3 functors as arguments.
       * Before each descend down(node) is called.
       * After each ascend up(node) is called.
       * For each node func(node) is called.
       */
      template<typename F, typename D, typename U>
      void
      post_order_traverse(F& func, D& down, U& up) const {
        if (m_first_child) {
          down(this);
          m_first_child->post_order_traverse(func, down, up);
          up(this);
        }
        if (m_next_sibling) {
          m_next_sibling->post_order_traverse(func, down, up);
        }
        func(this);
      }

      /**
       * @brief Tree traversal template.
       * Visits all nodes in breadth-first-order.
       * The method accepts 1 functors as arguments.
       * For each node of the sub-tree, including the called object,
       * func(node) is called.
       */
      template<typename F>
      void
      breadth_first_traverse(F& func) const {
        std::deque<TreeNode const *> nodes;
        nodes.push_back(this);
        while (!nodes.empty()) {
          func(nodes.front());
          std::back_insert_iterator<std::deque<TreeNode const *> > iter(nodes);
          nodes.front()->copy_children(iter);
          nodes.pop_front();
        }
      }

    protected:
      /** Internal helper method: registers child node at parent. */
      void register_child(TreeNode * child) throw() {
        child->m_next_sibling = this->m_first_child;
        this->m_first_child = child;
      }
      /** Internal helper method: deregisters child node at parent. */
      void deregister_child(TreeNode * child) throw() {
        TreeNode * node = this->m_first_child;
        TreeNode ** prev = &this->m_first_child;
        while (node) {
          if (node == child) {
            *prev = child->m_next_sibling;
            child->m_next_sibling = NULL;
            break;
          }
          prev = &node->m_next_sibling;
          node = node->m_next_sibling;
        }
      }
    };
  }
}
#endif // vw_geometry_TreeNode_h
