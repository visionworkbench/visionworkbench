// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include "FrameStore.h"

#include "vw/Core/Exception.h"

#include <algorithm>

namespace vw
{
namespace geometry
{
  using namespace std;

  FrameHandle const FrameStore::NULL_HANDLE = FrameHandle(NULL);

  FrameStore::~FrameStore() throw()
  {
    // delete all frames
    FrameTreeNodeVector::const_iterator first, last = m_root_nodes.end();
    for (first = m_root_nodes.begin(); first != last; ++first) {
      (*first)->recursive_delete();
    }
  }

  void
  FrameStore::clone_tree(FrameTree & tree, FrameHandle rootFrame) const
  {
    RecursiveMutex::Lock lock(m_mutex);

    tree.clear();
    if (rootFrame.node == NULL) {
      int num_nodes = 0;
      FrameTreeNodeVector::const_iterator first, last = m_root_nodes.end();
      for (first = m_root_nodes.begin(); first != last; ++first) {
        num_nodes += (*first)->num_offsprings() + 1;
      }
      tree.reserve(num_nodes);
    }
    else {
      tree.reserve(rootFrame.node->num_offsprings() + 1);
    }

    if (rootFrame.node == NULL) {
      FrameTreeNodeVector::const_iterator first, last = m_root_nodes.end();
      for (first = m_root_nodes.begin(); first != last; ++first) {
        (*first)->clone_vec(NULL, tree);
      }
    }
    else {
      rootFrame.node->clone_vec(NULL, tree);
    }
  }

  FrameTreeNode *
  FrameStore::clone_tree(FrameHandle rootFrame) const
  {
    RecursiveMutex::Lock lock(m_mutex);

    VW_ASSERT (rootFrame.node != NULL,
               vw::LogicErr("NULL handle not allowed as parameter."));

    return rootFrame.node->clone();
  }

  std::string const&
  FrameStore::name(FrameHandle frame) const
  {
    RecursiveMutex::Lock lock(m_mutex);

    VW_ASSERT (frame.node != NULL,
               vw::LogicErr("NULL handle not allowed as parameter."));

    return frame.node->data().name();
  }

  std::string
  FrameStore::full_name(FrameHandle frame) const
  {
    RecursiveMutex::Lock lock(m_mutex);

    VW_ASSERT (frame.node != NULL,
               vw::LogicErr("NULL handle not allowed as parameter."));

    string name;
    FrameTreeNodeVector const& ancestry = frame.node->ancestry();
    FrameTreeNodeVector::const_iterator first, last = ancestry.end();
    for (first = ancestry.begin(); first != last; ++first) {
      name += "/" + (*first)->data().name();
    }

    return name;
  }

  FrameHandle
  FrameStore::lookup(std::string const& name, FrameHandle scope) const
  {
    string searchName = name;

    // if not explicitly state otherwise, we search for .../name
    if (!name.empty() && name[0] != '/') {
      if (name.length() < 4 || name.substr(0, 4) != ".../") {
        searchName = ".../" + name;
      }
    }

    if (scope == NULL ||
        (!searchName.empty() && searchName[0] == '/')) {
      FrameTreeNodeVector::const_iterator first, last = m_root_nodes.end();
      for (first = m_root_nodes.begin(); first != last; ++first) {
        FrameTreeNode * node = ::vw::geometry::lookup(*first, searchName);

        if (node != NULL)
          return FrameHandle(node);
      }
      return NULL;
    }

    return ::vw::geometry::lookup(scope.node, searchName);
  }


    FrameHandle
    FrameStore::add(std::string const& name, FrameHandle parent, Transform const& p)
    {
      string full_name;
      vector<FrameTreeNode *> const& ancestry = node->ancestry();
      vector<FrameTreeNode *>::const_iterator first, last = ancestry.end();
      for (first = ancestry.begin(); first != last; ++first) {
        full_name += "/" + (*first)->data().name();
        }

      s_.push_back(full_name);
    }
  };

  class Nothing
  {
  public:
    void operator() (FrameTreeNode const *) {}
  };
  }

  std::vector<std::string>
  FrameStore::frame_names() const
  {
    RecursiveMutex::Lock lock(m_mutex);

    vector<string> names;

    Nothing n;
    Append app(names);
    FrameTreeNodeVector::const_iterator first, last = m_root_nodes.end();
    for (first = m_root_nodes.begin(); first != last; ++first) {
      (*first)->pre_order_traverse(app, n, n);
    }

    return names;
  }


  FrameHandle
  FrameStore::parent(FrameHandle frame) const
  {
    RecursiveMutex::Lock lock(m_mutex);

    VW_ASSERT (frame.node != NULL,
               vw::LogicErr("NULL handle not allowed as parameter."));

    return frame.node->parent();
  }

  FrameHandle
  FrameStore::root(FrameHandle frame) const
  {
    RecursiveMutex::Lock lock(m_mutex);

      // don't do anything if old and new parent are the same
      if (frame.node->parent() == parent.node)
        return;

    return frame.node->root();
  }

      // if root node, delete there
      FrameTreeNodeVector::iterator node =
        find(m_root_nodes.begin(), m_root_nodes.end(), frame.node);
      if (node != m_root_nodes.end()) {
        m_root_nodes.erase(node);
      }

      frame.node->set_parent(parent.node);

      if (parent.node == 0) {
        m_root_nodes.push_back(frame.node);
      }
    }

    if (frame.node == NULL) {
      return FrameHandleVector(m_root_nodes.begin(), m_root_nodes.end());
    }

    FrameTreeNodeVector c = frame.node->children();
    return FrameHandleVector(c.begin(), c.end());
  }



  FrameHandle
  FrameStore::add(std::string const& name, FrameHandle parent, Location const& p)
  {
    VW_ASSERT (name.length() != 0,
               vw::LogicErr("None empty name required for frame."));

    RecursiveMutex::Lock lock(m_mutex);

    VW_ASSERT (parent.node == NULL || is_member(parent),
               vw::LogicErr("None member node not allowed as parent."));


    // check frame name uniqueness
    if (parent.node == NULL) {
      FrameTreeNodeVector::const_iterator first, last = m_root_nodes.end();
      for (first = m_root_nodes.begin(); first != last; ++ first) {
        VW_ASSERT((*first)->data().name() != name,
                  vw::LogicErr("Name not unique as root node."));
      }
    }

    Frame::Transform
    FrameStore::get_transform_of(FrameHandle frame, FrameHandle source, Transform const& trans)
    {
      RecursiveMutex::Lock lock(m_mutex);

      VW_ASSERT (frame.node != NULL && source.node != NULL,
                 vw::LogicErr("NULL handle not allowed as parameter."));

      return vw::geometry::get_transform_of(frame.node, source.node, trans);
    }

    FrameTreeNode * node = new FrameTreeNode(parent.node, Frame(name, p));
    if (parent.node == NULL)
      m_root_nodes.push_back(node);

    Vector3
    FrameStore::get_position_of(FrameHandle frame, FrameHandle source, Vector3 const& trans)
    {
      RecursiveMutex::Lock lock(m_mutex);

      VW_ASSERT (frame.node != NULL && source.node != NULL,
                 vw::LogicErr("NULL handle not allowed as parameter."));

      return vw::geometry::get_transform(frame.node, source.node) * trans;
    }

    Frame::Transform
    FrameStore::get_transform(FrameHandle frame, FrameHandle source)
    {
      RecursiveMutex::Lock lock(m_mutex);

      VW_ASSERT (source.node != NULL,
                 vw::LogicErr("NULL handle not allowed as frame parameter."));

      return vw::geometry::get_transform(frame.node, source.node);
    }
  }

  void
  FrameStore::set_parent(FrameHandle frame, FrameHandle parent)
  {
    RecursiveMutex::Lock lock(m_mutex);

    void
    FrameStore::set_transform(FrameHandle frame, FrameHandle wrt_frame, Transform const& update)
    {
      RecursiveMutex::Lock lock(m_mutex);

    frame.node->set_parent(parent.node);
  }

  bool
  FrameStore::is_root(FrameHandle frame) const
  {
    RecursiveMutex::Lock lock(m_mutex);

      vw::geometry::set_transform(frame.node, wrt_frame.node, update);
    }

    void
    FrameStore::set_transform_rel(FrameHandle frame, Transform const& update)
    {
      RecursiveMutex::Lock lock(m_mutex);

    return frame.node->is_ancestor_of(pop.node);
  }

      vw::geometry::set_transform(frame.node, NULL, update);
    }

    VW_ASSERT (frame.node != NULL && source.node != NULL,
               vw::LogicErr("NULL handle not allowed as parameter."));

    return vw::geometry::get_location_of(frame.node, source.node, loc);
  }

  Frame::Location
  FrameStore::get_location(FrameHandle frame, FrameHandle source)
  {
    RecursiveMutex::Lock lock(m_mutex);

    FrameTreeNode * src = source.node;
    if (src == NULL)
      src = frame.node->parent();

    VW_ASSERT (frame.node != NULL,
               vw::LogicErr("NULL handle not allowed as frame parameter."));

    return vw::geometry::get_location(frame.node, src);
  }

  void
  FrameStore::set_location(FrameHandle frame, FrameHandle wrt_frame, Location const& update)
  {
    RecursiveMutex::Lock lock(m_mutex);

    bool
    FrameStore::merge_tree(FrameTreeNode * tree, FrameHandle start_frame)
    {
      VW_ASSERT (!is_member(tree),
                 vw::LogicErr("Merged tree must not yet be member of the FrameStore."));

      // if a start node is given for mergin
      if (start_frame != NULL_HANDLE) {
        VW_ASSERT (tree->data().name() == start_frame.node->data().name(),
                   vw::LogicErr("Tree root node does not match start node."));
        vw::geometry::merge_frame_trees(tree, start_frame.node);
        tree->recursive_delete();
        return true;
      }

    vw::geometry::set_location(frame.node, wrt_frame.node, update);
  }

  bool
  FrameStore::is_member(FrameHandle frame) const throw() {
    if (frame.node != NULL) {
      FrameTreeNodeVector::const_iterator first, last = m_root_nodes.end();
      for (first = m_root_nodes.begin(); first != last; ++first) {
        if ((*first)->data().name() == tree->data().name()) {
          vw::geometry::merge_frame_trees(*first, tree);
          tree->recursive_delete();
          return true;
        }
      }
      return false;
    }
    return false;
  }

    void
    FrameStore::get_frame_transforms(FrameHandleVector const& frames, TransformVector& transforms) const
    {
      transforms.clear();
      transforms.reserve(frames.size());

      RecursiveMutex::Lock lock(m_mutex);
      FrameHandleVector::const_iterator first, last = frames.end();
      for (first = frames.begin(); first != last; ++first) {
        transforms.push_back(first->node->data().transform());
      }
    }

    void
    FrameStore::set_frame_transforms(FrameHandleVector const& frames, TransformVector const& transforms)
    {
      VW_ASSERT(frames.size() == transforms.size(),
                vw::LogicErr("Parameter vectors not of same size."));

      RecursiveMutex::Lock lock(m_mutex);
      FrameHandleVector::const_iterator first, last = frames.end();
      TransformVector::const_iterator trans = transforms.begin();
      for (first = frames.begin(); first != last; ++first, ++trans) {
        first->node->data().set_transform(*trans);
      }
    }

    bool
    FrameStore::is_member(FrameTreeNode * node) const throw()
    {
      if (node != NULL) {
        FrameTreeNodeVector::const_iterator first, last = m_root_nodes.end();
        for (first = m_root_nodes.begin(); first != last; ++first) {
          if ((*first) == node ||
                      (*first)->is_ancestor_of(node))
            return true;
        }
      }
      return false;
    }

    // check name uniqueness
    void
    FrameStore::assert_unique(std::string const& name, FrameTreeNode * parent) const
    {
      char const * const childError = "Name not unique as child node.";
      char const * const rootError = "Name not unique as root node.";
      char const * err = NULL;

      FrameTreeNodeVector children;
      FrameTreeNodeVector::const_iterator first, last;

      if (parent == NULL) {
        first = m_root_nodes.begin();
        last = m_root_nodes.end();
        err = rootError;
      }
      else {
        children = parent->children();
        first = children.begin();
        last = children.end();
        err = childError;
      }

      for (; first != last; ++ first) {
        if ((*first)->data().name() == name)
          vw_throw(vw::LogicErr(err));
      }
    }

  }

}
}
