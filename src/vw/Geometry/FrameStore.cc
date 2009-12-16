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
  FrameStore::name(FrameHandle frame) const throw()
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
  FrameStore::lookup(std::string const& name, FrameHandle scope)
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


  namespace 
  {
  class Append
  {
    vector<string>& s_;
  public:
    Append(vector<string>& s) :
      s_(s) 
    {}
    void operator() (FrameTreeNode const * node) 
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
  FrameStore::parent(FrameHandle frame) const throw()
  {
    RecursiveMutex::Lock lock(m_mutex);

    VW_ASSERT (frame.node != NULL,
	       vw::LogicErr("NULL handle not allowed as parameter."));

    return frame.node->parent();
  }

  /** Get the list of direct children of a frame.
   * 
   * @param frame
   */
  FrameStore::FrameHandleVector
  FrameStore::children(FrameHandle frame) const throw()
  {
    RecursiveMutex::Lock lock(m_mutex);

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
    else {
      FrameTreeNodeVector children = parent.node->children();
      FrameTreeNodeVector::const_iterator first, last = children.end();
      for (first = children.begin(); first != last; ++ first) {
	VW_ASSERT((*first)->data().name() != name,
		  vw::LogicErr("Name not unique as root node."));
      }
    }
    
    FrameTreeNode * node = new FrameTreeNode(parent.node, Frame(name, p));
    if (parent.node == NULL)
      m_root_nodes.push_back(node);

    return node;
  }
  
  void
  FrameStore::del(FrameHandle frame, bool recursive)
  {
    RecursiveMutex::Lock lock(m_mutex);

    VW_ASSERT (frame.node != NULL,
	       vw::LogicErr("NULL handle not allowed as parameter."));
    
    // erase the node from root node vector, if root node
    FrameTreeNodeVector::iterator elem =
      std::find(m_root_nodes.begin(), m_root_nodes.end(), frame.node);
    if (elem != m_root_nodes.end())
      m_root_nodes.erase(elem);

    if (frame.node != NULL) {
      if (recursive) {
	// delete subtree
	frame.node->recursive_delete();
	// decrease global node count
      }
      else {
	// safe children
	FrameTreeNodeVector children = frame.node->children();
	// insert child nodes into root-node vector
	m_root_nodes.insert(m_root_nodes.end(), children.begin(), children.end());
	// delete frame
	delete frame.node;
	// we deleted only one node
      }
    }
  }
  
  void
  FrameStore::set_parent(FrameHandle frame, FrameHandle parent) 
  {
    RecursiveMutex::Lock lock(m_mutex);

    VW_ASSERT (frame.node != NULL,
	       vw::LogicErr("NULL handle not allowed as frame parameter."));
    VW_ASSERT (parent.node == NULL || is_member(parent),
	       vw::LogicErr("None member node not allowed as parent."));

    frame.node->set_parent(parent.node);
  }
  
  bool
  FrameStore::is_root(FrameHandle frame)
  {
    RecursiveMutex::Lock lock(m_mutex);

    VW_ASSERT (frame.node != NULL,
	       vw::LogicErr("NULL handle not allowed as parameter."));
    
    return frame.node->parent() == NULL;
  }
  
  bool
  FrameStore::is_ancestor_of(FrameHandle frame, FrameHandle pop)
  {
    RecursiveMutex::Lock lock(m_mutex);

    VW_ASSERT (frame.node != NULL && pop.node != NULL,
	       vw::LogicErr("NULL handle not allowed as parameter."));

    return frame.node->is_ancestor_of(pop.node);
  }

  Frame::Location
  FrameStore::get_location_of(FrameHandle frame, FrameHandle source, Location const& loc)
  {
    RecursiveMutex::Lock lock(m_mutex);

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

    VW_ASSERT (frame.node != NULL,
	       vw::LogicErr("NULL handle not allowed as parameter."));

    vw::geometry::set_location(frame.node, wrt_frame.node, update);
  }

  bool
  FrameStore::is_member(FrameHandle frame) const throw() {
    if (frame.node != NULL) {
      FrameTreeNodeVector::const_iterator first, last = m_root_nodes.end();
      for (first = m_root_nodes.begin(); first != last; ++first) {
	if ((*first) == frame.node ||
	    (*first)->is_ancestor_of(frame.node))
	  return true;
      }
    }
    return false;
  }

  void
  FrameStore::update_frames(UpdateVector const& updates) 
  {
    UpdateVector::const_iterator first, last = updates.end();
    for (first = updates.begin(); first != last; ++first) {
      // update !!!
    }
  }

}
}
