#include "FrameTreeNode.h"

#include <algorithm>

namespace vw
{
namespace geometry
{
  using namespace std;

  Frame::Location
  get_location(FrameTreeNode const * wrtFrame, FrameTreeNode const * source)
  {
    Frame::Location loc = vw::identity_matrix(4);

    if (source != NULL) {
    
      if (wrtFrame == NULL)
	return source->data().location();
      
      if (wrtFrame != source) {
	
	FrameTreeNode const * ancestor = wrtFrame->last_common_ancestor(source);
	
	if (ancestor != NULL) {
	  // from the origin frame to the ancestor
	  {
	    FrameTreeNode::NodeVector const& nodes = wrtFrame->ancestry(true);
	    FrameTreeNode::NodeVector::const_iterator iter =
	      std::find(nodes.begin(), nodes.end(), ancestor);
	    
	    assert (iter != nodes.end());
	    
	    ++iter;
	    if (iter != nodes.end()) {
	      for (; iter != nodes.end(); ++iter) {
		loc *= (*iter)->data().location();
	      }
	      loc = inverse(loc);
	    }
	  }
	  
	  
	  // from the last common ancestor to the wrt frame
	  {	
	    FrameTreeNode::NodeVector const& nodes = source->ancestry(true);
	    FrameTreeNode::NodeVector::const_iterator iter =
	      std::find(nodes.begin(), nodes.end(), ancestor);
	    
	    assert (iter != nodes.end());
	    
	    for (; iter != nodes.end(); ++iter) {
	      loc *= (*iter)->data().location();
	    }
	  }
	}
      }
   }
   
   return loc;
  }

  namespace 
  {
    struct FtnLess
    {
      bool operator () (FrameTreeNode const * lhs, FrameTreeNode const * rhs)
      {
	return lhs->data().name() < rhs->data().name();
      }
    };
  }

  void 
  merge_frame_trees(FrameTreeNode * target_tree, FrameTreeNode * source_tree)
  {
    if (target_tree->data().name() != source_tree->data().name())
      return;

    FrameTreeNode::NodeVector src_children = target_tree->children();

    if (src_children.size() > 0) {
      FrameTreeNode::NodeVector tgt_children = source_tree->children();
      
      FtnLess less;
      sort(src_children.begin(), src_children.end(), less);
      sort(tgt_children.begin(), tgt_children.end(), less);

      FrameTreeNode::NodeVector::const_iterator first, last = src_children.end();
      FrameTreeNode::NodeVector::const_iterator tgt_iter = tgt_children.end();
      for (first = src_children.begin(); first != last; ++first) {
	while (tgt_iter != tgt_children.end() &&
	       (*first)->data().name() > (*tgt_iter)->data().name()) {
	  ++tgt_iter;
	}

	if (tgt_iter == tgt_children.end() ||
	    (*first)->data().name() < (*tgt_iter)->data().name()) {
	  (*first)->set_parent(target_tree);
	}
	else if ((*first)->data().name() == (*tgt_iter)->data().name()) {
	  merge_frame_trees(*tgt_iter, *first);
	}
      }

    }
  }
}
}
