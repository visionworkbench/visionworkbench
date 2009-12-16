// -*- C++ -*-
#ifndef vw_geometry_FrameTreeNode_h
#define vw_geometry_FrameTreeNode_h

#include "TreeNode.h"
#include "Frame.h"

namespace vw
{
namespace geometry
{
  //! A tree node with a Frame as payload.
  /**
   * A frame-tree is expected to adhere to the following requirements:
   *  * Frame names of children of a node are unique.
   *  * Frame names do not contain the character '/'.
   *  * Frames are not named ".", "..", nor "...".
   *  * The affinte transform has an invertable matrix.
   *
   * These requirements are not checked explicitly in the FrameTreeNode interface.
   * Only the FrameStore interface explicitily checks these requirements and
   * flags errors if frames are added to the frame-store, that do not meet
   * these assumptions.
   */
  typedef TreeNode<Frame> FrameTreeNode;

  //! Get location of the source frame relative to the wrt_frame (wrt = with respect to). 
  /**
   * @param wrt_frame The coordinate frame to convert to. If NULL is passed as source frame, source->parent() is assumed
   * and source->data().location() is returned.
   * @param source The source coordinate frame. If NULL, the identity-matrix is returned.
   */
  Frame::Location get_location(FrameTreeNode const * wrt_frame, FrameTreeNode const * source);
  //! Get location of location, expressed relative to the source frame, relative to the wrt_frame (wrt = with respect to). 
  Frame::Location get_location_of(FrameTreeNode const * wrt_frame, FrameTreeNode const * source,
				  Frame::Location const& location);
  //! Set location of the frame to the location specified realtive to the source frame.
  void set_location(FrameTreeNode * frame, FrameTreeNode const * source,
		    Frame::Location const& location);

  //! Lookup frame by name
  /**
   * See the FrameTreeNode typedef for the requirements on Frame naming.
   * @paraam start_frame The node relative to which the search starts.
   * @param path An expression describing the path to the searched frame.
   *  * '/' is the frame-name delimiter.
   *  Frame names can not contain this character.
   *  * A leading "/" starts the lookup from the start_frame.root().
   *  * ".", "..", and "..." are treated as following
   *    * ".": this frame - being ignored in practice
   *    * "..": parent frame
   *    * "...": 0 to n frames down in breadth-first order
   */
  FrameTreeNode * lookup(FrameTreeNode * start_frame, std::string const& path);

  inline
  Frame::Location 
  get_location_of(FrameTreeNode const * wrtFrame, FrameTreeNode const * source,
		  Frame::Location const& loc)
  {
    return get_location(wrtFrame, source) * loc;
  }

  inline
  void 
  set_location(FrameTreeNode * frame, FrameTreeNode const * source,
	       Frame::Location const& loc)
  {
    if (frame != NULL)
      frame->data().set_location(get_location_of(frame->parent(), source, loc));
  }

  //! Merging source_tree into target tree.
  /**
   * Source and target tree need to start with the same root node.
   * Frame tree nodes with the same name are considered a match.
   * The transformation matrix for a matched node stays in the target tree stays untouched.
   * New branches are swallowed by the target tree.
   *
   * @WARNING: Don't merge trees of nodes from different allocation categories (heap vs stack).
   */
  void merge_frame_trees(FrameTreeNode * target_tree, FrameTreeNode * source_tree);
}
}

#endif // vw_geometry_FrameTreeNode_h
