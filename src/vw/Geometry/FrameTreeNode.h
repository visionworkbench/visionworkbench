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
#ifndef vw_geometry_FrameTreeNode_h
#define vw_geometry_FrameTreeNode_h

#include <vw/Geometry/TreeNode.h>
#include <vw/Geometry/Frame.h>

namespace vw
{
  namespace geometry {
    /**
     * @brief A tree node with a Frame as payload.
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

    /**
     * @brief Get transform of the source frame relative to the traget frane.
     * @param target The coordinate frame to convert to. If NULL is passed as source frame, source->parent() is assumed
     * and source->data().transform() is returned.
     * @param source The source coordinate frame. If NULL, the identity-matrix is returned.
     */
    Frame::Transform get_transform(FrameTreeNode const * target, FrameTreeNode const * source);
    /**
     * Get transform of transform, expressed relative to the source frame, relative to the target frame.
     */
    Frame::Transform get_transform_of(FrameTreeNode const * target, FrameTreeNode const * source,
        Frame::Transform const& transform);
    /**
     * Set transform of the frame to the transform specified realtive to the source frame.
     */
    void set_transform(FrameTreeNode * frame, FrameTreeNode const * source,
                                       Frame::Transform const& transform);

    /**
     * @brief Lookup frame by name
     * See the FrameTreeNode typedef for the requirements on Frame naming.
     * @paraam start_frame The node relative to which the search starts.
     * @param path An expression describing the path to the searched frame.
     *  * '/' is the frame-name delimiter. '.' is used for directory wildcards.
     *  Frame names can not contain these characters.
     *  * A leading "/" starts the lookup from the start_frame.root().
     *  * ".", "..", and "..." are treated as following
     *    * ".": this frame - being ignored in practice
     *    * "..": parent frame
     *    * "...": 0 to n frames down in breadth-first order
     */
    FrameTreeNode * lookup(FrameTreeNode * start_frame, std::string const& path);

    inline
    Frame::Transform
    get_transform_of(FrameTreeNode const * target, FrameTreeNode const * source,
                     Frame::Transform const& trans)
    {
      return get_transform(target, source) * trans;
    }

    inline
    void
    set_transform(FrameTreeNode * frame, FrameTreeNode const * source,
                  Frame::Transform const& trans)
    {
      if (frame != NULL) {
        if (source == NULL) {
          frame->data().set_transform(trans);
        }
        else {
          frame->data().set_transform(get_transform_of(frame->parent(), source, trans));
        }
      }
    }

    /**
     * @brief Merging source_tree into target tree.
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
