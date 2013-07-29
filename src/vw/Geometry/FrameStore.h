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
#ifndef vw_geometry_FrameStore_h
#define vw_geometry_FrameStore_h

#include <vw/Geometry/FrameHandle.h>
#include <vw/Geometry/FrameTreeNode.h>
#include <vw/Core/Thread.h>


namespace vw
{
namespace geometry
{
  // forward declaration
  class FrameStore;

    /**
     * @brief Thread-safe coordinate-frame tree class.
     * The FrameStore class implements a thread-safe interface to a
     * multi-tree. (Tree with multiple root nodes). It pretty much
     * mirrors the FrameTreeNode interface. The FrameStore adds locking
     * primitives for thread-safety to the feature-set..
     *
     * Frame tree nodes in the FrameStore are referenced by the
     * FrameHandle data structure. This indirection is done to avoid
     * accidential use of FrameTreeNode nodes from the FrameStore
     * outside the class interface.
     *
     * The frame trees are managed as heap-allocated structures. Error
     * checking is done to ensure the integrity of the tree structure,
     * as well as to prevent segementation faults on null-handles.  So
     * for instance, when re-parenting a frame, we check that the frame
     * is actually part of the FrameStore instance and when adding
     * frames we check for uniqueness of the name within the set of
     * siblings.  Error checking is omitted for methods, where passing
     * erronous arguments does only affect the the return value, but not
     * the integrity of the tree. For instance, is_root() does not check
     * for membership of the frame to the FrameStore.

   */
  class FrameStore {

  public:
    //! @{ Public data types.

      /**
       * The affine transform type used to denote a location
       */
      typedef vw::ATrans3 Transform;

    //! A vector of frame handles.
    typedef std::vector<FrameHandle> FrameHandleVector;
    //! A vector of FrameTreeNode's.
    /**
     * Those are actual nodes, that can be manipulated using the
     * TreeNode<Frame> interface.
     */
    typedef std::vector<FrameTreeNode> FrameTree;

      /**
       * @brief A vector of Transforms's.
       */
      typedef std::vector<Transform> TransformVector;

      //! Handle to the frame to update.
//      FrameHandle handle;
      //! Axis to update
//      ParameterType axis;
      //! New axis value.
//      double value;
//    };
    //! A vector of frame updates used in update_frames().
//    typedef std::vector<FrameUpdate> UpdateVector;

    //! @}

    //! Destructor
    /** Deletes all frames owned by this FrameStore instance.
     */
    ~FrameStore() throw();

    //! Get a copy of the frame tree.
    /**
     *
     * The vector holds a set of FrameTreeNode objects, which describe
     * all tree(s) of the FrameStore.  The FrameStore is copied in
     * pre-order, so the first element of the vector is the root node
     * of the first tree.
     *
     * Note that the FrameTreeNodes have a completely different
     * interface than the FrameStore. FrameTreeNodes don't have
     * locking, so this is a static snapshot of the tree.
     */
    void clone_tree(FrameTree & tree, FrameHandle rootFrame = NULL) const;

    //! Get a copy of the frame tree.
    /**
     * The rootFrame is required to be non-NULL, otherwise
     * vw::LogicErr is thrown.
     *
     * Note that the FrameTreeNodes have a completely different
     * interface than the FrameStore. FrameTreeNodes don't have
     * locking, so this is a static snapshot of the tree.
     */
    FrameTreeNode * clone_tree(FrameHandle rootFrame) const;


    //! Get name of frame.
    /**
     * @param frame
     * The frame argument is required to be non-NULL, otherwise
     * vw::LogicErr is thrown.
     */
    std::string const& name(FrameHandle frame) const;

    //! Get fully qualified name of frame, including path of all parent frames.
    /**
     * @param frame
     * The frame argument is required to be non-NULL, otherwise
     * vw::LogicErr is thrown.
     */
    std::string full_name(FrameHandle frame) const;

    /**
     * Get cloned extras object from Frame.
     * The caller takes ownership of the returned object.
     * If no extras are associated with the Frame, a NULL-pointer
     * is returned.
     * The frame argument is required to be non-NULL, otherwise
     * vw::LogicErr is thrown.
     */
    Frame::Extras * get_extras(FrameHandle frame) const;

    /**
     * Set  extras object for Frame.
     * The frame argument is required to be non-NULL, otherwise
     * vw::LogicErr is thrown.
     * The FrameStore takes ownership of the passed Extras object,
     * a NULL-pointer will delete the stored extras object.
     */
    void set_extras(FrameHandle frame, Frame::Extras * extras);

    /**
     * @return list of fully qualified names of all frames.
     */
    std::vector<std::string> frame_names() const;

    //! Return the parent Frame
    /**
     * Returns the NULL-handle if root frame.
     *
     * @param frame
     * The frame argument is required to be non-NULL, otherwise
     * vw::LogicErr is thrown.
     */
    FrameHandle parent(FrameHandle frame) const;

    //! Get the list of direct children of a frame.
    /**
     * @param frame
     * The frame argument is required to be non-NULL, otherwise
     * vw::LogicErr is thrown.
     */
    FrameHandleVector children(FrameHandle frame = NULL_HANDLE) const;

    //! Lookup a frame by name.
    /**
     * @param name
     *
     * Note our special lookup naming semantics: Frame names are
     * alphanumeric only. . and / are reserved characters.  A frame
     * name can be specified by giving it's absolute or relative path
     * in Unix file-system convention: /rootNode/myNode or
     * ../../grandParent/uncleFrame.
     *
     * A special wild-card is the ... which starts a breadth-first
     * expansion of the tree.  So .../myNode will return the first
     * node named myNode beneath the scope-node in a bread-first
     * expansion.  As the ordering of children is not defined, it is
     * not guaranteed which node is returned if multiple nodes with
     * the same name are specified at the same depth level.
     *
     * @param scope
     * If a non-NULL scope frame is passed as second parameter, the
     * search is restricted to the sub-tree spawned by this frame.
     */
    FrameHandle lookup(std::string const& name, FrameHandle scope = NULL) const;

    //! Adding a frame to the frame store.
    /**
     * @param name
     * @param parent
     * @param p
     */
    FrameHandle add(std::string const& name, FrameHandle parent, Transform const& p);

      /**
       * Adding a sub-tree to the frame store
       * @param node
       * The FrameStore takes ownership of the passed sub-tree.
       * The tree must not be member of a FrameStore already. Otherwise,
       * vw::LogicErr is thrown.
       *
       * @param parent
       */
      void add(FrameTreeNode * node, FrameHandle parent);

      /**
       * Merging a tree with the the frame store
       * @param tree
       * On successful merging the FrameStore takes ownership of the passed sub-tree.
       * If vw::LogicErr is thrown, the passed tree is still owned by the caller, to
       * allow forensics.
       *  * The tree must not be member of a FrameStore already. Otherwise,
       *  vw::LogicErr is thrown.
       *  * If a NULL_HANDLE is passed as start_frame the tree is merged into the forrest.
       *  * If start_frame is not the NULL_HANDLE, the names of the start_frame and the
       *  name of the tree root_node have to match.
       *
       * @param start_frame
       * The start-node for the merge operation. The startFrame is required to have the same
       * name as the node.
       * @return True, if the tree was merged,
       * false if the tree was added as a new tree to the forest.
       */
      bool merge_tree(FrameTreeNode * tree, FrameHandle start_frame = NULL_HANDLE);

       /**
       * Delete frame from tree.
       * @param frame
       *
       * @param recursive If recursive is set to false, all children of
       * the frame will be added as root-frames to the FrameStore.
       */
      void del(FrameHandle frame, bool recursive = true);

   //! Reparent a frame.
    /**
     * @param frame
     * @param parent
     */
    void set_parent(FrameHandle frame, FrameHandle parent);

    //! Return root node of specified frame.
    /**
     * The frame-store can hold multiple-root nodes.
     * @param frame
     */
    FrameHandle root(FrameHandle frame) const;

    //! @{ Public predicates.

    //! Test if frame is a base frame.
    /**
     * That is, does not have a parent.
     * @param frame
     */
    bool is_root(FrameHandle frame) const;

      /**
       * Test if frame is a leaf frame.
       * That is, does not have any children.
       * @param frame
       */
      bool is_leaf(FrameHandle frame) const;

    //! Test if frame is somewhere up in the chain of parents of pop.
    /**
     * @param frame
     * @param pop
     */
    bool is_ancestor_of(FrameHandle frame, FrameHandle pop) const;

    //! Test if the frame belongs to this FrameStore instance.
    bool is_member(FrameHandle frame) const throw();

    //! @}

      /**
       * Return the transform of source expressed relative to frame.
       * @param frame
       * @param source
       */
      Transform get_transform(FrameHandle frame, FrameHandle source);

      /**
       * Return the transform loc, which is expressed relative to source
       * with respect to frame.
       * @param frame
       * @param source
       * @param loc
       */
      Transform get_transform_of(FrameHandle frame, FrameHandle source, Transform const& loc);

      /**
       * Return the position loc, which is expressed relative to source
       * with respect to frame.
       * @param frame
       * @param source
       * @param loc
       */
      Vector3 get_position_of(FrameHandle frame, FrameHandle source, Vector3 const& loc);

      /**
       * Set the transform of frame to update, which is expressed relative to wrt_frame.
       * @param frame
       * @param wrt_frame
       * @param update
       */
      void set_transform(FrameHandle frame, FrameHandle wrt_frame, Transform const& update);

      /**
       * Update the transform of frame to loc, expressed relative to current transform.
       * @param frame
       * @param loc
       */
      void set_transform_rel(FrameHandle frame, Transform const& loc);

      void get_frame_transforms(FrameHandleVector const& handles, TransformVector& transforms) const;
      void set_frame_transforms(FrameHandleVector const& handles, TransformVector const& transforms);

      /**
       * A static instance to a NULL handle, for reference.
       */
      static FrameHandle const NULL_HANDLE;

    protected:
      void assert_unique(std::string const& name, FrameTreeNode * parent) const;
      /** Test if the frame belongs to this FrameStore instance. */
      bool is_member(FrameTreeNode * node) const throw();

      /** Vector of FrameTreeNode pointers. */
      typedef std::vector<FrameTreeNode *> FrameTreeNodeVector;

      /** Mutex to ensure exclusive access to framestore operations. */
      mutable RecursiveMutex m_mutex;
      /** The vector of root nodes. */
      FrameTreeNodeVector m_root_nodes;
    };
  }
}


#endif // vw_geometry_FrameStore_h
