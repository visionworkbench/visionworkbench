// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// -*- C++ -*-
#ifndef vw_geometry_FrameStore_h
#define vw_geometry_FrameStore_h

#include "FrameTreeNode.h"

#include "vw/Core/Thread.h"


namespace vw
{
namespace geometry
{
  // forward declaration
  class FrameStore;

  //! Handle to a frame tree node stored in a frame-store.
  class FrameHandle
  {
    FrameTreeNode * node;
  public:
    FrameHandle() {}
    FrameHandle(FrameTreeNode * n) : 
      node(n)
    {}
    bool operator==(FrameHandle const& rhs) { return this->node == rhs.node; }
    bool operator!=(FrameHandle const& rhs) { return this->node != rhs.node; }

    friend class FrameStore;
  };

  //! Thread-safe coordinate-frame tree class.
  /**
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

    //! The affirm transform type used to denote a location
    typedef vw::ATrans3 Location;

    //! A vector of frame handles.
    typedef std::vector<FrameHandle> FrameHandleVector;
    //! A vector of FrameTreeNode's.
    /**
     * Those are actual nodes, that can be manipulated using the
     * TreeNode<Frame> interface.
     */
    typedef std::vector<FrameTreeNode> FrameTree;

    //! Struct holding a single axis frame update.
    /** 
     * As frame-updates from robot joints are a frequent operation, we
     * try to optimize this operation, by minimizing lock operations.
     */
    struct FrameUpdate
    {
      //! Update axis denominator.
      enum ParameterType { X, Y, Z, Roll, Pitch, Yaw };

      //! Handle to the frame to update.
      FrameHandle handle;
      //! Axis to update 
      ParameterType axis;
      //! New axis value.
      double value;
    };
    //! A vector of frame updates used in @update_frames().
    typedef std::vector<FrameUpdate> UpdateVector;

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
     * @vw::LogicErr is thrown.
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
     * @vw::LogicErr is thrown.
     */
    std::string const& name(FrameHandle frame) const;

    //! Get fully qualified name of frame, including path of all parent frames.
    /**
     * @param frame
     * The frame argument is required to be non-NULL, otherwise
     * @vw::LogicErr is thrown.
     */
    std::string full_name(FrameHandle frame) const;

    //! Return list of fully qualified names of all frames.
    std::vector<std::string> frame_names() const;
    
    //! Return the parent Frame
    /** 
     * Returns the NULL-handle if root frame.
     * 
     * @param frame
     * The frame argument is required to be non-NULL, otherwise
     * @vw::LogicErr is thrown.
     */
    FrameHandle parent(FrameHandle frame) const;

    //! Get the list of direct children of a frame.
    /**
     * @param frame
     * The frame argument is required to be non-NULL, otherwise
     * @vw::LogicErr is thrown.
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
     * @param sope
     * If a non-NULL scope frame is passed as second parameter, the 
     * search is restricted to the sub-tree spawned by this frame.
     */
    FrameHandle lookup(std::string const& name, FrameHandle scope = NULL) const;
    
    //! Adding a frame to the frame store.
    /** 
     * @param frame
     * @param parent
     * @param l
     */
    FrameHandle add(std::string const& name, FrameHandle parent, Location const& p);

    //! Adding a sub-tree to the frame store
    /** 
     * @param node
     * The FrameStore takes ownership of the passed sub-tree.
     * The tree must not be member of a FrameStore already. Otherwise,
     * vw::LogicErr is throwsn
     *
     * @param parent
     */
     FrameHandle add(FrameTreeNode * node, FrameHandle parent);

    //! Merging a tree with the the frame store
    /** 
     * @param node
     * The FrameStore takes ownership of the passed sub-tree.
     * The tree must not be member of a FrameStore already. Otherwise,
     * vw::LogicErr is throwsn
     *
     * @param startFrame
     * The start-node for the merge operation. The startFrame is required to have the same
     * name as the node.
     */
     FrameHandle merge_tree(FrameTreeNode * tree, FrameHandle startFrame);

    //! Delete frame from tree.
    /** 
     * @param tree
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

    //! Test if @frame is somewhere up in the chain of parents of @pop.
    /** 
     * @param frame
     * @param pop
     */
    bool is_ancestor_of(FrameHandle frame, FrameHandle pop) const;

    //! Test if the frame belongs to this FrameStore instance.    
    bool is_member(FrameHandle frame) const throw();

    //! @}

    /** Return the location of @source with respect to @frame.
     * 
     * @param frame
     * @param wrt_frame
     */
    Location get_location(FrameHandle frame, FrameHandle source = NULL);

    /** Return the position @loc, which is expressed relative to @source 
     * with respect to @frame.
     * 
     * @param frame
     * @param wrt_frame
     * @param loc
     */
    Location get_location_of(FrameHandle frame, FrameHandle source, Location const& loc);

    /** Set the location of @frame to @update, which is expressed relative to @wrt_frame.
     * 
     * @param frame
     * @param wrt_frame
     * @param update
     */
    void set_location(FrameHandle frame, FrameHandle wrt_frame, Location const& update);
    
    /** Update the location of @frame to @loc, expressed relative to current location.
     * 
     * @param frame
     * @param loc
     */
    void set_location_rel(FrameHandle frame, Location const& loc);

    //! Update a set of frames at once.
    void update_frames(UpdateVector const& updates);


    //! A static instance to a NULL handle, for reference.
    static FrameHandle const NULL_HANDLE;

  protected:
    //! Vector of FrameTreeNode pointers.
    typedef std::vector<FrameTreeNode *> FrameTreeNodeVector;

    //! Mutex to ensure exclusive access to framestore operations.
    mutable RecursiveMutex m_mutex;
    //! The vector of root nodes.
    FrameTreeNodeVector m_root_nodes;
  };
}
}

#endif // vw_geometry_FrameStore_h
