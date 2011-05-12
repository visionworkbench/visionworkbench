// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// -*- C++ -*-
#ifndef vw_geometry_FrameHandle_h
#define vw_geometry_FrameHandle_h

namespace vw
{
  namespace geometry
  {
    // forward declaration
    class FrameStore;
    class Frame;
    template<typename P> class TreeNode;
    typedef TreeNode<Frame> FrameTreeNode;

    /**
     * Handle to a frame tree node stored in a frame-store.
     */
    class FrameHandle
    {
      FrameTreeNode * node;
    public:
      FrameHandle() : node(0) {}
      FrameHandle(FrameTreeNode * n) :
        node(n) {}
      bool operator==(FrameHandle const& rhs) const throw() {
        return this->node == rhs.node;
      }
      bool operator!=(FrameHandle const& rhs) const throw() {
        return this->node != rhs.node;
      }

      friend class FrameStore;
    };
  }
}

#endif // vw_geometry_FrameHandle_h
