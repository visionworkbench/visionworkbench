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
