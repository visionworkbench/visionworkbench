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
#ifndef vw_Math_Frame_h
#define vw_Math_Frame_h

#include <vw/Geometry/ATrans.h>

#include <string>

namespace vw
{
  namespace geometry {

    /**
     * Class representing a named coordinate transform.
     *
     * This class is mostly used in combination with the TreeNode
     * class to create a frame tree.
     */
    class Frame
    {
    public:
      typedef vw::ATrans3 Transform;

      class Extras
      {
      public:
        virtual ~Extras() throw() {}
        virtual Extras * clone() = 0;
      };

      /**
       * Default constructor.
       */
      Frame() : m_extras(NULL) {}
      /**
       * Initializing constructor.
       */
      Frame(std::string const& name, Transform const& trans = identity_matrix<4>()) :
          m_name(name),
          m_trans(trans),
          m_extras(NULL)
      {}
      Frame(Frame const& rhs) :
        m_name(rhs.m_name),
        m_trans(rhs.m_trans),
        m_extras((rhs.m_extras == NULL)? NULL : rhs.m_extras->clone())
      {}

      ~Frame() throw()
      {
        delete m_extras;
      }

      Frame& operator = (Frame const& rhs) {
        if (&rhs != this) {
          m_name = rhs.m_name;
          m_trans = rhs.m_trans;
          delete m_extras;
          m_extras = (rhs.m_extras == NULL)? NULL : rhs.m_extras->clone();
        }
        return *this;
      }

      /// @{ Accessor methods

      /** Access name field. */
      std::string const& name() const throw() {
        return m_name;
      }
      /** Set name field. */
      void set_name(std::string const& name) {
        m_name = name;
      }
      /** Access transform field. */
      Transform const& transform() const throw() {
        return m_trans;
      }
      /** Mutable access to transform field. */
      Transform& transform() throw() {
        return m_trans;
      }
      /** Set transform field. */
      void set_transform(Transform const& trans) {
        m_trans = trans;
      }

      Extras * extras() const throw() {
        return m_extras;
      }

      void set_extras(Extras * extras) {
        delete m_extras;
        m_extras = extras;
      }
      /// @}

    protected:
      std::string m_name;
      Transform m_trans;
      Extras * m_extras;
    };
  }
}
#endif // vw_Math_Frame_h
