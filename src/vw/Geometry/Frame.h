// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// -*- C++ -*-
#ifndef vw_Math_Frame_h
#define vw_Math_Frame_h

#include "ATrans.h"

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

      /**
       * Default constructor.
       */
      Frame() {}
      /**
       * Initializing constructor.
       */
      Frame(std::string const& name, Transform const& trans = identity_matrix<4>()) :
          m_name(name),
          m_trans(trans) {}

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
      /// @}

    protected:
      std::string m_name;
      Transform m_trans;
    };
  }
}
#endif // vw_Math_Frame_h
