// -*- C++ -*-
#ifndef vw_Math_Frame_h
#define vw_Math_Frame_h

#include "ATrans.h"

#include <string>

namespace vw
{
namespace geometry
{

  /** Class representing a named coordinate transform.
   *
   * This class is mostly used in combination with the TreeNode 
   * class to create a frame tree.
   */
  class Frame 
  {
  public:
    typedef vw::ATrans3 Location;

    /// Default constructor.
    Frame() {}
    /// Initializing constructor.
    Frame(std::string const& name, Location const& loc = identity_matrix<4>()) :
      m_name(name),
      m_loc(loc)
    {}

    /// @{ Accessor methods

    /// Access name field.
    std::string const& name() const throw() { return m_name; }
    /// Set name field.
    void set_name(std::string const& name) { m_name = name; }
    /// Access location field.
    Location const& location() const throw() { return m_loc; }
    /// Set location field.
    void set_location(Location const& loc) { m_loc = loc; }
    /// @}

  protected:
    std::string m_name;
    Location m_loc;
  };
}
}
#endif // vw_Math_Frame_h
