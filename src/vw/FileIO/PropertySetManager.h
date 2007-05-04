// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__
#ifndef __VW_FILEIO_PROPERTY_SET_MANAGER_H__
#define __VW_FILEIO_PROPERTY_SET_MANAGER_H__

#include <vw/Core/Exception.h>

#include <set>
#include <map>

// Boost
#include <boost/algorithm/string.hpp>
 
namespace vw {

  /// Class PropertySetManager manages property sets.
  template<class KeyT,class PropertyT>
  class PropertySetManager
  {
  public:
    typedef std::set<PropertyT> SetT;
    
  protected:
    typedef std::map<KeyT,SetT* > MapT;
    typedef std::pair<KeyT,SetT* > ValueT;

    MapT m;

  public:
    /// Set property prop in property set key.
    void set_property(const KeyT& key, const PropertyT& prop)
    {
      SetT* s = property_set_private(key, true);
      s->insert(prop);
    }
  
    /// Return whether property prop is in property set key.
    bool property_is_set(const KeyT& key, const PropertyT& prop)
    {
      SetT* s = property_set_private(key);
      if(!s)
        return false;
      return (s->find(prop) != s->end());
    }
    
    /// Return property set key, creating it if desired and required.
    const SetT* property_set(const KeyT& key, bool create = false)
    {
      return property_set_private(key, create);
    }

    /// Constructor.
    PropertySetManager()
    {
    }

    /// Destructor.
    ~PropertySetManager()
    {
      typename MapT::iterator i;
      for(i = m.begin(); i != m.end(); i++)
      {
        delete (*i).second;
      }
    }

  protected:
    /// Return property set key, creating it if desired and required.
    SetT* property_set_private(const KeyT& key, bool create = false)
    {
      typename MapT::iterator i = m.find(key);
      if(i == m.end())
      {
        if(create)
        {
          SetT* s = new SetT();
          m.insert(std::make_pair(key, s));
          return s;
        }
        else
          return 0;
      }
      return (*i).second;
    }

  };

} // namespace vw

#endif // __VW_FILEIO_PROPERTY_SET_MANAGER_H__
