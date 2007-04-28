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
#ifndef __VW_FILEIO_PROPERTY_MULTI_MAP_H__
#define __VW_FILEIO_PROPERTY_MULTI_MAP_H__

#include <vw/Core/Exception.h>

#include <set>
#include <map>

// Boost
#include <boost/algorithm/string.hpp>
 
namespace vw {

  /// Class PropertyMultiMap is a map where each key can be associated
  /// with multiple values, each having given properties. The find()
  /// operation chooses between the set of items registered with the
  /// given key based on the properties with which each item is
  /// associated, the requested properties, and the chosen scoring
  /// strategy.
  template<class KeyT,class DataT,class PropertyT>
  class PropertyMultiMap
  {
  protected:
    typedef std::list<std::pair<DataT,const std::set<PropertyT>* > > MapValueType;
    typedef std::map<KeyT,MapValueType> MapType;

    MapType property_map;

  public:
    typedef std::pair<KeyT,DataT> ValueT;
    typedef int (*ScoringStrategy)(const std::set<PropertyT>* prop, const std::list<PropertyT>* prop_list);

    /// Register an item and its set of supported properties.
    void insert(const ValueT& item, const std::set<PropertyT>* prop = 0) {
      MapValueType v;
      typename MapType::iterator i = property_map.find(item.first);
      if(i == property_map.end()) {
        v.push_back(std::make_pair(item.second, prop));
        property_map.insert(std::make_pair(item.first, v));
      }
      else {
        (*i).second.push_back(std::make_pair(item.second, prop));
      }
    }

    /// Return an item (and whether an item was found) based on a given key,
    /// a list of requested properties, and a scoring strategy.
    bool find(ValueT& item, const KeyT& key, const std::list<PropertyT>* prop_list = 0, ScoringStrategy scoring_strategy = &score_all_required) const {
      int score;
      DataT best_data;
      int best_score = 0;
      typename MapType::const_iterator i;
      typename MapValueType::const_iterator vi;
      i = property_map.find(key);
      if(i == property_map.end())
        return false;
      if(!prop_list) {
        if((*i).second.empty())
          return false;
        else {
          item = std::make_pair(key, (*i).second.front().first);
          return true;
        }
      }
      for(vi = (*i).second.begin(); vi != (*i).second.end(); vi++) {
        score = scoring_strategy((*vi).second, prop_list);
        if(score > best_score) {
          best_score = score;
          best_data = (*vi).first;
        }
      }
      if(best_score > 0) {
        item = std::make_pair(key, best_data);
        return true;
      }
      return false;
    }

    PropertyMultiMap() {}
    ~PropertyMultiMap() {}

    /// The available scoring strategies:
    /// score_all_required(): only choose an item if it supports all requested properties
    /// score_most_present(): choose the item that supports the greatest number of
    ///   requested properties
    /// score_priority(): if all items are missing support for at least one requested
    ///   property, choose the item where the first missing property is of lowest priority
    ///   (where priority is given by closeness to the beginning of the list of requested
    ///   properties)
    /// In all cases, if more than one item is equally suitable, the tie is broken by which
    ///   item was inserted first.
    
    static int score_all_required(const std::set<PropertyT>* prop, const std::list<PropertyT>* prop_list) {
      typename std::list<PropertyT>::const_iterator i;
      for(i = prop_list->begin(); i != prop_list->end(); i++) {
        if(!(prop && prop->find(*i) != prop->end()))
          return 0;
      }
      return 1;
    }
    
    static int score_most_present(const std::set<PropertyT>* prop, const std::list<PropertyT>* prop_list) {
      int score = 1;
      typename std::list<PropertyT>::const_iterator i;
      for(i = prop_list->begin(); i != prop_list->end(); i++) {
        if(prop && prop->find(*i) != prop->end())
          score++;
      }
      return score;
    }
    
    static int score_priority(const std::set<PropertyT>* prop, const std::list<PropertyT>* prop_list) {
      int score = 1;
      typename std::list<PropertyT>::const_iterator i;
      for(i = prop_list->begin(); i != prop_list->end(); i++) {
        if(prop && prop->find(*i) != prop->end())
          score++;
        else
          return score;
      }
      return score;
    }
    
  };

} // namespace vw

#endif // __VW_FILEIO_PROPERTY_MULTI_MAP_H__
