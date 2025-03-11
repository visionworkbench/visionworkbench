
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef VW_OPENMVG_MATCHING_IND_MATCH_HPP
#define VW_OPENMVG_MATCHING_IND_MATCH_HPP

#include <vw/BundleAdjustment/openMVG_types.h>

#include <iostream>
#include <map>
#include <set>
#include <vector>

namespace VwOpenMVG {
namespace matching {

/// Structure in order to save pairwise indexed references.
/// A sort operator exist in order to remove duplicates of IndMatch series.
struct IndMatch
{
  IndMatch(IndexT i = 0, IndexT j = 0) : m_left(i), m_right(j)  {}

  friend bool operator==(const IndMatch& m1, const IndMatch& m2)  {
    return (m1.m_left == m2.m_left && m1.m_right == m2.m_right);
  }

  friend bool operator!=(const IndMatch& m1, const IndMatch& m2)  {
    return !(m1 == m2);
  }

  // Lexicographical ordering of matches. Used to remove duplicates
  friend bool operator<(const IndMatch& m1, const IndMatch& m2)  {
    return (m1.m_left < m2.m_left || (m1.m_left == m2.m_left && m1.m_right < m2.m_right));
  }

  /// Remove duplicates ((m_left, m_right) that appears multiple times)
  static bool getDeduplicated(std::vector<IndMatch> & vec_match)  {

    const size_t sizeBefore = vec_match.size();
    const std::set<IndMatch> set_deduplicated( vec_match.begin(), vec_match.end());
    vec_match.assign(set_deduplicated.begin(), set_deduplicated.end());
    return sizeBefore != vec_match.size();
  }

  IndexT m_left, m_right;  // Left, right index
};

inline std::ostream& operator<<(std::ostream & out, const IndMatch & obj) {
  return out << obj.m_left << " " << obj.m_right;
}

inline std::istream& operator>>(std::istream & in, IndMatch & obj) {
  return in >> obj.m_left >> obj.m_right;
}

using IndMatches = std::vector<matching::IndMatch>;

/// Pairwise matches (indexed matches for a pair <I,J>)
/// The interface used to store corresponding point indexes per images pairs
class PairWiseMatchesContainer {
public:
  virtual ~PairWiseMatchesContainer() {}
  virtual void insert(std::pair<Pair, IndMatches>&& pairWiseMatches) = 0;
};

//--
/// Pairwise matches (indexed matches for a pair <I,J>)
/// A structure used to store corresponding point indexes per images pairs
struct PairWiseMatches :
  public PairWiseMatchesContainer,
  public std::map<Pair, IndMatches> {
  void insert(std::pair<Pair, IndMatches> && pairWiseMatches)override {
    std::map< Pair, IndMatches >::insert(
      std::forward<std::pair<Pair, IndMatches>>(pairWiseMatches));
  }

};

inline Pair_Set getPairs(const PairWiseMatches & matches)
{
  Pair_Set pairs;
  for ( const auto & cur_pair : matches )
    pairs.insert(cur_pair.first);
  return pairs;
}

}  // namespace matching
}  // namespace VwOpenMVG

#endif // VW_OPENMVG_MATCHING_IND_MATCH_HPP
