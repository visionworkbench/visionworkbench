// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef VW_OPENMVG_TYPES_HPP
#define VW_OPENMVG_TYPES_HPP

#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <vector>

/**
* @brief Main namespace of openMVG API
*/
namespace VwOpenMVG
{

/// Portable type used to store an index
using IndexT = uint32_t;

/// Standard Pair of IndexT
using Pair = std::pair<IndexT, IndexT>;

/// Set of Pair
using Pair_Set = std::set<Pair>;

} // namespace VwOpenMVG

#endif  // VW_OPENMVG_TYPES_HPP
