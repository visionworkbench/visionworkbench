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


/// \file RandomSet.h
///
/// Utilities involving random numbers.
///

#ifndef __VW_MATH_RANDOM_SET_H__
#define __VW_MATH_RANDOM_SET_H__

#include <vector>

namespace vw {
namespace math {

  // Given a vector of integers, pick a random subset. Must have
  // output_len <= input_len. The output elements are not sorted.
  // Note: The complexity of this is O(input_len * log(input_len)), which
  // may be not be good enough for some applications.
  void pick_random_subset(int output_len, std::vector<int> const& input, 
                          std::vector<int> & output);

  // Pick unsorted random indices in [0, ..., input_len - 1]
  void pick_random_indices_in_range(int input_len, int output_len,
                                    std::vector<int> & output);

}} // namespace vw::math

#endif // __VW_MATH_RANDOM_SET_H__
