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

#include <vw/Math/RandomSet.h>
#include <vw/Core/Exception.h>

#include <algorithm>
#include <random>

namespace vw{
namespace math{

  // Given a vector of integers, pick a random subset. 
  
  // Note: The complexity of this is O(input_len * log(input_len)), which
  // may be not be good enough for some applications.
  
  void pick_random_subset(int output_len, std::vector<int> const& input, 
                          std::vector<int> & output) {

    // output_len must be non-negative
    if (output_len < 0)
      vw::vw_throw(vw::ArgumentErr() 
        << "pick_random_subset: The specified output length"
        << " is negative.");
      
    // Shuffle the input
    std::vector<int> shuffled = input;
    std::mt19937 g; // Each time this is run same random numbers should be produced
    std::shuffle(shuffled.begin(), shuffled.end(), g);

    int input_len = shuffled.size();
    if (output_len >= input_len) {
      output = shuffled;
      return;
    }

    // Wipe the output
    output.clear();

    // Copy the first several randomly shuffled elements to the output
    output.resize(output_len);
    for (int it = 0; it < output_len; it++) 
      output[it] = shuffled[it];
  }

  // Pick unsorted random indices in [0, ..., input_len - 1]
  void pick_random_indices_in_range(int input_len, int output_len,
                                    std::vector<int> & output) {
   
    // input and output len must be non-negative
    if (input_len < 0 || output_len < 0)
      vw::vw_throw(vw::ArgumentErr() 
        << "pick_random_indices_in_range: The specified input or output length"
        << " is negative.");

    std::vector<int> input(input_len);
    for (int it = 0; it < input_len; it++)
      input[it] = it;
    
    pick_random_subset(output_len, input, output);
  }
  
}} // End namespace vw::math

