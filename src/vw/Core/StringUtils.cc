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

#include <vw/Core/StringUtils.h>

#include <algorithm>

namespace vw {

size_t string_replace(std::string &s, std::string const& find, std::string const& replace) {
  // Keep replacing the first match found until we don't find any more.
  size_t count = 0;

  size_t pos = 0;
  while (true) {
    pos = s.find(find, pos);
    if (pos == std::string::npos)
      return count;
    s.replace(pos, find.length(), replace);
    pos++; // to not get stuck in an infinite loop
    count++;
  }
  return count;
}

std::vector<double> str_to_std_vec(std::string const& str, std::string separators) {

  std::istringstream iss;
  if (separators == "") {
    iss = std::istringstream(str);
  } else {
    // Replace the separators with space before tokenizing
    std::string proc_str = str;
    for (size_t it = 0; it < separators.size(); it++) {
      std::string find;
      find += separators[it]; // this is now a string with just one character
      std::string replace = " ";
      string_replace(proc_str, find, replace);
    }
    iss = std::istringstream(proc_str);
  }
  
  std::vector<double> vec;
  double val = 0;
    while (iss >> val) 
      vec.push_back(val);
  
  return vec;
} 

// Parses a string containing a list of numbers separated by commas or spaces
// TODO(oalexan1): Replace this with str_to_std_vec. Needs testing.
void split_number_string(const std::string &input, std::vector<double> &output) {
  
  // Get a space delimited string
  std::string delimiter = " ";
  std::string s = input;
  std::replace(s.begin(), s.end(), ',', ' ');

  double val;
  std::stringstream stream(s);
  while (stream >> val)
    output.push_back(val);
  
  // If the input is non-empty but the output is empty, that means
  // an invalid string was passed.
  if (!input.empty() && output.empty())
    vw_throw(ArgumentErr() << "Invalid value for the DEM spacing: " << input << "\n");
  
}

} // End namespace vw
