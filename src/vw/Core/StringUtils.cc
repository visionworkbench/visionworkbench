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
//#include <boost/algorithm/string.hpp>

namespace vw {


size_t string_replace(std::string &s, std::string const& find, std::string const& replace) {
  // Keep replacing the first match found until we don't find any more.
  size_t count = 0;
  while (true) {
    size_t pos = s.find(find);
    if (pos == std::string::npos)
      return count;
    s.replace(pos, find.length(), replace);
    ++count;
  }
}


} // End namespace vw
