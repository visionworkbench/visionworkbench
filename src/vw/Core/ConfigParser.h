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

/// \file ConfigParser.h Contains config parsing functions.


#ifndef __VW_CORE_CONFIGPARSER_H__
#define __VW_CORE_CONFIGPARSER_H__

#include <istream>

namespace vw {
  class Settings;

  // Parse a stream containing a config file,
  // and sets options through given settings object
  // throws on error, prints to cerr on warning
  void parse_config(std::basic_istream<char>& stream, vw::Settings&);

  // Parse a config file, and sets options through given settings object
  // throws on error, prints to cerr on warning
  void parse_config_file(const char* fn, vw::Settings&);
}

#endif
