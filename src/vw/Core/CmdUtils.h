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

/// \file CmdUtils.h
/// Utilities for invoking system commands

#ifndef __VW_CORE_CMDUTILS_H__
#define __VW_CORE_CMDUTILS_H__

namespace vw {

  // Execute a command and capture its output in a string
  std::string exec_cmd(const char* cmd);

  // Find the full path to a program in a given search directory based on the
  // path to the calling program. If not there, find it using the system PATH.
  std::string program_path(std::string const& prog_name, std::string const& curr_exec_path);
  
  // Look up the full path to prog_name in PATH
  std::string find_executable_in_path(std::string const& prog_name);

}

#endif
