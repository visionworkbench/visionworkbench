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

#include <cstdio>
#include <string>
#include <stdexcept>
#include <iostream>
#include <vw/Core/CmdUtils.h>

namespace vw {

  // Execute a command and capture its output in a string
  std::string exec_cmd(const char* cmd) {
    const int buf_len = 1024;
    char buffer[buf_len];
    std::string result = "";
    FILE* pipe = popen(cmd, "r");
    if (!pipe)
      throw std::runtime_error("popen() failed!");
    try {
      while (!feof(pipe)) {
        if (fgets(buffer, buf_len, pipe) != NULL)
          result += buffer;
      }
    } catch (...) {
      pclose(pipe);
      throw;
    }
    pclose(pipe);
    return result;
  }

  // Look up the full path to prog_name in PATH.
  // TODO: There is a boost function for this, but in a header
  // which we don't ship.
  std::string find_executable_in_path(std::string const& prog_name){

    // Use the 'which' command
    std::string cmd = "which " + prog_name;
    std::string ans = exec_cmd(cmd.c_str());

    // Wipe the newline
    if (!ans.empty() && ans[ans.length()-1] == '\n') 
      ans.erase(ans.length()-1); 
    
    if ( ans.find("not found") != std::string::npos )
      throw std::runtime_error("Could not find path to " + prog_name);      
    
    return ans;
  }
}  
