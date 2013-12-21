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


/// \file Core/Debugging.h
///
/// Types and functions to assist in debugging code.
///
#ifndef __VW_CORE_DEBUGGING_H__
#define __VW_CORE_DEBUGGING_H__

#include <vw/Core/Log.h>

#include <string>

#include <boost/current_function.hpp>

#define VW_CURRENT_FUNCTION BOOST_CURRENT_FUNCTION

namespace vw {


  // *******************************************************************
  // Timing types and functions
  // *******************************************************************

  class Timer {
    std::string m_desc;
    MessageLevel m_level;
    std::string m_log_namespace;
#ifdef WIN32
    __int64 m_begin;
#else
    timeval m_begin;
#endif
  public:
    Timer( std::string const& desc, MessageLevel level=InfoMessage, std::string const& log_namespace = "console" );

    ~Timer();
  };

  inline void call_top(std::string prog_name, std::string tag){
    // Call the top command to find the memory usage at the
    // given time.
    std::string cmd = "echo " + tag + " $(top -u $(whoami) -b -n 1 | grep "
      + prog_name + ")";
    std::system(cmd.c_str());
  }

} // namespace vw

#endif  // __VW_CORE_DEBUGGING_H__
