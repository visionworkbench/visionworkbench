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

/// * A RunOnce class, implementing run-once semantics.  This must be
///   a POD class with an initializer macro VW_RUNONCE_INIT.  It must
///   implement a single method, run( void (*func)() ), which runs the
///   given function exactly once no matter how many times it is
///   called.  (This behavior is only defined for RunOnce objects
///   that are statically allocated at global or namespace scope and
///   statically initialized to VW_RUNONCE_INIT.)

#ifndef __VW_CORE_RUNONCE_H__
#define __VW_CORE_RUNONCE_H__

#include <boost/thread/once.hpp>

namespace vw {

  // --------------------------------------------------------------
  //                            RUNONCE
  // --------------------------------------------------------------

#define VW_RUNONCE_INIT { BOOST_ONCE_INIT }

  // A special POD class to enable safe library initialization.  You
  // should only define these objects at global or namespace scope,
  // and statically initialize them to VW_RUNONCE_INIT.
  struct RunOnce {
    boost::once_flag m_flag;

    inline void run( void (*func)() ) {
      boost::call_once( func, m_flag );
    }
  };

}

#endif//__VW_CORE_RUNONCE_H__
