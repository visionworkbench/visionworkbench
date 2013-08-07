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

/// \file System.h
/// All the Core singletons are consolidated here, to make their initialization
/// more deterministic (since some of them depend on each other)

#ifndef __VW_CORE_SYSTEM_H__
#define __VW_CORE_SYSTEM_H__

namespace vw {

  class Cache;
  class Log;
  class Settings;
  class StopwatchSet;

  // This cache is used by default for all new BlockImageView<>'s such as
  // DiskImageView<>.
  Cache& vw_system_cache();

  // You should *always* use this method if you want to access Vision Workbench
  // system log, where all Vision Workbench log messages go.  For example:
  //     vw_log().console_log() << "Some text\n";
  Log& vw_log();

  // Global instance of Settings
  Settings& vw_settings();

  // Global instance of StopwatchSet
  StopwatchSet& vw_stopwatch_set();
}

#endif
