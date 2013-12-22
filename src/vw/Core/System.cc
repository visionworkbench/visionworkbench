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


#include <vw/Core/System.h>
#include <vw/Core/Cache.h>
#include <vw/Core/Log.h>
#include <vw/Core/Settings.h>
#include <vw/Core/Stopwatch.h>
#include <vw/Core/RunOnce.h>

namespace {
  vw::RunOnce settings_once      = VW_RUNONCE_INIT;
  vw::RunOnce resize_once        = VW_RUNONCE_INIT;
  vw::RunOnce stopwatch_set_once = VW_RUNONCE_INIT;
  vw::RunOnce system_cache_once  = VW_RUNONCE_INIT;
  vw::RunOnce log_once           = VW_RUNONCE_INIT;

  vw::Settings     *settings_ptr      = 0;
  vw::StopwatchSet *stopwatch_set_ptr = 0;
  vw::Cache        *system_cache_ptr  = 0;
  vw::Log          *log_ptr           = 0;

  void init_settings() {
    settings_ptr = new vw::Settings();
  }

  void resize_cache() {
    if (system_cache_ptr->max_size() == 0)
      system_cache_ptr->resize(settings_ptr->system_cache_size());
  }

  void init_system_cache() {
    system_cache_ptr = new vw::Cache(0);
  }

  void init_stopwatch_set() {
    stopwatch_set_ptr = new vw::StopwatchSet();
  }

  void init_log() {
    log_ptr = new vw::Log();
  }
}

vw::Settings &vw::vw_settings() {
  system_cache_once.run( init_system_cache );
  settings_once.run( init_settings );
  return *settings_ptr;
}

vw::Cache &vw::vw_system_cache() {
  system_cache_once.run( init_system_cache );
  settings_once.run( init_settings );
  settings_ptr->reload_config();
  resize_once.run(resize_cache);
  return *system_cache_ptr;
}

vw::StopwatchSet &vw::vw_stopwatch_set() {
  stopwatch_set_once.run( init_stopwatch_set );
  return *stopwatch_set_ptr;
}

vw::Log &vw::vw_log() {
  log_once.run( init_log );
  return *log_ptr;
}
