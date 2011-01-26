#include <vw/Core/System.h>
#include <vw/Core/Cache.h>
#include <vw/Core/Log.h>
#include <vw/Core/Settings.h>
#include <vw/Core/Stopwatch.h>

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
