#include <vw/Core/System.h>
#include <vw/Core/Cache.h>
#include <vw/Core/Log.h>
#include <vw/Core/Settings.h>
#include <vw/Core/Stopwatch.h>

namespace {
  vw::RunOnce core_once = VW_RUNONCE_INIT;
        vw::Cache*      cache_ptr = 0;
     vw::Settings*   settings_ptr = 0;
          vw::Log*        log_ptr = 0;
  vw::StopwatchSet* stopwatch_ptr = 0;

  void init_core() {
         settings_ptr = new vw::Settings();
              log_ptr = new vw::Log();
            cache_ptr = new vw::Cache(settings_ptr->system_cache_size());
        stopwatch_ptr = new vw::StopwatchSet();
  }
}

namespace vw {
  Cache& vw_system_cache() {
    core_once.run(init_core);
    return *cache_ptr;
  }
  Log& vw_log() {
    core_once.run(init_core);
    return *log_ptr;
  }
  Settings& vw_settings() {
    core_once.run(init_core);
    return *settings_ptr;
  }
  StopwatchSet& vw_stopwatch_set() {
    core_once.run(init_core);
    return *stopwatch_ptr;
  }
}
