#include <vw/Core/System.h>
#include <vw/Core/Cache.h>
#include <vw/Core/Log.h>
#include <vw/Core/Settings.h>
#include <vw/Core/Stopwatch.h>

#define FUNC(Name, Type, Code)                    \
  namespace {                                     \
    vw::RunOnce Name ## _once = VW_RUNONCE_INIT;  \
    Type* Name ## _ptr = 0;                       \
    void init_ ## Name() {                        \
      Name ## _ptr = new Type Code;               \
    }                                             \
  }                                               \
  Type& vw::vw_ ## Name() {                       \
    Name ## _once.run(init_ ## Name);             \
    return *Name ## _ptr;                         \
  }

FUNC(     settings, vw::Settings,     ());
FUNC(          log, vw::Log,          ());
FUNC( system_cache, vw::Cache,        (vw::vw_settings().system_cache_size()));
FUNC(stopwatch_set, vw::StopwatchSet, ());
