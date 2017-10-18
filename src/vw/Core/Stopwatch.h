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


#ifndef __VW_CORE_STOPWATCH_H__
#define __VW_CORE_STOPWATCH_H__

#include <vw/config.h>

#include <vw/Core/Thread.h>
#include <vw/Core/FundamentalTypes.h>

#include <map>
#include <string>

#include <boost/smart_ptr/shared_ptr.hpp>

/// \file Stopwatch.h Contains some simple stopwatch classes.

namespace vw {

  /// Stopwatch measures time elapsed between calls to start() and stop()

  class VW_API Stopwatch {
    struct data {
      uint64 m_total_elapsed; // in microseconds
      uint64 m_last_start;    // from Stopwatch::microtime
      uint32 m_startdepth;
      uint32 m_numstops;
      mutable Mutex m_mutex;
      data() :  m_total_elapsed(0),
                m_last_start(0),
                m_startdepth(0),
                m_numstops(0) {}
    };

    boost::shared_ptr<data> m_data;
    bool m_use_cpu_time;

  public:
    static uint64 microtime(bool use_cpu_time = false);

    Stopwatch(bool use_cpu_time = false) : m_data(new data()), m_use_cpu_time(use_cpu_time) {}

    void start() {
      Mutex::Lock lock(m_data->m_mutex);
      if (!(m_data->m_startdepth++)) {
        m_data->m_last_start= microtime(m_use_cpu_time);
      }
    }

    void stop() {
      Mutex::Lock lock(m_data->m_mutex);
      if (!--(m_data->m_startdepth)) {
        m_data->m_numstops++;
        m_data->m_total_elapsed += microtime(m_use_cpu_time) - m_data->m_last_start;
      }
    }

    bool is_running() const {
      return m_data->m_startdepth != 0;
    }

    uint64 elapsed_microseconds() const {
      Mutex::Lock lock(m_data->m_mutex);
      return m_data->m_total_elapsed;
    }

    double elapsed_seconds() const {
      return (double)elapsed_microseconds() / 1000000.;
    }

    uint32 num_stops() const {
      return m_data->m_numstops;
    }

  }; // class Stopwatch

  // StopwatchSet is a named set of Stopwatches
  class VW_API StopwatchSet {
    mutable Mutex m_mutex;
    uint64 m_construction_time;

    std::map<std::string, Stopwatch> m_stopwatches;

  public:
    StopwatchSet() :
      m_construction_time(Stopwatch::microtime()) {}

    // Find or create stopwatch named "name"
    Stopwatch get(const std::string &name) {
      Mutex::Lock lock(m_mutex);
      return m_stopwatches[name];
    }

    uint64 elapsed_microseconds_since_construction() const {
      return Stopwatch::microtime()-m_construction_time;
    }

    double elapsed_seconds_since_construction() const {
      return (double)elapsed_microseconds_since_construction() / 1000000.;
    }

    // Returns a copy of the map between names and stopwatches
    // (Copy is for thread safety)
    std::map<std::string, Stopwatch> get_stopwatches() const {
      Mutex::Lock lock(m_mutex);
      return m_stopwatches;
    }

    std::string report() const;
  }; // class StopwatchSet


  //
  // Simple functions for global StopwatchSet
  //

  // Find or create names stopwatch from the global StopwatchSet
  Stopwatch VW_API stopwatch_get(const std::string &name);

  // Start the named stopwatch from the global StopwatchSet
  inline void stopwatch_start(const std::string &name) { stopwatch_get(name).start(); }

  // Stop the named stopwatch from the global StopwatchSet
  inline void stopwatch_stop(const std::string &name) { stopwatch_get(name).stop(); }

  //
  // ScopedWatch starts a Stopwatch on construction and stops is
  //  on destruction
  //

  class VW_API ScopedWatch {
  private:
    Stopwatch m_stopwatch;
  public:
    ScopedWatch(const Stopwatch &stopwatch)
      : m_stopwatch(stopwatch) {
      m_stopwatch.start();
    }
    // Use the named Stopwatch from the global StopwatchSet
    ScopedWatch(const std::string &name)
      : m_stopwatch(stopwatch_get(name)) {
      m_stopwatch.start();
    }
    // This must be overloaded to prevent ambiguity
    ScopedWatch(const char * name)
      : m_stopwatch(stopwatch_get(name)) {
      m_stopwatch.start();
    }
    ~ScopedWatch() {
      m_stopwatch.stop();
    }
  }; // class ScopedWatch

} // namespace vw

#endif
