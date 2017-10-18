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


#include <vw/Core/Stopwatch.h>
#include <vw/Core/System.h>
#include <vw/Core/Exception.h>

// Time
#ifdef WIN32
#include <Windows.h>
#else // For Linux
#include <sys/time.h>
#endif

#include <iostream>
#include <vector>
#include <algorithm>

using std::endl;
using std::map;
using std::pair;
using std::string;
using std::vector;

namespace vw {

  // Stopwatch
  uint64 Stopwatch::microtime(bool use_cpu_time) {
#ifdef WIN32
    if (use_cpu_time) {
      // TODO: Find the analogous function for "clock()" for windows
      throw NoImplErr();
    } else {
      return ((uint64)GetTickCount() * 1000);
    }
#else
    if (use_cpu_time) {
      // XXX: This provides no useful information. On any linux system after
      // 2003 or so, the clock rate is not a constant.
      return ((uint64)clock() * 1000000 / CLOCKS_PER_SEC);
    } else {
      struct timeval tv;
      gettimeofday(&tv, NULL);
      return ((uint64)tv.tv_sec) * 1000000 + tv.tv_usec;
    }
#endif
  }

  bool pair_string_stopwatch_elapsed_gt(const pair<string, Stopwatch> &a,
                                        const pair<string, Stopwatch> &b)
  {
    return a.second.elapsed_microseconds() > b.second.elapsed_microseconds();
  }

  // StopwatchSet
  string StopwatchSet::report() const {
    Mutex::Lock lock(m_mutex);

    vector<pair<string, Stopwatch> > sorted(m_stopwatches.begin(), m_stopwatches.end());
    sort(sorted.begin(), sorted.end(), pair_string_stopwatch_elapsed_gt);

    std::ostringstream out;
    out << "StopwatchSet lifetime: " << elapsed_seconds_since_construction() << endl;
    out << "Elapsed seconds for all stopwatches:" << endl;

    for (size_t i= 0; i< sorted.size(); i++) {
      const std::string &name(sorted[i].first);
      const Stopwatch &sw(sorted[i].second);

      double elapsed = sw.elapsed_seconds();
      uint32 n = sw.num_stops();

      out << elapsed;
      if (n) out << " (avg " << elapsed / n << " x " << n << ")";
      out << ": " << name;

      if (sw.is_running()) out << " (still running!)";
      out << endl;
    }

    return out.str();
  }

  Stopwatch stopwatch_get(const std::string &name) {
    return vw_stopwatch_set().get(name);
  }
}; // namespace vw
