// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// System includes
#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>

// Time
#ifdef WIN32
#include <Windows.h>
#else // For Linux
#include <sys/time.h>
#endif 

// BOOST includes
#include <boost/thread/once.hpp>

// VisionWorkbench includes
#include <vw/Core/Debugging.h>

// Self
#include <vw/Core/Stopwatch.h>

using std::endl;
using std::map;
using std::pair;
using std::string;
using std::vector;

namespace vw {

  // Stopwatch
  unsigned long long Stopwatch::microtime() {
#ifdef WIN32
    return ((long long)GetTickCount() * 1000); 
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return ((long long) tv.tv_sec * (long long) 1000000 +
            (long long) tv.tv_usec);
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
    
    for (unsigned int i= 0; i< sorted.size(); i++) {
      const std::string &name(sorted[i].first);
      const Stopwatch &sw(sorted[i].second);

      double elapsed = sw.elapsed_seconds();
      int n = sw.num_stops();
      
      out << elapsed;
      if (n) out << " (avg " << elapsed / n << " x " << n << ")";
      out << ": " << name;
      
      if (sw.is_running()) out << " (still running!)";
      out << endl;
    }

    return out.str();
  }
  
  // Global StopwatchSet
  
  namespace GlobalStopwatchSet {
    static StopwatchSet *_g_global_stopwatch_set;
    void _create() {
      _g_global_stopwatch_set = new StopwatchSet();
    }
  };

  // Return the global StopwatchSet
  //
  // We want to be able to use this during globals and statics construction, so we cannot
  // assume that constructors in this file have already been called
  
  StopwatchSet *global_stopwatch_set() {
    static RunOnce once = VW_RUNONCE_INIT;
    once.run(GlobalStopwatchSet::_create);
    return GlobalStopwatchSet::_g_global_stopwatch_set;
  }

}; // namespace vw
