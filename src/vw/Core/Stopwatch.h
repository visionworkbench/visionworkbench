#ifndef __VW_CORE_STOPWATCH_H__
#define __VW_CORE_STOPWATCH_H__

// System includes
#include <map>
#include <string>

// BOOST includes
#include <boost/thread/mutex.hpp>
#include <boost/shared_ptr.hpp>


namespace vw {

  // Stopwatch measures time elapsed between calls to start() and stop()
  
  class Stopwatch {
    struct data {
      unsigned long long m_total_elapsed; // in microseconds
      unsigned long long m_last_start;    // from Stopwatch::microtime
      unsigned long m_startdepth;
      unsigned long m_numstops;
      mutable boost::mutex m_mutex;
      data() :  m_total_elapsed(0),
                m_last_start(0),
                m_startdepth(0),
                m_numstops(0) {}
    };
    
    boost::shared_ptr<data> m_data;
    
  public:
    static unsigned long long microtime();
    
    Stopwatch() : m_data(new data()) {}
    
    void start() {
      boost::mutex::scoped_lock lock(m_data->m_mutex);
      if (!(m_data->m_startdepth++)) {
        m_data->m_last_start= microtime();
      }
    }

    void stop() {
      boost::mutex::scoped_lock lock(m_data->m_mutex);
      if (!--(m_data->m_startdepth)) {
        m_data->m_numstops++;
        m_data->m_total_elapsed += microtime() - m_data->m_last_start;
      }
    }

    bool is_running() const {
      return m_data->m_startdepth != 0;
    }

    unsigned long long elapsed_microseconds() const {
      boost::mutex::scoped_lock lock(m_data->m_mutex);
      return m_data->m_total_elapsed;
    }

    double elapsed_seconds() const {
      return (double)elapsed_microseconds() / 1000000.;
    }

    unsigned long num_stops() const {
      return m_data->m_numstops;
    }
    
  }; // class Stopwatch

  // StopwatchSet is a named set of Stopwatches
  class StopwatchSet {
    mutable boost::mutex m_mutex;
    unsigned long long m_construction_time;

    std::map<std::string, Stopwatch> m_stopwatches;

  public:
    StopwatchSet() : m_construction_time(Stopwatch::microtime()) {}

    // Find or create stopwatch named "name"
    Stopwatch get(const std::string &name) {
      boost::mutex::scoped_lock lock(m_mutex);
      return m_stopwatches[name];
    }

    unsigned long long elapsed_microseconds_since_construction() const {
      return Stopwatch::microtime()-m_construction_time;
    }

    double elapsed_seconds_since_construction() const {
      return (double)elapsed_microseconds_since_construction() / 1000000.;
    }

    // Returns a copy of the map between names and stopwatches
    // (Copy is for thread safety)
    std::map<std::string, Stopwatch> get_stopwatches() const {
      boost::mutex::scoped_lock lock(m_mutex);
      return m_stopwatches;
    }

    std::string report() const;
  }; // class StopwatchSet


  // 
  // Simple functions for global StopwatchSet
  //

  // Return the global StopwatchSet
  StopwatchSet *global_stopwatch_set();
  
  // Find or create names stopwatch from the global StopwatchSet
  inline Stopwatch stopwatch_get(const std::string &name) { return global_stopwatch_set()->get(name); }
  
  // Start the named stopwatch from the global StopwatchSet
  inline void stopwatch_start(const std::string &name) { stopwatch_get(name).start(); }
  
  // Stop the named stopwatch from the global StopwatchSet
  inline void stopwatch_stop(const std::string &name) { stopwatch_get(name).stop(); }

  //
  // ScopedWatch starts a Stopwatch on construction and stops is
  //  on destruction
  //
  
  class ScopedWatch {
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
    ~ScopedWatch() {
      m_stopwatch.stop();
    }
  }; // class ScopedWatch
  
} // namespace vw

#endif



  
