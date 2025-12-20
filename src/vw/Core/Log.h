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


/// \file Core/Log.h
///
/// Aside from the standard basic ostream functionality, this set of
/// classes provides:
///
/// - Buffering of the log messages on a per-thread messages so that
///   messages form different threads are nicely interleaved in the
///   log.
///
/// Some notes on the behavior of the log.
///
/// - A new line in the logfile starts every time a newline character
///   appears at the end of a string of characters, or when you
///   explicitly add std::flush() to the stream of operators.

#ifndef __VW_CORE_LOG_H__
#define __VW_CORE_LOG_H__

#include <vw/Core/Features.h>
#include <vw/Core/Thread.h>
#include <vw/Core/System.h>

// Boost Headers
#include <boost/numeric/conversion/cast.hpp>

// STD Headers
#include <string>
#include <list>
#include <vector>
#include <map>
#include <iostream>

// STD Headers for defining our own ostream subclass
#include <streambuf>
#include <ostream>
#include <fstream>

// For stringify
#include <sstream>

namespace vw {

  // ----------------------------------------------------------------
  // Debugging output types and functions
  // ----------------------------------------------------------------

  // Lower number -> higher priority
  enum MessageLevel {
    NoMessage           = -1,
    ErrorMessage        =  0,
    WarningMessage      = 10,
    InfoMessage         = 20,
    DebugMessage        = 30,
    VerboseDebugMessage = 40,
    EveryMessage        = 100
  };

  /// \cond INTERNAL

  // Forward declarations for internal stream classes
  template<typename CharT, typename traits = std::char_traits<CharT> > class NullOutputStream;
  template<typename CharT, typename traits = std::char_traits<CharT> > class MultiOutputStream;
  template<class CharT, class traits = std::char_traits<CharT> > class PerThreadBufferedStream;

  // Some handy typedefs
  //
  // These are made to be lower case names to jive with the C++ std
  // library naming convention for streams (i.e. std::cin, std::cout,
  // etc.)
  typedef NullOutputStream<char> null_ostream;
  typedef MultiOutputStream<char> multi_ostream;

  /// \endcond

  template<typename T>
  inline std::string stringify(const T& x)
  {
    std::ostringstream o;
    if (!(o << x))
      return std::string("[failed to stringify ") + typeid(x).name() + "]";
    return o.str();
  }

  class LogRuleSet {
    // The ruleset determines what log messages are sent to the VW system log file
    typedef std::pair<int, std::string> rule_type;
    typedef std::list<rule_type> rules_type;
    rules_type m_rules;
    Mutex m_mutex;

    // Help functions
    bool has_leading_wildcard( std::string const& exp );
    std::string after_wildcard( std::string const& exp );

  public:
    // by default, the LogRuleSet is set up to pass "console" messages
    // at level vw::InfoMessage or higher priority.
    LogRuleSet();
    virtual ~LogRuleSet();

    // Ensure Copyable semantics (Mutex is not copyable, so we create our own
    // mutex and copy the internal rules)
    LogRuleSet( LogRuleSet const& copy_log);
    LogRuleSet& operator=( LogRuleSet const& copy_log);

    void add_rule(int log_level, std::string const& log_namespace);
    void clear();

    // You can overload this method from a subclass to change the
    // behavior of the LogRuleSet.
    virtual bool operator() (int log_level, std::string const& log_namespace);
  };


  // -------------------------------------------------------
  //                         LogInstance
  // -------------------------------------------------------
  //
  class LogInstance : private boost::noncopyable {
    PerThreadBufferedStream<char>* m_log_stream;
    std::ostream *m_log_ostream_ptr;
    bool m_prepend_infostamp;
    LogRuleSet m_rule_set;

  public:

    // Initialize a log from a filename.  A new internal ofstream is
    // created to stream log messages to disk.
    LogInstance(std::string const& log_filename, bool prepend_infostamp = true);

    // Initialize a log using an already open stream.  Warning: The
    // log stores the stream by reference, so you MUST delete the log
    // object _before_ closing and de-allocating the stream.
    LogInstance(std::ostream& log_ostream, bool prepend_infostamp = true);

    ~LogInstance();

    /// This method return an ostream that you can write a log message
    /// to if the rule_set matches the log level and namespace
    /// provided.  Otherwise, a null ostream is returned.
    std::ostream& operator() (int log_level, std::string const& log_namespace="console");

    /// Access the rule set for this log object.
    LogRuleSet& rule_set() { return m_rule_set; }
  };


  // -------------------------------------------------------
  //                         Log
  // -------------------------------------------------------

  /// The system log class manages logging to the console and to files
  /// on disk.  It supports multiple open log streams, each with their
  /// own LogRuleSet.
  ///
  /// Important Note: You should access the system log using the vw::vw_log
  /// free function, which is a singleton instance of the system log class.  You
  /// should not need to create a log object yourself.
  class Log : private boost::noncopyable {

    // Pointers to various log instances that are currently being
    // managed by the system log.
    std::vector<boost::shared_ptr<LogInstance> > m_logs;
    boost::shared_ptr<LogInstance> m_console_log;

    // Member variables
    Mutex m_system_log_mutex;
    Mutex m_multi_ostreams_mutex;

    // The multi_ostream creates a single stream that delegates to its
    // child streams. We store one multi_ostream per thread, since
    // each thread will have a different set of output streams it is
    // currently accessing.
    //
    // TODO: Rather than manage this as a simple map (which will grow
    // quite large and hurt performance if there are many threads), we
    // should really use some sort of cache of shared_ptr's to
    // ostreams here.  The tricky thing is that once we return the
    // ostream, we don't know how long the thread will use it before
    // it can be safely de-allocated.
    std::map<vw::uint64, boost::shared_ptr<multi_ostream> > m_multi_ostreams;

  public:

    /// You should probably not create an instance of Log on your own
    /// using this constructor.  Instead, you can access a global
    /// instance of the log class using the static Log::system_log()
    /// method below.
    Log();
    ~Log();

    /// The call operator returns a subclass of the basic_ostream
    /// object, which is suitable for use with the C++ << operator.
    /// The returned stream object proxy's for the various log streams
    /// being managed by the system log that match the log_level and
    /// log_namespace.
    std::ostream& operator() (int log_level, std::string const& log_namespace="console");

    /// Add a stream to the Log manager.  You may optionally specify a
    /// LogRuleSet.
    void add(std::ostream &stream, LogRuleSet rule_set = LogRuleSet(), bool prepend_infostamp = true);

    // Add an already existing LogInstance to the system log manager.
    void add(boost::shared_ptr<LogInstance> log);

    /// Reset the System Log; closing all of the currently open Log
    /// streams.
    void clear();

    /// Return a reference to the console LogInstance.
    LogInstance& console_log();

    /// Set the output stream and LogRuleSet for the console log
    /// instance.  This can be used to redirect the console output to
    /// a file, for example.
    void set_console_stream(std::ostream& stream, LogRuleSet rule_set = LogRuleSet(), bool prepend_infostamp = true);

    /// A mostly non-locking check to determine if this Log object has any
    /// stream that is open for a requested log level and namespace.
    ///
    /// LogRuleSets still lock on access of their operator().
    /// We're also using shared_ptr, which lock on access.
    bool is_enabled( int log_level = vw::InfoMessage,
                     std::string const& log_namespace="console" );
  };

  /// The vision workbench logging operator.  Use this to generate a
  /// message in the system log using the given log_level and
  /// log_namespace.
  std::ostream& vw_out( int log_level = vw::InfoMessage,
                        std::string const& log_namespace = "console" );

  /// A macro for vw_out that allows for lazy evaluation. It achieves
  /// this by putting a conditional around the stream operation.
#define VW_OUT(...) if(::vw::vw_log().is_enabled(__VA_ARGS__)) ::vw::vw_out(__VA_ARGS__)

  /// Deprecated: Set the debug level for the system console log.  You
  /// can exercise much more fine grained control over the system log
  /// by manipulating the Log::system_log().console_log().rule_set().
  void set_debug_level( int log_level ) VW_DEPRECATED;

  /// Deprecated: Set the output stream for the system console log to
  /// an arbitrary C++ ostream.  You should use
  /// Log::system_log().set_console_stream() instead.
  void set_output_stream( std::ostream& stream ) VW_DEPRECATED;

} // namespace vw

#endif // __VW_CORE_LOG_H__
