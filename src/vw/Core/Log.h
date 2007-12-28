// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__

/// \file Core/Log.h
/// 
/// Aside from the standard basic ostream functionality, this set of
/// classes provides:
///
/// - Buffering of the log messages on a per-thread messages so that 
///   messages form different threads are nicelely interleaved in the 
///   log.
///
/// Some notes on the behavior of the log.
///
/// - A new line in the logfile starts every time a newline character
///   appears at the end of a string of characters, or when you
///   exlicitly add std::flush() to the stream of operators.

#ifndef __VW_CORE_LOG_H__
#define __VW_CORE_LOG_H__

#include <vw/Core/Thread.h>
#include <vw/Core/Exception.h>

// Boost Headers
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string.hpp>

// STD Headers
#include <string>
#include <list>
#include <vector>

// STD Headers for defininog our own ostream subclass
#include <streambuf>
#include <ostream>
#include <fstream>

namespace vw {

  // ----------------------------------------------------------------
  // Debugging output types and functions
  // ----------------------------------------------------------------

  enum MessageLevel {
    ErrorMessage = 0,
    WarningMessage = 10,
    InfoMessage = 20,
    DebugMessage = 30,
    VerboseDebugMessage = 40,
    EveryMessage = 100
  };
  
  // ----------------------------------------------------------------
  //                     Utility Streams
  // ----------------------------------------------------------------
  //
  // These classes provide a basic NULL ostream and an ostream that
  // can take data and re-stream it to multiple sub-streams.

  namespace {
    // Null Output Stream
    template<typename CharT, typename traits = std::char_traits<CharT> >
    class NullOutputBuf : public std::basic_streambuf<CharT, traits> {
      typedef typename std::basic_streambuf<CharT, traits>::int_type int_type;
      virtual int_type overflow(int_type c) { return traits::not_eof(c); }
    };
    
    template<typename CharT, typename traits>
    class NullOutputStreamInit {
      NullOutputBuf<CharT, traits> m_buf;      
    public:
      NullOutputBuf<CharT, traits>* buf() { return &m_buf; }
    };

    template<typename CharT, typename traits = std::char_traits<CharT> >
    class NullOutputStream : private virtual NullOutputStreamInit<CharT, traits>,
                             public std::basic_ostream<CharT, traits> {
    public:
      NullOutputStream() : NullOutputStreamInit<CharT, traits>(),
                           std::basic_ostream<CharT,traits>(NullOutputStreamInit<CharT, traits>::buf()) {}
    };

    // Multi output stream
    template<typename CharT, typename traits = std::char_traits<CharT> >
    class MultiOutputBuf : public std::basic_streambuf<CharT, traits> {
      typedef typename std::basic_streambuf<CharT, traits>::int_type int_type;
      typedef std::vector<std::basic_ostream<CharT, traits>* > stream_container;
      typedef typename stream_container::iterator stream_iterator;
      stream_container m_streams;
      
    protected:
      virtual std::streamsize xsputn(const CharT* sequence, std::streamsize num) {
        stream_iterator current = m_streams.begin();
        stream_iterator end = m_streams.end();
        for(; current != end; ++current) 
          (*current)->write(sequence, num);
        return num;
      }
      
      virtual int_type overflow(int_type c) {
        stream_iterator current = m_streams.begin();
        stream_iterator end = m_streams.end();
        
        for(; current != end; ++current) 
          (*current)->put(c);
        return c;
      }

    public:
      void add(std::basic_ostream<CharT, traits>& stream) { m_streams.push_back(&stream); }
      void remove(std::basic_ostream<CharT, traits>& stream) {
        stream_iterator pos = std::find(m_streams.begin(),m_streams.end(), &stream);
        if(pos != m_streams.end()) 
          m_streams.erase(pos);
      }
      void clear() { m_streams.clear(); }
    };

    template<typename CharT, typename traits>
    class MultiOutputStreamInit {
      MultiOutputBuf<CharT, traits> m_buf;
    public:
      MultiOutputBuf<CharT, traits>* buf() { return &m_buf; }
    };

    template<typename CharT, typename traits = std::char_traits<CharT> >
    class MultiOutputStream : private MultiOutputStreamInit<CharT, traits>, 
                              public std::basic_ostream<CharT, traits> {
    public:
      MultiOutputStream() : MultiOutputStreamInit<CharT, traits>(), 
                            std::basic_ostream<CharT, traits>(MultiOutputStreamInit<CharT, traits>::buf()) {}
      void add(std::basic_ostream<CharT, traits>& str) { MultiOutputStreamInit<CharT, traits>::buf()->add(str); }
      void remove(std::basic_ostream<CharT, traits>& str) { MultiOutputStreamInit<CharT, traits>::buf()->remove(str); }
      void clear() { MultiOutputStreamInit<CharT, traits>::buf()->clear(); }
    };
    
    // Some handy typedefs
    //
    // These are made to be lower case names to jive with the C++ std
    // library naming convertion for streams (i.e. std::cin, std::cout,
    // etc.)
    typedef NullOutputStream<char> null_ostream;
    typedef MultiOutputStream<char> multi_ostream;
  }

  // In order to create our own C++ streams compatible ostream object,
  // we must first define a subclass of basic_streambuf<>, which
  // handles stream output on a character by character basis.  This is
  // not the most elegent block of code, but this seems to be the
  // "approved" method for defining custom behaviour in a subclass of
  // basic_ostream<>.
  //
  template<class CharT, class traits = std::char_traits<CharT> >
  class LogStreamBuf : public std::basic_streambuf<CharT, traits> {

    typedef typename std::basic_streambuf<CharT, traits>::int_type int_type;

    // Characters are buffered is vectors until a newline appears at
    // the end of a line of input or flush() is called.  These vectors
    // are indexed in a std::map by thread id.
    // 
    // TODO: This map could grow quite large if a program spawns (and
    // logs to) many, many threads.  We should think carefully about
    // cleaning up this map structure from time to time.
    typedef std::vector<CharT> buffer_type;
    typedef std::map<int, buffer_type> lookup_table_type;
    lookup_table_type m_buffers;

    std::basic_streambuf<CharT, traits>* m_out;
    Mutex m_mutex;

    // This method is called when a single character is fed to the
    // streambuf.  In practice, characters are fed in batches using
    // xputn() below.
    virtual int_type overflow(int_type c) {
      Mutex::Lock lock(m_mutex);
      if(!traits::eq_int_type(c, traits::eof())) {
        m_buffers[ Thread::id() ].push_back(c);
      }
      return traits::not_eof(c);
    }

    virtual std::streamsize xsputn(const CharT* s, std::streamsize num) {
      {
        Mutex::Lock lock(m_mutex);
        std::copy(s, s + num, std::back_inserter<buffer_type>( m_buffers[ Thread::id() ] ));
      }

      // This is a bit of a hack that forces a sync whenever the
      // character string *ends* with a newline, thereby flushing the
      // buffer and printing a line to the log file.
      if ( *(s+num-1) == '\n' || *(s+num-1) == '\r' )
        sync();

      return num;
    }

    virtual int sync() {
      Mutex::Lock lock(m_mutex);
      if(!m_buffers[ Thread::id() ].empty() && m_out) {
        m_out->sputn(&m_buffers[ Thread::id() ][0], static_cast<std::streamsize>(m_buffers[ Thread::id() ].size()));
        m_buffers[ Thread::id() ].clear();
      }
      return 0;
    }

  public:
    LogStreamBuf() : m_out(NULL), m_buffers() {}
    ~LogStreamBuf() { sync(); }

    void init(std::basic_streambuf<CharT,traits>* out) { m_out = out; }
  };

  // The order with which the base classes are initialized in
  // LogStream is not fully defined unless we inherit from this as a
  // pure virtual base class.  This gives us the extra wiggle room we
  // need to connect two properly initialized LogStreamBuf and
  // LogStream objects.
  template<class CharT, class traits = std::char_traits<CharT> >
  class LogStreamBufInit {
    LogStreamBuf<CharT, traits> m_buf;    
  public:
    LogStreamBuf<CharT, traits>* buf() { return &m_buf; }
  };

  // Aside from some tricky initialization semantics, this subclass of
  // basic_ostream is actaully fairly simple.  It passes along
  // characters to the LogStreamBuf, which does the actual interesting
  // stuff.
  //
  // You can use the ctor or the set_stream method to pass in any C++
  // ostream (e.g. std::cout or a std::ofstream) to be the ultimate
  // recipient of the characters that are fed through this stream,
  // which acts as an intermediary, queuing characters on a per thread
  // basis and ensuring thread safety.
  template<class CharT, class traits = std::char_traits<CharT> >
  class LogStream : private virtual LogStreamBufInit<CharT, traits>,
                    public std::basic_ostream<CharT, traits> {
  public:    
    // No stream specified.  Will swallow characters until one is set
    // using set_stream().
    LogStream() : LogStreamBufInit<CharT,traits>(),
                  std::basic_ostream<CharT, traits>(LogStreamBufInit<CharT,traits>::buf()) {
    }


    LogStream(std::basic_ostream<CharT, traits>& out) : LogStreamBufInit<CharT,traits>(),
                                                        std::basic_ostream<CharT, traits>(LogStreamBufInit<CharT,traits>::buf()) {
      LogStreamBufInit<CharT,traits>::buf()->init(out.rdbuf());
    }
    
    void set_stream(std::basic_ostream<CharT, traits>& out) {
      LogStreamBufInit<CharT,traits>::buf()->init(out.rdbuf());
    }
  };

  class LogRuleSet {
    // The ruleset determines what log messages are sent to the VW system log file
    typedef std::pair<int, std::string> rule_type;
    typedef std::list<rule_type> rules_type;
    rules_type m_rules;

  public:
    
    // By default, the LogRuleSet is set up to pass "console" messages
    // at level vw::InfoMessage or higher.
    LogRuleSet() {
      m_rules.push_back(rule_type(vw::InfoMessage, "console"));      
    }

    void add_rule(std::string log_namespace, int log_level) {
      m_rules.push_back(rule_type(log_level, boost::to_lower_copy(log_namespace)));
    }
        
    // You can overload this method froma subclass to change the
    // behavior of the LogRuleSet.
    virtual bool operator() (std::string log_namespace, int log_level) {
      for (rules_type::iterator it = m_rules.begin(); it != m_rules.end(); ++it) {
        
        std::string ns = boost::to_lower_copy((*it).second);
        // Pass through rule for complete wildcard
        if ( vw::EveryMessage == (*it).first && ns == "any" )
          return true;
        
        // Pass through if the level matches and the namespace is a wildcard
        if ( log_level <= (*it).first && ns == "any" )
          return true;

        // Pass through if the level and namepace match
        if ( log_level <= (*it).first && ns == log_namespace )
           return true;
      }

      // We reach this line if all of the rules have failed, in
      // which case we return a NULL stream, which will result in
      // nothing being logged.
      return false;
    }
  };

  class Log {
    std::ostream *m_log_ostream_ptr;
    LogStream<char> m_log_stream;
    bool m_prepend_infostamp;
    LogRuleSet m_rule_set;

    // Ensure non-copyable semantics
    Log( Log const& );
    Log& operator=( Log const& );

  public:

    // Initialize a log from a filename.  A new internal ofstream is
    // created to stream log messages to disk.
    Log(std::string log_filename, bool prepend_infostamp = true) : m_prepend_infostamp(prepend_infostamp) {
      m_log_ostream_ptr = new std::ofstream(log_filename.c_str());
      if (! static_cast<std::ofstream*>(m_log_ostream_ptr)->is_open())
        vw_throw(IOErr() << "Could not open log file " << log_filename << " for writing.");
      m_log_stream.set_stream(*m_log_ostream_ptr);
    }

    // Initialize a log using an already open stream.  Warning: The
    // log stores the stream by reference, so you MUST delete the log
    // object _before_ closing and de-allocating the stream.
    Log(std::ostream& log_ostream, bool prepend_infostamp = true) : m_log_stream(log_ostream), 
                                                                    m_log_ostream_ptr(NULL),
                                                                    m_prepend_infostamp(prepend_infostamp) {}

    ~Log() {
      m_log_stream.set_stream(std::cout);
      if (m_log_ostream_ptr) 
        delete static_cast<std::ofstream*>(m_log_ostream_ptr);
    }

    std::ostream& operator() (int level, std::string log_namespace="general") { 
      boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
      if (m_prepend_infostamp) 
        return m_log_stream << now << " {" << Thread::id() << "} [ " << log_namespace << " ] : ";
      else 
        return m_log_stream;
    }

    LogRuleSet& rule_set() { return m_rule_set; }
  };

  class SystemLog {

    // Member variables assoc. with periodically polling the log
    // configuration (logconf) file.
    long m_logconf_last_polltime;
    long m_logconf_last_modification;
    std::string m_logconf_filename;
    double m_logconf_poll_period;
    Mutex m_logconf_time_mutex;
    Mutex m_logconf_file_mutex;
    
    std::vector<boost::shared_ptr<Log> > m_logs;
    boost::shared_ptr<Log> m_console_log;

    // The multi_ostream creates a single stream that delegates to
    // its child streams.      
    multi_ostream m_multi_ostream;

    // Private methods for checking the logconf file
    void stat_logconf();
    void reload_logconf_rules();

    // Ensure non-copyable semantics
    SystemLog( SystemLog const& );
    SystemLog& operator=( Log const& );

  public:

    SystemLog() : m_logconf_last_polltime(0),
                  m_logconf_last_modification(0),
                  m_logconf_filename("/tmp/.vw_logconf"),
                  m_logconf_poll_period(5.0),
                  m_console_log(new Log(std::cout, false)) {}

    std::ostream& operator() (int level, std::string log_namespace="console");

    void set_logconf_filename(std::string filename) { m_logconf_filename = filename; }
    void set_logconf_poll_period(double period) { m_logconf_poll_period = period; }

    /// Add a stream to the SystemLog manager.  You may optionally specify a
    /// LogRuleSet.
    void add(std::ostream &stream, LogRuleSet rule_set = LogRuleSet()) {
      m_logs.push_back( boost::shared_ptr<Log>(new Log(stream)) );
    }

    // Add an already existing log to the system log manager.
    void add(boost::shared_ptr<Log> log) {
      m_logs.push_back( log );
    }

    /// Reset the System Log; closing all of the currently open Log
    /// streams.
    void clear() { m_logs.clear(); }
    
    void set_console_stream(std::ostream& stream, LogRuleSet rule_set = LogRuleSet()) {
      m_console_log = boost::shared_ptr<Log>(new Log(stream) );
      m_console_log->rule_set() = rule_set;
    }

    void set_console_rule_set(std::ostream& stream, LogRuleSet rule_set = LogRuleSet()) {
      m_console_log->rule_set() = rule_set;
    }

    // Static instance used for accessing a single instance of the
    // System log.
    static SystemLog& system_log();
  };


  // Declaration of the standard VW I/O routine.
  std::ostream& vw_out( int log_level, std::string log_namespace = "console" );
  SystemLog& system_log();
  void set_debug_level( int log_level );
  void set_output_stream( std::ostream& stream );



} // namespace vw

#endif // __VW_CORE_LOG_H__
