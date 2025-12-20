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

#include <vw/config.h>
#include <vw/Core/Exception.h>
#include <vw/Core/Log.h>
#include <vw/Core/Settings.h>
#include <vw/Core/System.h>
#include <vw/Core/Thread.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <ctime>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include <streambuf>
#include <fstream>
#include <boost/numeric/conversion/cast.hpp>

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/foreach.hpp>

#ifdef WIN32
#define stat _stat
typedef struct _stat struct_stat;
#else
typedef struct stat struct_stat;
#endif


inline std::string
current_posix_time_string()
{
  char time_string[2048];
  time_t t = time(0);
  struct tm* time_struct = localtime(&t);
  strftime(time_string, 2048, "%F %T", time_struct);
  return std::string(time_string);
}

namespace vw {
  // ----------------------------------------------------------------
  //                     Utility Streams
  // ----------------------------------------------------------------
  //
  // These classes provide a basic NULL ostream and an ostream that
  // can take data and re-stream it to multiple sub-streams.

  // Null Output Stream
  template<typename CharT, typename traits>
  class NullOutputBuf : public std::basic_streambuf<CharT, traits> {
    typedef typename std::basic_streambuf<CharT, traits>::int_type int_type;
    virtual int_type overflow(int_type c) { return traits::not_eof(c); }
    virtual std::streamsize xsputn(const CharT* /*sequence*/, std::streamsize num) { return num; }
  };

  template<typename CharT, typename traits>
  class NullOutputStreamInit {
    NullOutputBuf<CharT, traits> m_buf;
  public:
    NullOutputBuf<CharT, traits>* buf() { return &m_buf; }
  };

  template<typename CharT, typename traits>
  class NullOutputStream : private virtual NullOutputStreamInit<CharT, traits>,
                           public std::basic_ostream<CharT, traits> {
  public:
    NullOutputStream() : NullOutputStreamInit<CharT, traits>(),
                         std::basic_ostream<CharT,traits>(NullOutputStreamInit<CharT, traits>::buf()) {}
  };

  // Multi output stream
  template<typename CharT, typename traits>
  class MultiOutputBuf : public std::basic_streambuf<CharT, traits> {
    typedef typename std::basic_streambuf<CharT, traits>::int_type int_type;
    typedef std::vector<std::basic_ostream<CharT, traits>* > stream_container;
    typedef typename stream_container::iterator stream_iterator;
    stream_container m_streams;
    Mutex m_mutex;

  protected:
    virtual std::streamsize xsputn(const CharT* sequence, std::streamsize num) {
      Mutex::Lock lock(m_mutex);
      stream_iterator current = m_streams.begin();
      stream_iterator end = m_streams.end();
      for(; current < end; ++current)
        (*current)->write(sequence, num);
      return num;
    }

    virtual int_type overflow(int_type c) {
      Mutex::Lock lock(m_mutex);
      stream_iterator current = m_streams.begin();
      stream_iterator end = m_streams.end();

      for(; current < end; current++)
        (*current)->put(static_cast<CharT>(c));
      return c;
    }

    // The sync function simply propogates the sync() request to the
    // child ostreams.
    virtual int sync() {
      Mutex::Lock lock(m_mutex);
      stream_iterator current = m_streams.begin();
      stream_iterator end = m_streams.end();

      for(; current < end; current++)
        (*current)->rdbuf()->pubsync();
      return 0;
    }

  public:
    void add(std::basic_ostream<CharT, traits>& stream) {
      Mutex::Lock lock(m_mutex);
      m_streams.push_back(&stream);
    }
    void remove(std::basic_ostream<CharT, traits>& stream) {
      Mutex::Lock lock(m_mutex);
      stream_iterator pos = std::find(m_streams.begin(),m_streams.end(), &stream);
      if(pos != m_streams.end())
        m_streams.erase(pos);
    }
    void clear() {
      Mutex::Lock lock(m_mutex);
      m_streams.clear();
    }
  };

  template<typename CharT, typename traits>
  class MultiOutputStreamInit {
    MultiOutputBuf<CharT, traits> m_buf;
  public:
    MultiOutputBuf<CharT, traits>* buf() { return &m_buf; }
  };

  template<typename CharT, typename traits>
  class MultiOutputStream : private MultiOutputStreamInit<CharT, traits>,
                            public std::basic_ostream<CharT, traits> {
  public:
    MultiOutputStream() : MultiOutputStreamInit<CharT, traits>(),
                            std::basic_ostream<CharT, traits>(MultiOutputStreamInit<CharT, traits>::buf()) {}
    void add(std::basic_ostream<CharT, traits>& str) { MultiOutputStreamInit<CharT, traits>::buf()->add(str); }
    void remove(std::basic_ostream<CharT, traits>& str) { MultiOutputStreamInit<CharT, traits>::buf()->remove(str); }
    void clear() { MultiOutputStreamInit<CharT, traits>::buf()->clear(); }
  };

  // In order to create our own C++ streams compatible ostream object,
  // we must first define a subclass of basic_streambuf<>, which
  // handles stream output on a character by character basis.  This is
  // not the most elegant block of code, but this seems to be the
  // "approved" method for defining custom behaviour in a subclass of
  // basic_ostream<>.
  template<class CharT, class traits>
  class PerThreadBufferedStreamBuf : public std::basic_streambuf<CharT, traits> {

    typedef typename std::basic_streambuf<CharT, traits>::int_type int_type;

    // Characters are buffered is vectors until a newline appears at
    // the end of a line of input or flush() is called.  These vectors
    // are indexed in a std::map by thread id.
    //
    // TODO: This map could grow quite large if a program spawns (and
    // logs to) many, many threads.  We should think carefully about
    // cleaning up this map structure from time to time.
    typedef std::vector<CharT> buffer_type;
    typedef std::map<vw::uint64, buffer_type> lookup_table_type;
    lookup_table_type m_buffers;

    std::basic_streambuf<CharT, traits>* m_out;
    Mutex m_mutex;

    // This method is called when a single character is fed to the
    // streambuf.  In practice, characters are fed in batches using
    // xputn() below.
    virtual int_type overflow(int_type c) {
      Mutex::Lock lock(m_mutex);
      buffer_type& buffer = m_buffers[ Thread::id() ];

      if(!traits::eq_int_type(c, traits::eof())) {
        buffer.push_back(static_cast<CharT>(c));
      }

      // If the last character is a newline or carrage return, then
      // we force a call to sync().
      if ( c == '\n' || c == '\r' )
        locked_sync(buffer);
      return traits::not_eof(c);
    }

    virtual std::streamsize xsputn(const CharT* s, std::streamsize num) {
      Mutex::Lock lock(m_mutex);
      buffer_type& buffer = m_buffers[ Thread::id() ];

      std::copy(s, s + num, std::back_inserter<buffer_type>( buffer ));

      // This is a bit of a hack that forces a sync whenever the
      // character string *ends* with a newline, thereby flushing the
      // buffer and printing a line to the log file.
      if ( buffer.size() > 0 ) {
        size_t last_char_position = buffer.size()-1;

        if ( buffer[last_char_position] == '\n' ||
             buffer[last_char_position] == '\r' )
          locked_sync(buffer);
      }
      return num;
    }

    // You must call this with the lock already held!
    int locked_sync(buffer_type& buffer) {
      if(!buffer.empty() && m_out ) {
        m_out->sputn(&buffer[0], boost::numeric_cast<std::streamsize>(buffer.size()));
        m_out->pubsync();
        buffer.clear();
      }
      return 0;
    }

    virtual int sync() {
      Mutex::Lock lock(m_mutex);
      if (m_buffers.count(Thread::id()))
        return locked_sync(m_buffers[ Thread::id() ]);
      return 0;
    }

  public:
    PerThreadBufferedStreamBuf() : m_buffers(), m_out(NULL) {}
    ~PerThreadBufferedStreamBuf() { sync(); }

    void init(std::basic_streambuf<CharT,traits>* out) { m_out = out; }
  };

  // The order with which the base classes are initialized in
  // PerThreadBufferedStream is not fully defined unless we inherit from this as a
  // pure virtual base class.  This gives us the extra wiggle room we
  // need to connect two properly initialized PerThreadBufferedStreamBuf and
  // PerThreadBufferedStream objects.
  template<class CharT, class traits = std::char_traits<CharT> >
  class PerThreadBufferedStreamBufInit {
    PerThreadBufferedStreamBuf<CharT, traits> m_buf;
  public:
    PerThreadBufferedStreamBuf<CharT, traits>* buf() {
      return &m_buf;
    }
  };

  // Aside from some tricky initialization semantics, this subclass of
  // basic_ostream is actually fairly simple.  It passes along
  // characters to the PerThreadBufferedStreamBuf, which does the
  // actual interesting stuff.
  //
  // You can use the ctor or the set_stream method to pass in any C++
  // ostream (e.g. std::cout or a std::ofstream) to be the ultimate
  // recipient of the characters that are fed through this stream,
  // which acts as an intermediary, queuing characters on a per thread
  // basis and ensuring thread safety.
  template<class CharT, class traits>
  class PerThreadBufferedStream : private virtual PerThreadBufferedStreamBufInit<CharT, traits>,
                                  public std::basic_ostream<CharT, traits> {
  public:
    // No stream specified.  Will swallow characters until one is set
    // using set_stream().
    PerThreadBufferedStream() : PerThreadBufferedStreamBufInit<CharT,traits>(),
                                std::basic_ostream<CharT, traits>(PerThreadBufferedStreamBufInit<CharT,traits>::buf()) {}


    PerThreadBufferedStream(std::basic_ostream<CharT, traits>& out) : PerThreadBufferedStreamBufInit<CharT,traits>(),
                                                                      std::basic_ostream<CharT, traits>(PerThreadBufferedStreamBufInit<CharT,traits>::buf()) {
      PerThreadBufferedStreamBufInit<CharT,traits>::buf()->init(out.rdbuf());
    }

    void set_stream(std::basic_ostream<CharT, traits>& out) {
      PerThreadBufferedStreamBufInit<CharT,traits>::buf()->init(out.rdbuf());
    }

  };

  static vw::null_ostream g_null_ostream;
}

// ---------------------------------------------------
// Basic stream support
// ---------------------------------------------------
std::ostream& vw::vw_out( int log_level, std::string const& log_namespace ) {
  return vw_log()(log_level, log_namespace);
}

void vw::set_debug_level( int log_level ) {
  vw_log().console_log().rule_set().add_rule(log_level, "console");
}

void vw::set_output_stream( std::ostream& stream ) {
  vw_log().set_console_stream(stream);
}

// ---------------------------------------------------
// LogInstance Methods
// ---------------------------------------------------
vw::LogInstance::LogInstance(std::string const& log_filename, bool prepend_infostamp) : m_prepend_infostamp(prepend_infostamp) {
  // Open file and place the insertion pointer at the end of the file (ios_base::ate)
  m_log_ostream_ptr = new std::ofstream(log_filename.c_str(), std::ios::app);
  if (! static_cast<std::ofstream*>(m_log_ostream_ptr)->is_open())
    vw_throw(IOErr() << "Could not open log file " << log_filename << " for writing.");

  *m_log_ostream_ptr << "\n\n" << "Vision Workbench log started at " << current_posix_time_string() << ".\n\n";

  m_log_stream = new PerThreadBufferedStream<char>(*m_log_ostream_ptr);
}

vw::LogInstance::LogInstance(std::ostream& log_ostream, bool prepend_infostamp) : m_log_stream(new PerThreadBufferedStream<char>(log_ostream)),
                                                                                  m_log_ostream_ptr(NULL),
                                                                                  m_prepend_infostamp(prepend_infostamp) {}

vw::LogInstance::~LogInstance() {
  m_log_stream->set_stream(std::cout);
  if (m_log_ostream_ptr)
    delete static_cast<std::ofstream*>(m_log_ostream_ptr);
  delete m_log_stream;
}

std::ostream& vw::LogInstance::operator() (int log_level, std::string const& log_namespace) {
  if (m_rule_set(log_level, log_namespace)) {
    if (m_prepend_infostamp)
      *m_log_stream << current_posix_time_string() << " {" << Thread::id() << "} [ " << log_namespace << " ] : ";
    switch (log_level) {
    case ErrorMessage:   *m_log_stream << "Error: ";   break;
    case WarningMessage: *m_log_stream << "Warning: "; break;
    default: break;
    }
    return *m_log_stream;
  } else {
    return g_null_ostream;
  }
}

vw::Log::Log() : m_console_log(new LogInstance(std::cout, false)) { }
vw::Log::~Log() {}

std::ostream& vw::Log::operator() (int log_level, std::string const& log_namespace) {
  // First, check to see if the rc file has been updated.
  // Reload the rulesets if it has.
  vw_settings().reload_config();

  {
    Mutex::Lock multi_ostreams_lock(m_multi_ostreams_mutex);

    // Check to see if we have an ostream defined yet for this thread.
    if(m_multi_ostreams.find( Thread::id() ) == m_multi_ostreams.end())
      m_multi_ostreams[ Thread::id() ] = boost::shared_ptr<multi_ostream>(new multi_ostream);

    boost::shared_ptr<multi_ostream>& stream = m_multi_ostreams[ Thread::id() ];

    // Reset and add the console log output...
    stream->clear();
    stream->add(m_console_log->operator()(log_level, log_namespace));

    // ... and the rest of the active log streams.
    std::vector<boost::shared_ptr<LogInstance> >::iterator iter = m_logs.begin();
    for (;iter != m_logs.end(); ++iter)
      stream->add((*iter)->operator()(log_level,log_namespace));

    return *stream;
  }
}

void
vw::Log::add(std::ostream &stream, LogRuleSet rule_set, bool prepend_infostamp) {
  Mutex::Lock lock(m_system_log_mutex);
  boost::shared_ptr<LogInstance> li(new LogInstance(stream, prepend_infostamp));
  li->rule_set() = rule_set;
  m_logs.push_back(li);
}

void
vw::Log::add(boost::shared_ptr<LogInstance> log) {
  Mutex::Lock lock(m_system_log_mutex);
  m_logs.push_back( log );
}

void
vw::Log::clear() {
  Mutex::Lock lock(m_system_log_mutex);
  m_logs.clear();
}

vw::LogInstance&
vw::Log::console_log() {
  Mutex::Lock lock(m_system_log_mutex);
  return *m_console_log;
}

void
vw::Log::set_console_stream(std::ostream& stream, LogRuleSet rule_set, bool prepend_infostamp ) {
  Mutex::Lock lock(m_system_log_mutex);
  m_console_log = boost::shared_ptr<LogInstance>(new LogInstance(stream, prepend_infostamp) );
  m_console_log->rule_set() = rule_set;
}

bool vw::Log::is_enabled( int log_level,
                          std::string const& log_namespace ) {
  // Early exit option before iterating through m_logs
  if ( m_console_log->rule_set()(log_level, log_namespace) )
    return true;
  BOOST_FOREACH( boost::shared_ptr<LogInstance> log_instance, m_logs ) {
    if ( log_instance->rule_set()(log_level, log_namespace) )
      return true;
  }
  return false;
}

vw::LogRuleSet::LogRuleSet( LogRuleSet const& copy_log) {
  m_rules = copy_log.m_rules;
}

vw::LogRuleSet& vw::LogRuleSet::operator=( LogRuleSet const& copy_log) {
  m_rules = copy_log.m_rules;
  return *this;
}


vw::LogRuleSet::LogRuleSet() { }
vw::LogRuleSet::~LogRuleSet() { }


// TODO: Clean up so valgrind does not report a memory leak
void vw::LogRuleSet::add_rule(int log_level, std::string const& log_namespace) {
  ssize_t count = std::count(log_namespace.begin(), log_namespace.end(), '*');
  if (count > 1)
    vw::vw_throw(vw::ArgumentErr() << "Illegal log rule: only one wildcard is supported.");

  if (count == 1 && *(log_namespace.begin()) != '*'
      && *(log_namespace.end()-1) != '*')
    vw::vw_throw(vw::ArgumentErr() << "Illegal log rule: wildcards must be at the beginning or end of a rule");

  Mutex::Lock lock(m_mutex);
  m_rules.push_front(rule_type(log_level, boost::to_lower_copy(log_namespace)));
}

void vw::LogRuleSet::clear() {
  Mutex::Lock lock(m_mutex);
  m_rules.clear();
}

namespace {
  bool wildcard_match(const std::string& pattern, const std::string& str) {
    // Rules:
    // *   matches anything
    // *.a matches [first.a, second.a]
    // a   matches just a (ie, no namespace)
    // a.* matches [a, a.first, a.first.second]
    //

    if (pattern == "*")
      return true;

    // If there's no wildcard, just do a comparison.
    size_t idx = pattern.find("*");
    if (idx == std::string::npos)
      return (pattern == str);

    // There's a wildcard. Try to expand it.
    if (idx == 0) {
      // leading *. it's a suffix rule.
      return boost::ends_with(str, pattern.substr(1));
    } else {
      // add_rule above verifies that the wildcard is first or last, so this
      // one must be last.
      if (pattern.size() > 1 && pattern[idx-1] == '.')
        if (str == pattern.substr(0, idx-1))
          return true;
      return boost::starts_with(str, pattern.substr(0, idx));
    }

    return false;
  }
}

// You can overload this method from a subclass to change the
// behavior of the LogRuleSet.
bool vw::LogRuleSet::operator() (int log_level, std::string const& log_namespace) {
  Mutex::ReadLock lock(m_mutex);

  std::string lower_namespace = boost::to_lower_copy(log_namespace);

  for (rules_type::const_iterator it = m_rules.begin(); it != m_rules.end(); ++it) {
    const int&         rule_lvl = it->first;
    const std::string& rule_ns  = it->second;

    // first rule that matches the namespace spec is applied.

    if (!wildcard_match(rule_ns, lower_namespace) )
      continue;

    if (rule_lvl == vw::EveryMessage)
      return true;

    return log_level <= rule_lvl;
  }

  if (log_level <= vw::InfoMessage)
    if (log_namespace == "console" || wildcard_match("*.progress", lower_namespace))
      return true;
  if (log_level <= vw::WarningMessage)
    return true;

  // We reach this line if all of the rules have failed, in
  // which case we return a NULL stream, which will result in
  // nothing being logged.
  return false;
}
