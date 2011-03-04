// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vector>

#include <vw/Core/Log.h>
#include <vw/Core/Exception.h>
#include <vw/Core/Settings.h>

// C Standard Library headers ( for stat(2) and getpwuid() )
#include <sys/types.h>
#include <sys/stat.h>
#include <ctime>

#include <vw/config.h>

#ifdef VW_HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef VW_HAVE_PWD_H
#include <pwd.h>
#endif

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

namespace {
  static vw::null_ostream g_null_ostream;
}

// ---------------------------------------------------
// Basic stream support
// ---------------------------------------------------
std::ostream& vw::vw_out( int log_level, std::string log_namespace ) {
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
vw::LogInstance::LogInstance(std::string log_filename, bool prepend_infostamp) : m_prepend_infostamp(prepend_infostamp) {
  // Open file and place the insertion pointer at the end of the file (ios_base::ate)
  m_log_ostream_ptr = new std::ofstream(log_filename.c_str(), std::ios::app);
  if (! static_cast<std::ofstream*>(m_log_ostream_ptr)->is_open())
    vw_throw(IOErr() << "Could not open log file " << log_filename << " for writing.");

  *m_log_ostream_ptr << "\n\n" << "Vision Workbench log started at " << current_posix_time_string() << ".\n\n";

  m_log_stream.set_stream(*m_log_ostream_ptr);
}

vw::LogInstance::LogInstance(std::ostream& log_ostream, bool prepend_infostamp) : m_log_stream(log_ostream),
                                                                                  m_log_ostream_ptr(NULL),
                                                                                  m_prepend_infostamp(prepend_infostamp) {}

std::ostream& vw::LogInstance::operator() (int log_level, std::string log_namespace) {
  if (m_rule_set(log_level, log_namespace)) {
    if (m_prepend_infostamp)
      m_log_stream << current_posix_time_string() << " {" << Thread::id() << "} [ " << log_namespace << " ] : ";
    switch (log_level) {
      case ErrorMessage:   m_log_stream << "Error: ";   break;
      case WarningMessage: m_log_stream << "Warning: "; break;
      default: break;
    }
    return m_log_stream;
  } else {
    return g_null_ostream;
  }
}


std::ostream& vw::Log::operator() (int log_level, std::string log_namespace) {
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

vw::LogRuleSet::LogRuleSet( LogRuleSet const& copy_log) {
  m_rules = copy_log.m_rules;
}

vw::LogRuleSet& vw::LogRuleSet::operator=( LogRuleSet const& copy_log) {
  m_rules = copy_log.m_rules;
  return *this;
}


vw::LogRuleSet::LogRuleSet() { }
vw::LogRuleSet::~LogRuleSet() { }

void vw::LogRuleSet::add_rule(int log_level, std::string log_namespace) {
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
bool vw::LogRuleSet::operator() (int log_level, std::string log_namespace) {
  Mutex::Lock lock(m_mutex);

  std::string lower_namespace = boost::to_lower_copy(log_namespace);

  for (rules_type::iterator it = m_rules.begin(); it != m_rules.end(); ++it) {
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
