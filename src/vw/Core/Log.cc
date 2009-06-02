// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vector>
#include <deque>

#include <vw/Core/Log.h>
#include <vw/Core/Exception.h>
#include <vw/Core/Settings.h>

// Boost headers
#include <boost/thread/xtime.hpp>

// C Standard Library headers ( for stat(2) and getpwuid() )
#include <sys/types.h>
#include <sys/stat.h>
#include <ctime>

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

// ---------------------------------------------------
// Create a single instance of the SystemLog
// ---------------------------------------------------
namespace {
  static vw::null_ostream g_null_ostream;
  vw::RunOnce system_log_once = VW_RUNONCE_INIT;
  boost::shared_ptr<vw::Log> system_log_ptr;
  void init_system_log() {
    system_log_ptr = boost::shared_ptr<vw::Log>(new vw::Log());
  }
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
      return m_log_stream << current_posix_time_string() << " {" << Thread::id() << "} [ " << log_namespace << " ] : ";
    else
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

    // Reset and add the console log output...
    m_multi_ostreams[ Thread::id() ]->clear();
    m_multi_ostreams[ Thread::id() ]->add(m_console_log->operator()(log_level, log_namespace));

    // ... and the rest of the active log streams.
    std::vector<boost::shared_ptr<LogInstance> >::iterator iter = m_logs.begin();
    for (;iter != m_logs.end(); ++iter)
      m_multi_ostreams[ Thread::id() ]->add((*iter)->operator()(log_level,log_namespace));

    return *m_multi_ostreams[ Thread::id() ];
  }
}

vw::Log& vw::vw_log() {
  system_log_once.run( init_system_log );
  return *system_log_ptr;
}
