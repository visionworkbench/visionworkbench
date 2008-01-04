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
#include <vw/Core/Log.h>

// Boost header
#include <boost/thread/xtime.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

// C Standard Library headers ( for stat(2) )
#include <sys/types.h>
#include <sys/stat.h>

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
  return Log::system_log()(log_level, log_namespace);
}

void vw::set_debug_level( int log_level ) {
  Log::system_log().console_log().rule_set().add_rule(log_level, "console");
}

void vw::set_output_stream( std::ostream& stream ) {
  Log::system_log().set_console_stream(stream);
}

// ---------------------------------------------------
// LogInstance Methods
// ---------------------------------------------------
vw::LogInstance::LogInstance(std::string log_filename, bool prepend_infostamp) : m_prepend_infostamp(prepend_infostamp) {
  // Open file and place the insertion pointer at the end of the file (ios_base::ate)
  m_log_ostream_ptr = new std::ofstream(log_filename.c_str(), std::ios::app);
  if (! static_cast<std::ofstream*>(m_log_ostream_ptr)->is_open())
    vw_throw(IOErr() << "Could not open log file " << log_filename << " for writing.");
  
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  *m_log_ostream_ptr << "\n\n" << "Vision Workbench log started at " << now << ".\n\n";
      
  m_log_stream.set_stream(*m_log_ostream_ptr);
}

vw::LogInstance::LogInstance(std::ostream& log_ostream, bool prepend_infostamp) : m_log_stream(log_ostream), 
                                                                                  m_log_ostream_ptr(NULL),
                                                                                  m_prepend_infostamp(prepend_infostamp) {}

std::ostream& vw::LogInstance::operator() (int log_level, std::string log_namespace) {
  if (m_rule_set(log_level, log_namespace)) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    if (m_prepend_infostamp) 
      return m_log_stream << now << " {" << Thread::id() << "} [ " << log_namespace << " ] : ";
    else 
      return m_log_stream;
  } else {
    return g_null_ostream;
  }
}

// ---------------------------------------------------
// Log Methods
// ---------------------------------------------------
void vw::Log::reload_logconf_rules() {
  vw_out(InfoMessage, "log") << "Reloading log configuration file: " << m_logconf_filename << ".\n";
  
  std::ifstream f(m_logconf_filename.c_str());
  
  if (f.is_open()) {
    boost::shared_ptr<LogInstance> current_log;
    while (!f.eof()) {
      char c_line[2048];
      f.getline(c_line, 2048);
      std::string line = boost::trim_copy(std::string(c_line));

      // Check to see if this line is empty or if it starts with '#',
      // which we ignore as a comment.
      if ( line.size() != 0 && line[0] != '#' ) {
        std::vector<std::string> tokens;
        boost::split(tokens, line, boost::is_any_of(" "));

        // All lines in the file should contain exactly two tokens
        // seperated by a space.  If not, we ignore the line.
        if (tokens.size() == 2) {

          // Handle the wildcard "*" for log level.
          int log_level;
          if (tokens[0] == "*")
            log_level = vw::EveryMessage;
          else 
            log_level = atoi(tokens[0].c_str());            

          if (boost::to_lower_copy(tokens[0]) == "logfile") {
            // If the first token is the "logfile" string, we start a
            // new log file stream with the supplied filename.
            if (tokens[1] == "console") {
              vw_out(DebugMessage, "log") << "Adding rules for console.\n";
              current_log.reset();
            } else {
              vw_out(DebugMessage, "log") << "Adding rules for log file: " << tokens[1] << ".\n";
              current_log = boost::shared_ptr<LogInstance>( new LogInstance(tokens[1]) );
              this->add(current_log);
            }        
          } else {
            // Otherwise, we need to add the rule for the current log.
            vw_out(DebugMessage, "log") << "Adding rule: " << tokens[0] << "   " << tokens[1] << "\n";
            if (current_log) 
              current_log->rule_set().add_rule(log_level, tokens[1]);
            else 
              m_console_log->rule_set().add_rule(log_level, tokens[1]);
          }
        }
      }
    }
  }
}

// Every m_log_settings_period seconds, this method polls the
// m_log_settings_file to see if it exists and to see if it has been
// recently modified.  If so, we reload the log ruleset from the file.
void vw::Log::stat_logconf() {
  boost::xtime xt;
  boost::xtime_get(&xt, boost::TIME_UTC);
  bool needs_reloading = false;

  // Every five seconds, we attempt to open the log config file to see
  // if there have been any changes.  The mutex locking for querying
  // the time is handled seperately from reading the file so that only
  // one thread takes the performance hit of reading the logconf file
  // during any given reload.
  {
    Mutex::Lock time_lock(m_logconf_time_mutex);
    if (xt.sec - m_logconf_last_polltime > m_logconf_poll_period) {
      m_logconf_last_polltime = xt.sec;
      needs_reloading = true;
    }
  }
  
  if (needs_reloading) {
    Mutex::Lock lock(m_logconf_file_mutex);
    FILE *f;
    if ( (f = fopen(m_logconf_filename.c_str(), "r")) ) {
      fclose(f);

      // Check to see if the file has changed.  If so, re-read the
      // settings.
      struct stat stat_struct;
      if (stat(m_logconf_filename.c_str(), &stat_struct) == 0) {
        if (stat_struct.st_mtimespec.tv_sec > m_logconf_last_modification) {
          m_logconf_last_modification = stat_struct.st_mtimespec.tv_sec;
          reload_logconf_rules();
        }
      }
    }
  }
}

std::ostream& vw::Log::operator() (int log_level, std::string log_namespace) { 
  // First, check to see if the logconf file has been updated.
  // Reload the rulesets if it has.
  stat_logconf();

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

vw::Log& vw::Log::system_log() {
  system_log_once.run( init_system_log );
  return *system_log_ptr;
}
