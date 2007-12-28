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
  boost::shared_ptr<vw::SystemLog> system_log_ptr;
  void init_system_log() {
    system_log_ptr = boost::shared_ptr<vw::SystemLog>(new vw::SystemLog());
  }
}

// ---------------------------------------------------
// Basic stream support
// ---------------------------------------------------
std::ostream& vw::vw_out( int log_level, std::string log_namespace ) {
  return system_log()(log_level, log_namespace);
}

void vw::set_debug_level( int log_level ) {
  system_log().console_log().rule_set().add_rule("console", log_level);
}

void vw::set_output_stream( std::ostream& stream ) {
  system_log().set_console_stream(stream);
}

vw::SystemLog& vw::system_log() {
  system_log_once.run( init_system_log );
  return *system_log_ptr;
}

// ---------------------------------------------------
// Log Methods
// ---------------------------------------------------
std::ostream& vw::Log::operator() (int level, std::string log_namespace) {
  if (m_rule_set(level, log_namespace)) {
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
// SystemLog Methods
// ---------------------------------------------------
void vw::SystemLog::reload_logconf_rules() {
  std::ifstream f(m_logconf_filename.c_str());
  
  if (f.is_open()) {
    while (!f.eof()) {
      char c_line[2048];
      f.getline(c_line, 2048);
      std::string line = boost::trim_copy(std::string(c_line));
      std::vector<std::string> split_vec;
      boost::split(split_vec, line, boost::is_any_of(" "));
      if (split_vec.size() != 2) {
        // failed
      }
      int log_level = atoi(split_vec[1].c_str());
      //      m_file_log_ruleset.add_rule(split_vec[0], log_level);
    }
  }
}

// Every m_log_settings_period seconds, this method polls the
// m_log_settings_file to see if it exists and to see if it has been
// recently modified.  If so, we reload the log ruleset from the file.
void vw::SystemLog::stat_logconf() {
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
    if ( f = fopen(m_logconf_filename.c_str(), "r") ) {
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

std::ostream& vw::SystemLog::operator() (int level, std::string log_namespace) { 
  // First, check to see if the logconf file has been updated.
  // Reload the rulesets if it has.
  stat_logconf();
  
  // Add the console log output...
  m_multi_ostream.clear();
  m_multi_ostream.add(m_console_log->operator()(level, log_namespace));
  
  // ... and the rest of the active log streams.
  std::vector<boost::shared_ptr<Log> >::iterator iter = m_logs.begin();
  for (;iter != m_logs.end(); ++iter) 
    m_multi_ostream.add((*iter)->operator()(level,log_namespace));

  return m_multi_ostream;
}


