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

#include <vw/config.h>
#include <vw/Core/Thread.h>
#include <vw/Core/Settings.h>

// Boost headers
#include <boost/algorithm/string.hpp>
#include <boost/thread/xtime.hpp>
// Posix time is not fully supported in the version of Boost for RHEL
// Workstation 4
#ifdef __APPLE__
#include <boost/date_time/posix_time/posix_time.hpp>
#else
#include <ctime>
#endif

// C Standard Library headers ( for stat(2) and getpwuid() )
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <pwd.h>
#include <fstream>

// ---------------------------------------------------
// Create a single instance of the Settings class
// ---------------------------------------------------
namespace {
  vw::RunOnce settings_once = VW_RUNONCE_INIT;
  boost::shared_ptr<vw::Settings> system_settings_ptr;
  void init_system_settings() {
    system_settings_ptr = boost::shared_ptr<vw::Settings>(new vw::Settings());
  }
}

// ---------------------------------------------------
// Settings Methods
// ---------------------------------------------------
void vw::Settings::reload_vwrc() {
  std::ifstream f(m_vwrc_filename.c_str());

  if (f.is_open()) {
    m_log_settings.clear();
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
          if (boost::to_upper_copy(tokens[0]) == "DEFAULT_NUM_THREADS") {
            m_default_num_threads = atoi(tokens[1].c_str());

          } else if (boost::to_upper_copy(tokens[0]) == "LOGFILE") {
            LogSetting l;
            l.filename = tokens[1];

            // LOGFILE <name> is following by entries contained in
            // {}'s .  The closing bracket must be at the beginning of
            // a line of its own.
            while (line[0] != '}') {
              f.getline(c_line, 2048);
              line = std::string(c_line) + "\n";
              if (line[0] != '{' && line[0] != '}') {
                l.rules += line;
              }
            }
            m_log_settings.push_back(l);
          }
        }
      }
    }
  }
}

// Every m_vwrc_poll_period seconds, this method polls the
// m_vwrc_filename to see if it exists and to see if it has been
// recently modified.  If so, we reload the log ruleset from the file.
void vw::Settings::stat_vwrc() {
  boost::xtime xt;
  boost::xtime_get(&xt, boost::TIME_UTC);
  bool needs_reloading = false;

  // Every five seconds, we attempt to open the log config file to see
  // if there have been any changes.  The mutex locking for querying
  // the time is handled seperately from reading the file so that only
  // one thread takes the performance hit of reading the vwrc file
  // during any given reload.
  {
    Mutex::Lock time_lock(m_vwrc_time_mutex);
    if (xt.sec - m_vwrc_last_polltime > m_vwrc_poll_period) {
      m_vwrc_last_polltime = xt.sec;
      needs_reloading = true;
    }
  }
  
  if (needs_reloading) {
    Mutex::Lock lock(m_vwrc_file_mutex);
    FILE *f;
    if ( (f = fopen(m_vwrc_filename.c_str(), "r")) ) {
      fclose(f);

      // Check to see if the file has changed.  If so, re-read the
      // settings.
      struct stat stat_struct;
      if (stat(m_vwrc_filename.c_str(), &stat_struct) == 0) {
#ifdef __APPLE__
        if (stat_struct.st_mtimespec.tv_sec > m_vwrc_last_modification) {
          m_vwrc_last_modification = stat_struct.st_mtimespec.tv_sec;
          reload_vwrc();
        }
#else // Linux
        if (stat_struct.st_mtime > m_vwrc_last_modification) {
          m_vwrc_last_modification = stat_struct.st_mtime;
          reload_vwrc();
        }
#endif
      }
    }
  }
}


vw::Settings::Settings() : m_vwrc_last_polltime(0),
                           m_vwrc_last_modification(0),
                           m_vwrc_poll_period(5.0) {
  struct passwd *pw;
  pw = getpwuid( getuid() );
  std::string homedir = pw->pw_dir; 
  m_vwrc_filename = homedir + "/.vwrc";

  // Set defaults
  m_default_num_threads = VW_NUM_THREADS;

  // By default, the .vwrc file has precedence, but the user can
  // override these settings by explicitly changing them using the
  // system_settings() API.
  m_default_num_threads_override = false;
}

vw::Settings& vw::vw_settings() {
  settings_once.run( init_system_settings );
  return *system_settings_ptr;
}


// -----------------------------------------------------------------
//                        Settings API
// -----------------------------------------------------------------    

int vw::Settings::default_num_threads() { 
  if (!m_default_num_threads_override)
    stat_vwrc(); 
  Mutex::Lock lock(m_settings_mutex);
  return m_default_num_threads; 
}

void vw::Settings::set_default_num_threads(int num) { 
  Mutex::Lock lock(m_settings_mutex);
  m_default_num_threads_override = true;
  m_default_num_threads = num; 
}


std::vector<vw::Settings::LogSetting> vw::Settings::log_settings() { 
  stat_vwrc(); 
  Mutex::Lock lock(m_settings_mutex);
  return m_log_settings;
}
