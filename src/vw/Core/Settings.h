// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Core/Settings.h
/// 
/// This file provides a singleton object that provides access to
/// Vision Workbench system-wide settings.  These can be twiddled
/// programmatically by interacting with this object, or they can be
/// set using the ~/.vwrc file in the user's home directory.  That
/// file is reloaded on every call to access settings, so the user can
/// modify the contents of that file in real time as the program is
/// running.

#ifndef __VW_CORE_SETTINGS_H__
#define __VW_CORE_SETTINGS_H__

#include <vw/Core/Thread.h>

namespace vw {

  // -------------------------------------------------------
  //                    Settings
  // -------------------------------------------------------

  /// The system settings class manages runtime configuration for the
  /// Vision Workbench.
  ///
  /// Important Note: You should access the system settings using the
  /// Settings::system_settings() static method, which access a
  /// singleton instance of the system settings class.  You should not
  /// need to create a settings object yourself!!
  class Settings {

    struct LogSetting {
      std::string filename;
      std::string rules;
    };

    // Vision Workbench Global Settings
    int m_default_num_threads;
    std::vector< LogSetting > m_log_settings;

    bool m_default_num_threads_override;
    
    // Member variables assoc. with periodically polling the log
    // configuration (logconf) file.
    long m_vwrc_last_polltime;
    long m_vwrc_last_modification;
    std::string m_vwrc_filename;
    double m_vwrc_poll_period;
    Mutex m_vwrc_time_mutex;
    Mutex m_vwrc_file_mutex;
    Mutex m_settings_mutex;
    
    // Private methods for checking the logconf file
    void stat_vwrc();
    void reload_vwrc();

    // Ensure non-copyable semantics
    Settings( Settings const& );
    Settings& operator=( Settings const& );

  public:

    /// You should not create an instance of Settings on your own
    /// using this constructor.  Instead, you can access a global
    /// instance of the settings class using the static
    /// Setting::system_settings() method below.
    Settings();   

    /// Change the vwrc filename (default: ~/.vwrc)
    void set_vwrc_filename(std::string filename) { 
      m_vwrc_filename = filename; 
      m_vwrc_last_polltime = 0;
      stat_vwrc();
    }

    /// Change the vwrc file poll period.  (default: 5 seconds)
    /// Note -- this sets the *minimum* poll time for the file.  The
    /// actual file is only polled when a setting is requested.
    void set_vwrc_poll_period(double period) { 
      m_vwrc_poll_period = period; 
      m_vwrc_last_polltime = 0;
      stat_vwrc();
    }

    // -----------------------------------------------------------------
    //                        Settings API
    // -----------------------------------------------------------------    

    int default_num_threads();
    void set_default_num_threads(int num);

    std::vector<LogSetting> log_settings();
    // Note: you can set the log settings programmatically using the
    // Log::system_log() object in Log.h.  It will ignore the settings
    // in .vwrc if you have done so.
  };
  
  // Static method to access the singleton instance of the system
  // log.  You should *always* use this method if you want to access
  // Vision Workbench system log, where all Vision Workbench log
  // messages go.
  Settings& vw_settings();

} // namespace vw

#endif // __VW_CORE_SETTINGS_H__
