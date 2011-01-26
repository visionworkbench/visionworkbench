// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
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

#include <vector>

#include <vw/Core/System.h>
#include <vw/Core/Thread.h>

namespace vw {

  // -------------------------------------------------------
  //                    Settings
  // -------------------------------------------------------

  /// The system settings class manages runtime configuration for the
  /// Vision Workbench.
  ///
  /// Important Note: You should access the system settings using the
  /// static vw_settings() function below, which returns a singleton
  /// instance of the system settings class.  You should _not_ need to
  /// create a settings object yourself!!!!
  class Settings : private boost::noncopyable {

#define VW_DECLARE_SETTING(Name, Type)\
    private:\
      Type m_ ## Name; \
      bool m_ ## Name ## _override; \
    public: \
      Type Name(); \
      void set_ ## Name(const Type& x);

    // The default number of threads used in block processing operations.
    VW_DECLARE_SETTING(default_num_threads, uint32);

    // The current system cache size (in bytes). The system cache is shared by
    // all BlockRasterizeView<>'s, including DiskImageView<>'s.
    VW_DECLARE_SETTING(system_cache_size, size_t);

    // Write cache is only used in block writing. This is the number of threads
    // that can be blocked on IO before the code stops creating more jobs (to
    // let the writes catch up).
    VW_DECLARE_SETTING(write_pool_size, uint32);

    // The default tile size (in pixels) used for block processing ops.
    VW_DECLARE_SETTING(default_tile_size, uint32);

    // The directory used to store temporary files.
    VW_DECLARE_SETTING(tmp_directory, std::string);

#undef VW_DECLARE_SETTING

    // Member variables assoc. with periodically polling the log
    // configuration (logconf) file.
    long m_rc_last_polltime;
    long m_rc_last_modification;
    std::string m_rc_filename;
    float m_rc_poll_period;
    RecursiveMutex m_rc_time_mutex;
    RecursiveMutex m_rc_file_mutex;
    RecursiveMutex m_settings_mutex;

  public:

    /// You should not create an instance of Settings on your own
    /// using this constructor.  Instead, you can access a global
    /// instance of the settings class using the static
    /// vw_settings() method below.
    Settings();

    /// Change the rc filename (default: ~/.vwrc)
    void set_rc_filename(std::string filename, bool parse_now = true);

    /// Change the rc file poll period.  (default: 5 seconds)
    /// Note -- this sets the *minimum* poll time for the file.  The
    /// actual file is only polled when a setting is requested.
    void set_rc_poll_period(float period);

    // -----------------------------------------------------------------
    //                        Settings API
    // -----------------------------------------------------------------

    // Reload the config. When this returns, the config will be reloaded.
    // It may happen in a different thread than requested, however.
    void reload_config();
  };

} // namespace vw

#endif // __VW_CORE_SETTINGS_H__
