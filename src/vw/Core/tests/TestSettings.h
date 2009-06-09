// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <iostream>
#include <fstream>
#include <sstream>

#include <cxxtest/TestSuite.h>
#include <vw/Core/Settings.h>
#include <vw/Core/ConfigParser.h>
#include <vw/config.h>

using namespace vw;

class TestSettings : public CxxTest::TestSuite
{
public:


  void test_vwrc() {
#if defined(VW_ENABLE_CONFIG_FILE) && VW_ENABLE_CONFIG_FILE == 1
    const char *conf = "\n\
      # Comment 1                         \n\
      # Comment 2                         \n\
                                          \n\
      [general]                           \n\
                                          \n\
      nonexistent_entry = 1               \n\
      default_num_threads = 20            \n\
      system_cache_size = 623             \n\
      # Comment                           \n\
                                          \n\
      [logfile console]                   \n\
      10 = *           # all level 10     \n\
      30 = WMS.foo     # WMS.foo level 30 \n\
                                          \n\
      [logfile log.txt]                   \n\
      *  = * wee                          \n\
      10 = WMS.foo                        \n\
                                          \n\
      [logfile console]                   \n\
      VerboseDebugMessage = *             \n\
      ";

    std::ofstream ostr(TEST_SRCDIR"/test_vwrc");
    ostr << conf;
    ostr.close();

    // Test to see if the settings were correctly read in
    vw_settings().set_rc_filename(TEST_SRCDIR"/test_vwrc");
    TS_ASSERT_EQUALS(vw_settings().default_num_threads(), 20);
    TS_ASSERT_EQUALS(vw_settings().system_cache_size(), 623);
#else
    TS_TRACE("Config files disabled: skipping read test!");
#endif

    // Test to make sure that the API overrides the contents of vwrc
    vw_settings().set_default_num_threads(5);
    vw_settings().set_system_cache_size(223);
    TS_ASSERT_EQUALS(vw_settings().default_num_threads(), 5);
    TS_ASSERT_EQUALS(vw_settings().system_cache_size(), 223);
  }

  void test_old_vwrc() {
    const char *conf = "\n\
logfile console\n\
40 *\n";

    Settings s;
    std::istringstream stream(conf);

    TS_TRACE("This should print a warning:");
    parse_config(stream, s);
  }
};
