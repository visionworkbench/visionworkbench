// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <iostream>
#include <fstream>
#include <sstream>

#include <gtest/gtest.h>
#include <vw/Core/Settings.h>
#include <vw/Core/ConfigParser.h>
#include <test/Helpers.h>

using namespace vw;
using namespace vw::test;

TEST(Settings, HAS_CONFIG_FILE(VWrc)) {

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

  UnlinkName file("test_vwrc");
  std::ofstream ostr(file.c_str());
  ASSERT_TRUE(ostr.is_open()) << "Could not open test config file for writing";
  ostr << conf;
  ostr.close();

  // Test to see if the settings were correctly read in
  vw_settings().set_rc_filename(file);
  EXPECT_EQ( 20u, vw_settings().default_num_threads() );
  EXPECT_EQ( 623u, vw_settings().system_cache_size() );

  // Test to make sure that the API overrides the contents of vwrc
  vw_settings().set_default_num_threads(5);
  vw_settings().set_system_cache_size(223);
  EXPECT_EQ( 5u, vw_settings().default_num_threads() );
  EXPECT_EQ( 223u, vw_settings().system_cache_size() );
}

TEST(Settings, Override) {
  // Test to make sure that the API overrides the contents of vwrc
  vw_settings().set_default_num_threads(5);
  vw_settings().set_system_cache_size(223);
  EXPECT_EQ( 5u, vw_settings().default_num_threads() );
  EXPECT_EQ( 223u, vw_settings().system_cache_size() );
}

TEST(SettingsDeathTest, OldVWrc) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  const char *conf = "\n\
                      logfile console\n\
                      40 *\n";

  Settings s;
  std::istringstream stream(conf);

  EXPECT_EXIT(parse_config(stream, s); exit(251);, ::testing::ExitedWithCode(251), "Could not parse config file");
}
