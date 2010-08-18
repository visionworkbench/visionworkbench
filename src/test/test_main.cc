// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>
#include <test/Helpers.h>
#include <vw/Core/Settings.h>

#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;

int main(int argc, char **argv) {
  // Disable the user's config file
  vw::vw_settings().set_rc_filename("");
  ::testing::InitGoogleTest(&argc, argv);

  // Default to the "threadsafe" style because we can't delete our singletons
  // yet; this style of test launches a new process, so the singletons are
  // fresh.
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  return RUN_ALL_TESTS();
}

vw::test::UnlinkName::UnlinkName(std::string base, std::string directory)
  : std::string(directory + "/" + base) {
  fs::remove_all(this->c_str());
}

vw::test::UnlinkName::UnlinkName(const char* base, std::string directory)
  : std::string(directory + "/" + base) {
  fs::remove_all(this->c_str());
}

vw::test::UnlinkName::~UnlinkName() {
  if (!this->empty())
    fs::remove_all(this->c_str());
}
