// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Core/TemporaryFile.h>
#include <vw/Core/Settings.h>
#include <gtest/gtest.h>
#include <test/Helpers.h>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <cerrno>

namespace fs = boost::filesystem;

#ifdef VW_HAVE_UNISTD_H
#  include <unistd.h>
#endif

#ifdef WIN32
#  define access _access
#endif

using namespace std;
using namespace vw;

string get_directory(const TemporaryFile& f) {
  return fs::path(f.filename()).parent_path().string();
}

string get_prefix(const TemporaryFile& f, size_t len) {
  return fs::path(f.filename()).leaf().substr(0, len);
}

string get_suffix(const TemporaryFile& f, size_t len) {
  string fn = fs::path(f.filename()).leaf();
  return fn.substr(fn.size()-len, len);
}

#define check1(a)           _check(a,__LINE__)
#define check2(a,b)         _check(a,__LINE__,b)
#define check3(a,b,c)       _check(a,__LINE__,b,c)
#define check4(a,b,c,d)     _check(a,__LINE__,b,c,d)
#define check5(a,b,c,d,e)   _check(a,__LINE__,b,c,d,e)

void _check(TemporaryFile& a, uint32 line, std::string dir = "", const std::string& prefix = "tmp", const std::string& suffix = "", bool can_read = true) {
  SCOPED_TRACE(::testing::Message() << "line " << line);

  if (dir.empty())
    dir = vw_settings().tmp_directory();

  const std::string& fn = a.filename();
  EXPECT_TRUE(fs::equivalent(dir, get_directory(a))) << "Expected " << dir << ", got " << get_directory(a);
  EXPECT_EQ(suffix, get_suffix(a, suffix.size()));
  EXPECT_EQ(prefix, get_prefix(a, prefix.size()));
  EXPECT_TRUE(fs::exists(fn));

  const char data[5] = "rawr";
  char buf[5] = {0};

  ASSERT_FALSE(a.bad());
  a.write(data, 5);
  ASSERT_FALSE(a.fail());

  a.seekg(0);
  ASSERT_TRUE(a.good());
  a.read(buf, 5);

  if (!can_read)
    ASSERT_TRUE(a.fail());
  else {
    ASSERT_FALSE(a.fail());
    EXPECT_RANGE_EQ(data+0, data+5, buf+0, buf+5);
  }
}

TEST(TemporaryFile, NoArg) {
  std::string fn;
  {
    TemporaryFile a;
    fn = a.filename();
    check1(a);
  }
  EXPECT_FALSE(fs::exists(fn));
}

TEST(TemporaryFile, NoDelete) {
  string fn;
  {
    TemporaryFile a(TEST_OBJDIR, false);
    fn = a.filename();
    check2(a, TEST_OBJDIR);
  }
  EXPECT_TRUE(fs::exists(fn));
  ::remove(fn.c_str());
}

TEST(TemporaryFile, Basic) {
    TemporaryFile a;
    check1(a);

    TemporaryFile b(TEST_OBJDIR);
    check2(b, TEST_OBJDIR);
    TemporaryFile c(TEST_OBJDIR, true);
    check2(c, TEST_OBJDIR);
    TemporaryFile d(TEST_OBJDIR, true, "prefix");
    check3(d, TEST_OBJDIR, "prefix");
    TemporaryFile e(TEST_OBJDIR, true, "prefix", "suffix");
    check4(e, TEST_OBJDIR, "prefix", "suffix");
    TemporaryFile f(TEST_OBJDIR, true, "prefix", "suffix", std::ios::out);
    check5(f, TEST_OBJDIR, "prefix", "suffix", false);
}
