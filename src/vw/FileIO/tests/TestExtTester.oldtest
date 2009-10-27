// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestBlockFileIO.h
#define CXXTEST_ABORT_TEST_ON_FAIL
#include <cxxtest/TestSuite.h>
#include <string>
#include <vw/Core/Log.h>
#include <vw/FileIO/DiskImageResource_internal.h>

using namespace vw;
using namespace vw::internal;

static void test_extension(std::string const& fn)
{
  TS_TRACE(vw::stringify("Testing ") + fn);
}

class TestExtTester : public CxxTest::TestSuite
{
public:

  void test_tester()
  {
    foreach_ext("rwtest", test_extension);
  }
}; // class TestBlockFileIO
