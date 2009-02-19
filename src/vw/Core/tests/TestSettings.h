// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestThreadPool.h
#include <cxxtest/TestSuite.h>

#include <vw/Core/Settings.h>
#include <vw/Core/Log.h>

#include <fstream>
#include <iostream>

using namespace vw;

class TestLog : public CxxTest::TestSuite
{
public:

  void test_vwrc() {
    std::ofstream ostr("./test_vwrc");
    ostr << "# Comment 1\n";
    ostr << "# Comment 2\n";
    ostr << "\n";
    ostr << "NONEXISTENT_ENTRY 1\n";
    ostr << "\n";
    ostr << "DEFAULT_NUM_THREADS 20\n";
    ostr << "SYSTEM_CACHE_SIZE 623\n";
    ostr << "# Comment \n";
    ostr << "\n";
    ostr << "LOGFILE console\n";
    ostr << "{\n";
    ostr << "  * any\n";    
    ostr << "  fileio 10\n";    
    ostr << "}\n";
    ostr.close();

    // Test to see if the settings were correctly read in
    vw_settings().set_vwrc_filename("./test_vwrc");
    TS_ASSERT_EQUALS(vw_settings().default_num_threads(), 20);
    TS_ASSERT_EQUALS(vw_settings().system_cache_size(), 623);

    // Test to make sure that the API overrides the contents of vwrc
    vw_settings().set_default_num_threads(5);
    vw_settings().set_system_cache_size(223);
    TS_ASSERT_EQUALS(vw_settings().default_num_threads(), 5);
    TS_ASSERT_EQUALS(vw_settings().system_cache_size(), 223);
  }

};
