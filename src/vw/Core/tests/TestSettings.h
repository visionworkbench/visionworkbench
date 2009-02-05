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

  void test_system_settings() {
    vw_out(0) << "\nTesting System Settings\n";

    vw_out(0) << "\t--> DEFAULT_NUM_THREADS = " << vw_settings().default_num_threads() << "\n";
    vw_settings().set_default_num_threads(3);
    TS_ASSERT_EQUALS(vw_settings().default_num_threads(), 3);

  }


  void test_vwrc() {
    std::ofstream ostr("./test_vwrc");
    ostr << "# Comment 1\n";
    ostr << "# Comment 2\n";
    ostr << "\n";
    ostr << "NONEXISTENT_ENTRY 1\n";
    ostr << "\n";
    ostr << "DEFAULT_NUM_THREADS 20\n";
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

    // Test to make sure that the API overrides the contents of vwrc
    vw_settings().set_default_num_threads(5);
    TS_ASSERT_EQUALS(vw_settings().default_num_threads(), 5);
  }

};
