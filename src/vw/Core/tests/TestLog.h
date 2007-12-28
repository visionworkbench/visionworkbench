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

#include <vw/Core/Log.h>
#include <vw/Core/Thread.h>

#include <iostream>
#include <string>

using namespace vw;

class TestLog : public CxxTest::TestSuite
{

  struct TestLogTask {
    bool m_terminate;
    Log& m_log;
    std::string m_ns;
    TestLogTask(Log &log, std::string ns) : m_terminate(false), m_log(log), m_ns(ns) {}

    void operator()() {
      m_log(0,m_ns) << "\tThread " << Thread::id() << " activated\n";

      while( !m_terminate ) {
        Thread::sleep_ms(100);
        m_log(0,m_ns) << "\tThread " << Thread::id() << " logging a test message\n";
      }

      m_log(0,m_ns) << "\tThread " << Thread::id() << " terminated\n";
    }

    void kill() { m_terminate = true; }
  };


public:

  void test_utility_ostreams() {
    
    vw::null_ostream null_strm;
    vw::multi_ostream multi_strm;

    multi_strm.add(std::cout);
    multi_strm.add(std::cout);
    multi_strm.add(null_strm);

    std::cout << "\nTesting utility ostreams.\n";
    null_strm << "\tYou should not see this message.\n";
    multi_strm << "\tYou should see this message twice.\n";
  }

  void test_basic_logging() {
    std::cout << "\n";

    Log stdout_log(std::cout);
    stdout_log(0,"log test") << "Testing logging to stdout\n";

    Log stderr_log(std::cerr);
    stderr_log(0,"log test") << "Testing logging to stderr\n";
  }

  void test_multithreaded_logging() {
    std::cout << "\n";
    Log log(std::cout);
    log(0, "log test") << "Testing logging from multiple threads\n";
    boost::shared_ptr<TestLogTask> task1( new TestLogTask(log,"log test") );
    boost::shared_ptr<TestLogTask> task2( new TestLogTask(log,"log test") );
    boost::shared_ptr<TestLogTask> task3( new TestLogTask(log,"log test") );
    Thread thread1(task1); 
    Thread thread2(task2); 
    Thread thread3(task3); 

    Thread::sleep_ms(100);
    
    task1->kill();
    task2->kill();
    task3->kill();

    thread1.join();
    thread2.join();
    thread3.join();
    log(0, "log test") << "Log test complete.\n";
  }

  void test_system_log() {
    std::cout << "\nTesting System Log\n";

    vw_out(0) << "\tTesting system log (first call)\n";
    vw_out(0,"test") << "\tTesting system log (second call)\n";
    system_log().add(std::cout);
    SystemLog::system_log()(0,"test") << "\tYou should see this message twice; once with the logging prefix and once without.\n";
    system_log().clear();
    SystemLog::system_log()(0,"test") << "\tYou should see this message once.\n";
  }

};
