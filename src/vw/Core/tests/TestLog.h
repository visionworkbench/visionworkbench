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
#include <vw/Core/ProgressCallback.h>

#include <iostream>
#include <fstream>
#include <string>

using namespace vw;

class TestLog : public CxxTest::TestSuite
{

  struct TestLogTask {
    bool m_terminate;
    LogInstance& m_log;
    std::string m_ns;
    TestLogTask(LogInstance &log, std::string ns) : m_terminate(false), m_log(log), m_ns(ns) {}

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

  void test_stringify() {
    TS_ASSERT_SAME_DATA(stringify("1234").c_str(), "1234", 4);
    TS_ASSERT_SAME_DATA(stringify(1234).c_str(), "1234", 4);
  }

  void test_utility_ostreams() {

    vw::null_ostream null_strm;
    vw::multi_ostream multi_strm;

    std::fstream log;

    log.open("log.txt", std::ios::out | std::ios::binary );
    multi_strm.add(log);
    multi_strm.add(log);
    multi_strm.add(null_strm);

    TS_TRACE("Testing utility ostreams.");
    null_strm << "==>You should not see this message.\n";
    multi_strm << "==>You should see this message twice.\n";

    multi_strm.remove(null_strm);
    multi_strm.remove(log);
    multi_strm.remove(log);

    log.close();

    log.open("log.txt", std::ios::in | std::ios::binary);

    char buf[38];

    log.read(buf, 38);
    TS_ASSERT_SAME_DATA(buf, "==>You should see this message twice.\n", 38); // 38 = string + newline

    log.read(buf, 38);
    TS_ASSERT_SAME_DATA(buf, "==>You should see this message twice.\n", 38); // 38 = string + newline

    memset(buf, 0, 38);

    log.read(buf, 38);
    TS_ASSERT_EQUALS(buf[0], '\0');
    TS_ASSERT(log.eof());

    log.close();
  }

  void test_log_rule_set() {
    LogRuleSet rs;
    rs.add_rule(vw::InfoMessage, "console");
    rs.add_rule(vw::VerboseDebugMessage, "foo");
    rs.add_rule(vw::EveryMessage, "Bar");

    TS_ASSERT_EQUALS( rs(vw::InfoMessage+1, "console"), false);
    TS_ASSERT_EQUALS( rs(vw::InfoMessage, "console"), true);
    TS_ASSERT_EQUALS( rs(vw::VerboseDebugMessage+1, "foo"), false);
    TS_ASSERT_EQUALS( rs(vw::VerboseDebugMessage, "foo"), true);
    TS_ASSERT_EQUALS( rs(vw::VerboseDebugMessage+1, "BAR"), true);
    TS_ASSERT_EQUALS( rs(vw::VerboseDebugMessage, "BAR"), true);
  }

  void test_basic_logging() {

    LogInstance stdout_log(std::cout);
    stdout_log(0,"log test") << "Testing logging to stdout" << std::endl;

    LogInstance stderr_log(std::cerr);
    stderr_log(0,"log test") << "Testing logging to stderr" << std::endl;
  }

  void test_multithreaded_logging() {
    LogInstance log(std::cout);
    log.rule_set().add_rule(vw::EveryMessage, "log test");
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
    Log::system_log().console_log().rule_set().add_rule(vw::EveryMessage, "test");

    vw_out(0) << "\tTesting system log (first call)\n";
    vw_out(0,"test") << "\tTesting system log (second call)\n";

    boost::shared_ptr<LogInstance> new_log(new LogInstance(std::cout));
    new_log->rule_set().add_rule(vw::EveryMessage, "test");
    Log::system_log().add(new_log);
    vw_out(0,"test") << "\tYou should see this message twice; once with the logging prefix and once without.\n";

    Log::system_log().clear();
    vw_out(0,"test") << "\tYou should see this message once.\n";
  }

  void test_progress_callback() {
     vw_out(0) << "\nTesting Logging with a progress callback\n";
    TerminalProgressCallback pc(vw::InfoMessage, "\tTesting: ");
    for (double i = 0; i < 1.0; i+=0.01) {
      pc.report_progress(i);
      Thread::sleep_ms(10);
    }
    pc.report_finished();
  }

  void test_flush_and_newline() {
    vw_out(0) << "\nTesting log termination operators.\n";
    vw_out(0) << "\tTesting log line terminated by std::flush..." << std::flush;
    Thread::sleep_ms(1000);
    vw_out(0) << " done.\n";
    vw_out(0) << "\tTesting log line terminated by std::endl... done." << std::endl;
    Thread::sleep_ms(1000);
    vw_out(0) << "Finished.\n";
  }

};
