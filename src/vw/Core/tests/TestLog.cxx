// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <iostream>
#include <sstream>
#include <string>
#include <map>

#include <vw/Core/Log.h>
#include <vw/Core/Thread.h>
#include <vw/Core/ProgressCallback.h>

#include <gtest/gtest.h>

using namespace vw;

struct TestLogTask {
  bool m_terminate;
  LogInstance& m_log;
  std::string m_ns;
  TestLogTask(LogInstance &log, std::string ns) : m_terminate(false), m_log(log), m_ns(ns) {}

  void operator()() {
    m_log(0,m_ns) << "Start " << Thread::id() << std::endl;

    while( !m_terminate ) {
      Thread::sleep_ms(100);
      m_log(0,m_ns) << "Tick " << Thread::id() << std::endl;
    }

    m_log(0,m_ns) << "Stop " << Thread::id() << std::endl;
  }

  void kill() { m_terminate = true; }
};


TEST(Log, Stringify) {
  EXPECT_EQ( "1234", stringify("1234") );
  EXPECT_EQ( "1234", stringify(1234) );
}

TEST(Log, MultiStream) {

  const static char no_string[]  = "==>You should not see this message.\n";
  const static char yes_string[] = "==>You should see this message twice.\n";

  vw::null_ostream null_strm;
  vw::multi_ostream multi_strm;

  std::ostringstream log;

  multi_strm.add(log);
  multi_strm.add(log);
  multi_strm.add(null_strm);

  null_strm << no_string;
  multi_strm << yes_string;

  multi_strm.remove(null_strm);
  multi_strm.remove(log);
  multi_strm.remove(log);

  std::string expect = std::string(yes_string) + yes_string;
  EXPECT_EQ( expect.size(), log.str().size() );
  EXPECT_EQ( expect, log.str() );
}

TEST(Log, RuleSet) {
  LogRuleSet rs;
  rs.add_rule(vw::InfoMessage, "console");
  rs.add_rule(vw::VerboseDebugMessage, "foo");
  rs.add_rule(vw::EveryMessage, "Bar");

  EXPECT_FALSE( rs(vw::InfoMessage+1, "console"));
  EXPECT_TRUE(  rs(vw::InfoMessage, "console"));
  EXPECT_FALSE( rs(vw::VerboseDebugMessage+1, "foo"));
  EXPECT_TRUE(  rs(vw::VerboseDebugMessage, "foo"));
  EXPECT_TRUE(  rs(vw::VerboseDebugMessage+1, "BAR"));
  EXPECT_TRUE(  rs(vw::VerboseDebugMessage, "BAR"));
}

TEST(Log, BasicLogging) {

  LogInstance stdout_log(std::cout);
  stdout_log(0,"log test") << "Testing logging to stdout" << std::endl;

  LogInstance stderr_log(std::cerr);
  stderr_log(0,"log test") << "Testing logging to stderr" << std::endl;
}

TEST(Log, MultiThreadLog) {
  std::ostringstream stream;

  LogInstance log(stream, false);
  log.rule_set().add_rule(vw::EveryMessage, "log test");
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

  std::istringstream rd;
  rd.str(stream.str());

  std::string typ;
  int id;

  typedef std::map<int, std::string> Map;
  Map status;

  while (rd >> typ && rd >> id) {
    if (typ == "Start") {
      EXPECT_FALSE(status.count(id));
      status[id] = "Start";
    }
    else if (typ == "Tick") {
      ASSERT_TRUE(status.count(id));
      EXPECT_EQ( "Start", status[id] );
    }
    else if (typ == "Stop") {
      ASSERT_TRUE(status.count(id));
      EXPECT_EQ( "Start", status[id] );
      status[id] = "Stop";
    }
    else
      FAIL() << "Unknown message [" << typ << "]";
  }

  EXPECT_EQ( 3, status.size() );

  for (Map::const_iterator i = status.begin(), end = status.end(); i != end; ++i) {
    EXPECT_EQ("Stop", i->second) << "Thread " << i->first << " out of order";
  }

  EXPECT_FALSE(rd);
}

TEST(Log, SystemLog) {
  vw_log().console_log().rule_set().add_rule(vw::EveryMessage, "test");

  vw_out(0) << "\tTesting system log (first call)\n";
  vw_out(0,"test") << "\tTesting system log (second call)\n";

  boost::shared_ptr<LogInstance> new_log(new LogInstance(std::cout));
  new_log->rule_set().add_rule(vw::EveryMessage, "test");
  vw_log().add(new_log);
  vw_out(0,"test") << "\tYou should see this message twice; once with the logging prefix and once without.\n";

  vw_log().clear();
  vw_out(0,"test") << "\tYou should see this message once.\n";
}

TEST(Log, ProgressCallback) {
  vw_out(0) << "\nTesting Logging with a progress callback\n";
  TerminalProgressCallback pc(vw::InfoMessage, "\tTesting: ");
  for (double i = 0; i < 1.0; i+=0.01) {
    pc.report_progress(i);
    Thread::sleep_ms(10);
  }
  pc.report_finished();
}

TEST(Log, HiresProgressCallback) {
  vw_out(0) << "\nTesting Logging with a progress callback\n";
  TerminalProgressCallback pc(vw::InfoMessage, "\tTesting: ", 2);
  for (int i = 0; i < 10000; ++i) {
    pc.report_progress(i/10000.0);
    if (i % 50 == 0)
      Thread::sleep_ms(10);
  }
  pc.report_finished();
}

TEST(Log, FlushAndNewline) {
  std::ostringstream stream;
  LogInstance log(stream);

  log(0) << "\nTesting log termination operators.\n";
  log(0) << "\tTesting log line terminated by std::flush..." << std::flush;
  Thread::sleep_ms(1000);
  log(0) << " done.\n";
  log(0) << "\tTesting log line terminated by std::endl... done." << std::endl;
  Thread::sleep_ms(1000);
  log(0) << "Finished.\n";
}
