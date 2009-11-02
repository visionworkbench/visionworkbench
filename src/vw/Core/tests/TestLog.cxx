// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
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

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/algorithm/string.hpp>

#include <gtest/gtest.h>

using namespace vw;

struct raii {
  typedef boost::function<void (void)> FuncT;
  FuncT m_leave;
  raii(FuncT enter, FuncT leave) : m_leave(leave) {enter();};
  ~raii() {m_leave();}
};

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

  EXPECT_EQ( 3u, status.size() );

  for (Map::const_iterator i = status.begin(), end = status.end(); i != end; ++i) {
    EXPECT_EQ("Stop", i->second) << "Thread " << i->first << " out of order";
  }

  EXPECT_FALSE(rd);
}

TEST(Log, SystemLog) {

  std::ostringstream sstr;

  raii fix(boost::bind(&Log::set_console_stream, boost::ref(vw_log()), boost::ref(sstr),      LogRuleSet(),  false),
           boost::bind(&Log::set_console_stream, boost::ref(vw_log()), boost::ref(std::cout), LogRuleSet(), false));

  vw_log().console_log().rule_set().add_rule(vw::EveryMessage, "test");

  vw_out(0) << "\tTesting system log (first call)\n";
  vw_out(0,"test") << "\tTesting system log (second call)\n";

  boost::shared_ptr<LogInstance> new_log(new LogInstance(sstr));
  new_log->rule_set().add_rule(vw::EveryMessage, "test");
  vw_log().add(new_log);
  vw_out(0,"test") << "\tYou should see this message twice; once with the logging prefix and once without.\n";

  vw_log().clear();
  vw_out(0,"test") << "\tYou should see this message once.\n";

  const std::string &out = sstr.str();
  std::vector<std::string> lines;
  boost::split(lines, out, boost::is_any_of("\n"));

  EXPECT_EQ(6u, lines.size());
  EXPECT_EQ("\tTesting system log (first call)", lines[0]);
  EXPECT_EQ("\tTesting system log (second call)", lines[1]);

  // Technically, these can arrive in any order, but in practice, it's the order the rules were added
  EXPECT_EQ("\tYou should see this message twice; once with the logging prefix and once without.", lines[2]);

  EXPECT_EQ("\tYou should see this message once.", lines[4]);
  EXPECT_EQ("", lines[5]);

}

TEST(Log, ProgressCallback) {
  std::ostringstream sstr;

  raii fix(boost::bind(&vw::set_output_stream, boost::ref(sstr)),
           boost::bind(&vw::set_output_stream, boost::ref(std::cout)));

  vw_out(0) << "\nTesting Logging with a progress callback\n";
  TerminalProgressCallback pc(vw::InfoMessage, "\tTesting: ");
  for (double i = 0; i < 1.0; i+=0.01) {
    pc.report_progress(i);
    Thread::sleep_ms(0);
  }
  pc.report_finished();

  const std::string &out = sstr.str();
  EXPECT_GT(out.size(), 0u);
  EXPECT_TRUE(boost::iends_with(out, std::string("\r\tTesting: [************************************************************] Complete!\n")));
}

TEST(Log, HiresProgressCallback) {
  std::ostringstream sstr;

  raii fix(boost::bind(&vw::set_output_stream, boost::ref(sstr)),
           boost::bind(&vw::set_output_stream, boost::ref(std::cout)));

  vw_out(0) << "\nTesting Logging with a progress callback\n";
  TerminalProgressCallback pc(vw::InfoMessage, "\tTesting: ", 2);
  for (int i = 0; i < 10000; ++i) {
    pc.report_progress(i/10000.0);
    if (i % 50 == 0)
      Thread::sleep_ms(0);
  }
  pc.report_finished();

  const std::string &out = sstr.str();
  EXPECT_GT(out.size(), 0u);
  EXPECT_TRUE(boost::iends_with(out, std::string("\r\tTesting: [************************************************************] Complete!\n")));
}

TEST(Log, FlushAndNewline) {
  std::ostringstream stream, end;
  LogInstance log(stream, false);

  // Capture the EOL character, since it differs per-platform
  end << std::endl;

  const std::string
              s1("\nTesting log termination operators.\n"),
              s2("\tTesting log line terminated by std::flush..."),
              s3(" done.\n"),
              s4("\tTesting log line terminated by std::endl... done."),
              s5(end.str()),
              s6("Finished.\n");

  // Make sure a newline flushes
  log(0) << s1;
  EXPECT_EQ(s1, stream.str());

  // Make sure no newline does not flush
  log(0) << s2;
  EXPECT_EQ(s1, stream.str());

  // Explicit flush
  log(0) << std::flush;
  EXPECT_EQ(s1 + s2, stream.str());

  // Another newline
  log(0) << s3;
  EXPECT_EQ(s1 + s2 + s3, stream.str());

  // No newline again
  log(0) << s4;
  EXPECT_EQ(s1 + s2 + s3, stream.str());

  // And try to flush with endl
  log(0) << std::endl;
  EXPECT_EQ(s1 + s2 + s3 + s4 + s5, stream.str());

  log(0) << s6;
  EXPECT_EQ(s1 + s2 + s3 + s4 + s5 + s6, stream.str());
}
