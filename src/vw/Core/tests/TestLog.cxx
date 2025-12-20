// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__

#include <gtest/gtest_VW.h>             // for EXPECT_EQ, TEST, etc
#include <vw/Core/FundamentalTypes.h>   // for uint64
#include <vw/Core/Log.h>                // for LogRuleSet, LogInstance, etc
#include <vw/Core/ProgressCallback.h>   // for TerminalProgressCallback
#include <vw/Core/Stopwatch.h>          // for Stopwatch
#include <vw/Core/Thread.h>             // for Thread
#include <vw/Core/System.h>             // for vw_log

#include <stddef.h>                     // for size_t
#include <stdlib.h>                     // for exit
#include <unistd.h>                     // for sleep
#include <iostream>                     // for ostringstream, operator<<, etc
#include <map>                          // for _Rb_tree_const_iterator, etc
#include <sstream>
#include <string>                       // for string, operator+, etc
#include <utility>                      // for pair, make_pair
#include <vector>                       // for vector

#include <boost/algorithm/string/classification.hpp>  // for is_any_of
#include <boost/algorithm/string/detail/classification.hpp>
#include <boost/algorithm/string/predicate.hpp>  // for iends_with
#include <boost/algorithm/string/split.hpp>  // for split
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/iterator/iterator_facade.hpp>  // for operator!=
#include <boost/ref.hpp>               // for reference_wrapper, ref
#include <boost/shared_ptr.hpp>  // for shared_ptr

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
    const std::string start("Start " + stringify(Thread::id()) + "\n"),
                        tick("Tick " + stringify(Thread::id()) + "\n"),
                        stop("Stop " + stringify(Thread::id()) + "\n");

    m_log(InfoMessage,m_ns) << start;

    while( !m_terminate ) {
      Thread::sleep_ms(100);
      m_log(InfoMessage,m_ns) << tick;
    }

    m_log(InfoMessage,m_ns) << stop;
  }

  void kill() { m_terminate = true; }
};


TEST(Log, Stringify) {
  EXPECT_EQ( "1234", stringify("1234") );
  EXPECT_EQ( "1234", stringify(1234) );
}

/*
TEST(Log, MultiStream) {

  const static char no_string[]  = "==>You should not see this message.\n";
  const static char yes_string[] = "==>You should see this message twice.\n";

  null_ostream null_strm;
  multi_ostream multi_strm;

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
*/

TEST(Log, RuleSet) {
  LogRuleSet rs;
  rs.clear();

  // *   matches anything
  // *.a matches [first.a, second.a]
  // a   matches just a (ie, no namespace)
  // a.* matches [a, a.first, a.first.second]

  rs.add_rule(InfoMessage, "console");
  rs.add_rule(VerboseDebugMessage, "verbose");
  rs.add_rule(EveryMessage, "every");

  rs.add_rule(EveryMessage, "unqual");
  rs.add_rule(EveryMessage, "prefix.*");
  rs.add_rule(EveryMessage, "*.suffix");

  const MessageLevel v = VerboseDebugMessage;

  EXPECT_FALSE( rs(InfoMessage+1, "console"));
  EXPECT_TRUE(  rs(InfoMessage, "console"));
  EXPECT_FALSE( rs(v+1, "verbose"));
  EXPECT_TRUE(  rs(v, "verbose"));
  EXPECT_TRUE(  rs(v+1, "every"));
  EXPECT_TRUE(  rs(v, "every"));

  EXPECT_TRUE( rs(v, "unqual"));
  EXPECT_FALSE(rs(v, "unqual.foo"));
  EXPECT_TRUE( rs(v, "prefix"));
  EXPECT_TRUE( rs(v, "prefix.foo"));
  EXPECT_FALSE(rs(v, "foo.prefix"));
  EXPECT_FALSE(rs(v, "foo.prefix.bar"));
  EXPECT_FALSE(rs(v, "suffix"));
  EXPECT_TRUE( rs(v, "foo.suffix"));
  EXPECT_TRUE( rs(v, "foo.bar.suffix"));
  EXPECT_FALSE(rs(v, "suffix.foo"));
  EXPECT_FALSE(rs(v, "bar.suffix.foo"));
  EXPECT_FALSE(rs(v, "foo.falsesuffix"));
  EXPECT_FALSE(rs(v,  "falseprefix"));
  EXPECT_FALSE(rs(v,  "falseprefix.foo"));
  EXPECT_FALSE(rs(v,  "prefixfalse.foo"));

  LogRuleSet all;
  all.add_rule(EveryMessage, "*");
  all.add_rule(ErrorMessage, "console");
  all.add_rule(WarningMessage, "off");

  EXPECT_FALSE(all(v, "console"));
  EXPECT_TRUE( all(v, "any"));
  EXPECT_FALSE(all(v, "off"));
  EXPECT_TRUE( all(v, "any.all"));
}

TEST(Log, RuleSetIllegal) {
  LogRuleSet rs;
  EXPECT_THROW(rs.add_rule(VerboseDebugMessage, "*.foo.*"), vw::ArgumentErr);
  EXPECT_THROW(rs.add_rule(VerboseDebugMessage, "foo.*.bar"), vw::ArgumentErr);
}

TEST(Log, BasicLogging) {

  LogInstance stdout_log(std::cout);
  stdout_log(InfoMessage,"log test") << "Testing logging to stdout" << std::endl;

  LogInstance stderr_log(std::cerr);
  stderr_log(InfoMessage,"log test") << "Testing logging to stderr" << std::endl;
}

TEST(Log, MultiThreadLog) {
  std::ostringstream stream;

  LogInstance log(stream, false);
  log.rule_set().add_rule(EveryMessage, "log test");

  typedef boost::shared_ptr<TestLogTask> TheTask;
  typedef boost::shared_ptr<Thread>  TheThread;
  std::vector<std::pair<TheTask, TheThread> > threads(100);

  for (size_t i = 0; i < threads.size(); ++i) {
    TheTask task(   new TestLogTask(log,"log test") );
    TheThread thread( new Thread(task) );
    threads[i] = std::make_pair(task, thread);
  }

  Thread::sleep_ms(100);

  for (size_t i = 0; i < threads.size(); ++i) {
    threads[i].first->kill();
  }

  for (size_t i = 0; i < threads.size(); ++i) {
    threads[i].second->join();
  }

  std::istringstream rd;
  std::string buf = stream.str();
  rd.str(buf);

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

  EXPECT_EQ( threads.size(), status.size() );

  for (Map::const_iterator i = status.begin(), end = status.end(); i != end; ++i) {
    EXPECT_EQ("Stop", i->second) << "Thread " << i->first << " out of order";
  }

  EXPECT_FALSE(rd);
}

TEST(Log, SystemLog) {

  const std::string
    s1("\tTesting system log (first call)"),
    s2("\tTesting system log (second call)"),
    s3("\tYou should see this message twice; once with the logging prefix and once without."),
    s4("\tYou should see this message once.");

  std::ostringstream sstr, sstr2;

  raii fix(boost::bind(&Log::set_console_stream, boost::ref(vw_log()), boost::ref(sstr),      LogRuleSet(),  false),
           boost::bind(&Log::set_console_stream, boost::ref(vw_log()), boost::ref(std::cout), LogRuleSet(), false));

  vw_log().console_log().rule_set().add_rule(EveryMessage, "test");

  vw_out() << s1 << "\n";
  vw_out(DebugMessage,"test") << s2 << "\n";

  boost::shared_ptr<LogInstance> new_log(new LogInstance(sstr));
  new_log->rule_set().add_rule(EveryMessage, "test");
  vw_log().add(new_log);
  vw_log().add(sstr2, new_log->rule_set());

  vw_out(DebugMessage,"test") << s3 << "\n";

  vw_log().clear();
  vw_out(DebugMessage,"test") << s4 << "\n";

  const std::string out1 = sstr.str(), out2 = sstr2.str();

  std::vector<std::string> lines, lines2;
  boost::split(lines,  out1, boost::is_any_of("\n"));
  boost::split(lines2, out2, boost::is_any_of("\n"));

  ASSERT_EQ(6u, lines.size());
  EXPECT_EQ(s1, lines[0]);
  EXPECT_EQ(s2, lines[1]);
  // Technically, 2 and 3 can arrive in any order, but in practice, it's the order the rules were added
  EXPECT_EQ(s3, lines[2]);
  EXPECT_EQ(s3, lines[3].substr(35)); // remove log prefix
  EXPECT_EQ(s4, lines[4]);
  EXPECT_EQ("", lines[5]);

  ASSERT_EQ(2u, lines2.size());
  EXPECT_EQ(s3, lines2[0].substr(35)); // remove log prefix
  EXPECT_EQ("", lines2[1]);
}

TEST(Log, ProgressCallback) {
  std::ostringstream sstr;

  raii fix(boost::bind(&Log::set_console_stream, boost::ref(vw_log()), boost::ref(sstr),      LogRuleSet(),  false),
           boost::bind(&Log::set_console_stream, boost::ref(vw_log()), boost::ref(std::cout), LogRuleSet(), false));

  TerminalProgressCallback pc( "test", "\tTesting: ");
  for (double i = 0; i < 1.0; i+=0.01) {
    pc.report_progress(i);
    Thread::sleep_ms(0);
  }
  pc.report_finished();

  const std::string &out = sstr.str();
  EXPECT_GT(out.size(), 0u);
  size_t last_line_idx = out.rfind("\r");
  EXPECT_EQ( 80u, out.size()-last_line_idx-2 );
  EXPECT_TRUE(boost::iends_with(out, std::string("***]\n")));
  EXPECT_THROW(TerminalProgressCallback("monkey","monkeymonkeymonkeymonkeymonkeymonkeymonkeymonkeymonkeymonkeymonkeymonkeymonkeymonkeymonkeymonkeymonkeymonkey"), ArgumentErr );
}

TEST(Log, HiresProgressCallback) {
  std::ostringstream sstr;

  raii fix(boost::bind(&Log::set_console_stream, boost::ref(vw_log()), boost::ref(sstr),      LogRuleSet(),  false),
           boost::bind(&Log::set_console_stream, boost::ref(vw_log()), boost::ref(std::cout), LogRuleSet(), false));

  TerminalProgressCallback pc( "test", "\tTesting: ", InfoMessage, 2);
  for (int i = 0; i < 10000; ++i) {
    pc.report_progress(i/10000.0);
    if (i % 50 == 0)
      Thread::sleep_ms(0);
  }
  pc.report_finished();

  const std::string &out = sstr.str();
  EXPECT_GT(out.size(), 0u);
  size_t last_line_idx = out.rfind("\r");
  EXPECT_EQ( 80u, out.size()-last_line_idx-2 );
  EXPECT_TRUE(boost::iends_with(out, std::string("***]\n")));
  EXPECT_THROW(TerminalProgressCallback("monkey","monkeymonkeymonkeymonkeymonkeymonkeymonkeymonkeymonkeymonkeymonkeymonkeymonkeymonkeymonkeymonkeymonkeymonkey"), ArgumentErr );
}

TEST(Log, ProgressHide) {

  std::ostringstream sstr;

  raii fix(boost::bind(&Log::set_console_stream, boost::ref(vw_log()), boost::ref(sstr),      LogRuleSet(),  false),
           boost::bind(&Log::set_console_stream, boost::ref(vw_log()), boost::ref(std::cout), LogRuleSet(), false));

  // Progress bars hidden
  vw_log().console_log().rule_set().add_rule(vw::NoMessage, "*.progress");

  // Can't see progress bar
  TerminalProgressCallback pc( "test", "Rawr:" );
  pc.report_progress(.1);
  pc.report_progress(.2);

  EXPECT_EQ( sstr.str().size(), 0u );

  vw_log().console_log().rule_set().clear();

  // Can see progress bar
  pc.report_progress(.5);
  pc.report_finished();

  EXPECT_GT( sstr.str().size(), 0u );
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
  log(InfoMessage) << s1;
  EXPECT_EQ(s1, stream.str());

  // Make sure no newline does not flush
  log(InfoMessage) << s2;
  EXPECT_EQ(s1, stream.str());

  // Explicit flush
  log(InfoMessage) << std::flush;
  EXPECT_EQ(s1 + s2, stream.str());

  // Another newline
  log(InfoMessage) << s3;
  EXPECT_EQ(s1 + s2 + s3, stream.str());

  // No newline again
  log(InfoMessage) << s4;
  EXPECT_EQ(s1 + s2 + s3, stream.str());

  // And try to flush with endl
  log(InfoMessage) << std::endl;
  EXPECT_EQ(s1 + s2 + s3 + s4 + s5, stream.str());

  log(InfoMessage) << s6;
  EXPECT_EQ(s1 + s2 + s3 + s4 + s5 + s6, stream.str());
}

TEST(LogDeathTest, MakeSureDeathTestsWork) {
  LogInstance log(std::cerr, false);
  EXPECT_EXIT(log(ErrorMessage) << "Rawr" << std::flush; exit(12), ::testing::ExitedWithCode(12), "Rawr");
}

// This test is tied to bug #199. Which is still a bug.
TEST(LogDeathTest, DISABLED_FlushOnExit) {
  LogInstance log(std::cerr, false);
  EXPECT_EXIT(log(ErrorMessage) << "Rawr"; exit(12), ::testing::ExitedWithCode(12), "Rawr");
}

const std::string& slow(int len, const std::string& msg) {
  sleep(len);
  return msg;
}

// This tests whether the loggers are lazy or not (whether they evaluate their
// arguments if they're turned off). This may never work, but here's a test for it!
TEST(Log, DISABLED_LazyLog) {
  std::ostringstream stream, ends;
  ends << std::endl;

  LogInstance log(stream, false);

  std::string s1("A"), s2("B"), s3("C"), end(ends.str());

  uint64 start = vw::Stopwatch::microtime();
  log(InfoMessage)         << s1          << std::endl;
  log(VerboseDebugMessage) << slow(1, s2) << std::endl;
  log(InfoMessage)         << s3          << std::endl;
  uint64 stop = vw::Stopwatch::microtime();

  ASSERT_EQ(s1 + end + s3 + end, stream.str());
  // the sleep inside slow() should not have run.
  EXPECT_LT(stop-start, 200000u);
}

TEST(Log, WarningsOnByDefault) {
  std::ostringstream sstr;
  LogInstance log(sstr, false);

  log(WarningMessage) << "foo\n";
  log(WarningMessage, "someothernamespace") << "bar\n";

  const std::string& x = sstr.str();
  ASSERT_FALSE(x.empty());
  EXPECT_EQ("Warning: foo\nWarning: bar\n", x);
}

TEST(Log, WarningsCanBeTurnedOff1) {
  std::ostringstream sstr;
  LogInstance log(sstr, false);
  log.rule_set().add_rule(ErrorMessage, "console");

  log(WarningMessage) << "foo\n";
  log(WarningMessage, "someothernamespace") << "bar\n";

  const std::string& x = sstr.str();
  ASSERT_FALSE(x.empty());
  EXPECT_EQ("Warning: bar\n", x);
}

TEST(Log, WarningsCanBeTurnedOff2) {
  std::ostringstream sstr;
  LogInstance log(sstr, false);
  log.rule_set().add_rule(ErrorMessage, "someothernamespace");

  log(WarningMessage) << "foo\n";
  log(WarningMessage, "someothernamespace") << "bar\n";

  const std::string& x = sstr.str();
  ASSERT_FALSE(x.empty());
  EXPECT_EQ("Warning: foo\n", x);
}

TEST(Log, WarningsCanBeTurnedOff3) {
  std::ostringstream sstr;
  LogInstance log(sstr, false);
  log.rule_set().add_rule(ErrorMessage, "someothernamespace");
  log.rule_set().add_rule(ErrorMessage, "console");

  log(WarningMessage) << "foo\n";
  log(WarningMessage, "someothernamespace") << "bar\n";

  const std::string& x = sstr.str();
  ASSERT_TRUE(x.empty());
}

TEST(Log, WarningsCanBeTurnedOff4) {
  std::ostringstream sstr;
  LogInstance log(sstr, false);
  log.rule_set().add_rule(ErrorMessage, "*");

  log(WarningMessage) << "foo\n";
  log(WarningMessage, "someothernamespace") << "bar\n";

  const std::string& x = sstr.str();
  ASSERT_TRUE(x.empty());
}

TEST(Log, LazyIsEnabled) {
  vw_log().console_log().rule_set().add_rule(InfoMessage, "test");

  EXPECT_TRUE(  vw_log().is_enabled() );
  EXPECT_FALSE( vw_log().is_enabled(DebugMessage) );
  EXPECT_TRUE(  vw_log().is_enabled(InfoMessage,"test") );
  EXPECT_FALSE( vw_log().is_enabled(DebugMessage,"test") );
  EXPECT_FALSE( vw_log().is_enabled(InfoMessage,"other") );

  std::ostringstream sstr;
  boost::shared_ptr<LogInstance> new_log(new LogInstance(sstr));
  new_log->rule_set().add_rule(EveryMessage, "other");
  vw_log().add(new_log);

  EXPECT_TRUE(  vw_log().is_enabled(InfoMessage,"other") );
  EXPECT_TRUE(  vw_log().is_enabled(VerboseDebugMessage,"other") );

  // How it should work with a macro
  std::ostringstream sstr2;
  boost::shared_ptr<LogInstance> new_log2(new LogInstance(sstr2, false));
  new_log2->rule_set().clear();
  new_log2->rule_set().add_rule(InfoMessage,"other2");
  vw_log().add(new_log2);

  VW_OUT(InfoMessage,"other2") << "Printed" << std::flush;
  VW_OUT(DebugMessage,"other2") << "NotPrinted" << std::flush;
  EXPECT_EQ( "Printed", sstr2.str() );

  // This would then be lazy
  uint64 start = vw::Stopwatch::microtime();
  VW_OUT(DebugMessage,"other2") << slow(1, "Not Printed") << std::flush;
  uint64 stop  = vw::Stopwatch::microtime();
  EXPECT_LT(stop-start, 200000u);

  // Make sure default options still work with the variadic
  // macro. This would fail to compile if it didn't work.
  raii fix(boost::bind(&Log::set_console_stream, boost::ref(vw_log()), boost::ref(sstr),      LogRuleSet(),  false),
           boost::bind(&Log::set_console_stream, boost::ref(vw_log()), boost::ref(std::cout), LogRuleSet(), false));
  VW_OUT() << "You should see me once\n";
}
