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

// TestFunctors.h
#include <cxxtest/TestSuite.h>

#include <vw/Core/Thread.h>

using namespace vw;


int once_value = 0;
void run_once_func() {
  ++once_value;
}
RunOnce once = VW_RUNONCE_INIT;

class TestThread : public CxxTest::TestSuite
{

  class TestTask {
  public:
    int value;
    bool terminate;
    Condition &m_condition;
    TestTask(Condition &cond) : value(0), terminate(false), m_condition(cond) {}

    void operator()() {
      value = 1;
      m_condition.notify_all();

      int count = 0;
      while( !terminate ) {
        Thread::yield();
        if (count++ > 100) {
          TS_FAIL("Test thread iterated 100 times... it shouldn't take this long.  Maybe deadlock occured?");
          exit(0);
        }
      }

      // Terminate, but wait 100ms first, to make sure our condition
      // synchronization is working.
      Thread::sleep_ms(100);
      value = 3;
      m_condition.notify_all();
    }

    void kill() {
      terminate = true;
    }
  };
      
public:
  void test_thread_and_mutex()
  {
    Mutex condition_mutex;
    Condition index_updated_event;
    boost::shared_ptr<TestTask> task( new TestTask(index_updated_event) );
    TS_ASSERT_EQUALS( task->value, 0 );

    // Okay, so this is not how you'd normally use a lock.
    Thread *thread = new Thread( task );

    {
      Mutex::Lock cond_lock(condition_mutex);
      index_updated_event.wait(cond_lock);
    }
    TS_ASSERT_EQUALS( task->value, 1 );

    task->kill();
    {
      Mutex::Lock cond_lock(condition_mutex);
      index_updated_event.wait(cond_lock);
    }
    TS_ASSERT_EQUALS( task->value, 3 );
    thread->join();
    delete thread;
  }

  void test_run_once()
  {
    TS_ASSERT_EQUALS( once_value, 0 );
    once.run( run_once_func );
    TS_ASSERT_EQUALS( once_value, 1 );
    once.run( run_once_func );
    TS_ASSERT_EQUALS( once_value, 1 );
  }
};
