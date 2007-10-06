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
  class TestRunnable : public Thread::Runnable {
  public:
    int value;
    bool terminate;
    Mutex mutex;
    TestRunnable() : value(0), terminate(false) {}
    void run() {
      value = 1;
      {
        Mutex::Lock lock(mutex);
        value = 2;
      }
      while( !terminate ) {
        Thread::yield();
      }
      value = 3;
    }
    void kill() {
      terminate = true;
    }
  };
      
public:
  void test_thread_and_mutex()
  {
    boost::shared_ptr<TestRunnable> runnable( new TestRunnable() );
    TS_ASSERT_EQUALS( runnable->value, 0 );

    // Okay, so this is not how you'd normally use a lock.
    Mutex::Lock *lock = new Mutex::Lock( runnable->mutex );
    Thread *thread = new Thread( runnable );
    // There is technically a race condition in this test, but this 
    // yield makes things reasonably deterministic.
    Thread::yield();
    TS_ASSERT_EQUALS( runnable->value, 1 );

    delete lock;
    // This is overkill: we're just testing both APIs
    Thread::yield();
    Thread::sleep_ms(10);
    TS_ASSERT_EQUALS( runnable->value, 2 );

    delete thread;
    TS_ASSERT_EQUALS( runnable->value, 3 );
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
