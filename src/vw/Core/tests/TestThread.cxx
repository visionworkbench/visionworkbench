// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>

#include <vw/Core/Thread.h>

using namespace vw;


int once_value = 0;
void run_once_func() {
  ++once_value;
}
RunOnce once = VW_RUNONCE_INIT;

class TestTask {
  public:
    int value;
    bool terminate;
    Mutex &m_mutex;
    Condition &m_condition;
    TestTask(Mutex &mutex, Condition &cond) : value(0), terminate(false), m_mutex(mutex), m_condition(cond) {}

    void operator()() {
      {
        Mutex::Lock lock(m_mutex);
        value = 1;
      }
      m_condition.notify_all();

      int count = 0;
      while( !terminate ) {
        Thread::sleep_ms(100);
        if (count++ > 50) {
          std::cerr << "Test thread iterated 100 times... it shouldn't take this long.  Maybe deadlock occured?" << std::endl;
          exit(-1);
        }
      }

      // Terminate, but wait 100ms first, to make sure our condition
      // synchronization is working.
      Thread::sleep_ms(10);
      {
        Mutex::Lock lock(m_mutex);
        value = 3;
      }
      m_condition.notify_all();
    }

    void kill() { Thread::sleep_ms(100); terminate = true; }
};

class TestThreadIdTask {
  public:
    bool terminate;
    TestThreadIdTask() : terminate(false) {}

    void operator()() {
      while( !terminate )
        Thread::sleep_ms(100);
    }

    void kill() { Thread::sleep_ms(100); terminate = true; }
};

TEST(Thread, ThreadAndMutex) {
  Mutex condition_mutex;
  Condition index_updated_event;
  boost::shared_ptr<TestTask> task( new TestTask(condition_mutex, index_updated_event) );
  ASSERT_EQ( 0, task->value );

  Thread *thread;

  {
    Mutex::Lock cond_lock(condition_mutex);
    thread = new Thread( task );
    index_updated_event.wait(cond_lock);
  }

  EXPECT_EQ( 1, task->value );

  task->kill();
  {
    Mutex::Lock cond_lock(condition_mutex);
    index_updated_event.wait(cond_lock);
  }

  EXPECT_EQ( 3, task->value );
  thread->join();
  delete thread;
}

TEST(Thread, ThreadId) {

  // Start by testing my own ID.
  int64 id;
  id = Thread::id();
  ASSERT_EQ( 0, id );

  // Make sure that it stays the same after repeated calls
  id = Thread::id();
  EXPECT_EQ( 0, id );
  id = Thread::id();
  EXPECT_EQ( 0, id );

  boost::shared_ptr<TestThreadIdTask> task1( new TestThreadIdTask() );
  boost::shared_ptr<TestThreadIdTask> task2( new TestThreadIdTask() );
  boost::shared_ptr<TestThreadIdTask> task3( new TestThreadIdTask() );
  Thread thread1(task1);
  Thread thread2(task2);
  Thread thread3(task3);

  task1->kill();
  task2->kill();
  task3->kill();

  thread1.join();
  thread2.join();
  thread3.join();
}

TEST(Thread, RunOnce) {
  EXPECT_EQ( 0, once_value );
  once.run( run_once_func );
  EXPECT_EQ( 1, once_value );
  once.run( run_once_func );
  EXPECT_EQ( 1, once_value );
}

struct ReturnFalse {
  bool operator()() { return false; }
};
struct ReturnTrue {
  bool operator()() { return true; }
};

TEST(Thread, TimedWait) {
  Mutex m;
  Mutex::Lock lock(m);
  Condition c;
  EXPECT_FALSE(c.timed_wait(m, 0));
  EXPECT_FALSE(c.timed_wait(m, 0, ReturnFalse()));
  EXPECT_TRUE (c.timed_wait(m, 0, ReturnTrue()));
}
