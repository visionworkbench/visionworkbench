// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestThreadPool.h
#include <cxxtest/TestSuite.h>

#include <vw/Core/ThreadPool.h>

#include <iostream>

using namespace vw;

class TestThread : public CxxTest::TestSuite
{

  class TestTask : public Task {

    TestTask(TestTask& copy) {}
    void operator=(TestTask& copy) {}
  public:
    int m_value;
    bool terminate;
    Mutex test;
    TestTask() : m_value(0), terminate(false) {}

    void operator()() {
      m_value = 1;

      int count = 0;
      while( !terminate ) {
        Thread::sleep_ms(10);
        if (count++ > 500) {
          TS_FAIL("Test thread iterated 100 times... it shouldn't take this long.  Maybe deadlock occured?");
          exit(0);
        }
      }
      m_value = 3;
    }

    int value() { return m_value; }

    void kill() {
      terminate = true;
    }
  };

public:
  void test_threadpool()
  {
    boost::shared_ptr<TestTask> task1 (new TestTask);
    boost::shared_ptr<TestTask> task2 (new TestTask);
    boost::shared_ptr<TestTask> task3 (new TestTask);

    TS_ASSERT_EQUALS( task1->value(), 0 );
    TS_ASSERT_EQUALS( task2->value(), 0 );
    TS_ASSERT_EQUALS( task3->value(), 0 );

    set_debug_level(VerboseDebugMessage);

    FifoWorkQueue queue(4);
    queue.add_task(task1);
    queue.add_task(task2);
    queue.add_task(task3);

    // This is slightly sloppy -- we just wait a few hundred ms for
    // the thread to start and to set the flag.  This should really be
    // done with a condition variable.
    Thread::sleep_ms(200);
    TS_ASSERT_EQUALS( task1->value(), 1 );
    TS_ASSERT_EQUALS( task2->value(), 1 );
    TS_ASSERT_EQUALS( task3->value(), 1 );

    Thread::sleep_ms(100);
    task3->kill();
    Thread::sleep_ms(100);
    TS_ASSERT_EQUALS( task1->value(), 1 );
    TS_ASSERT_EQUALS( task2->value(), 1 );
    TS_ASSERT_EQUALS( task3->value(), 3 );

    task2->kill();
    task1->kill();

    queue.join_all();
  }

  void test_threadpool_with_limited_threads()
  {
    boost::shared_ptr<TestTask> task1 (new TestTask);
    boost::shared_ptr<TestTask> task2 (new TestTask);
    boost::shared_ptr<TestTask> task3 (new TestTask);
    TS_ASSERT_EQUALS( task1->value(), 0 );
    TS_ASSERT_EQUALS( task2->value(), 0 );
    TS_ASSERT_EQUALS( task3->value(), 0 );

    set_debug_level(VerboseDebugMessage);
    FifoWorkQueue queue(2);
    queue.add_task(task1);
    queue.add_task(task2);
    queue.add_task(task3);

    // This is slightly sloppy -- we just wait a few hundred ms for
    // the thread to start and to set the flag.  This should really be
    // done with a condition variable.
    Thread::sleep_ms(200);
    TS_ASSERT_EQUALS( task1->value(), 1 );
    TS_ASSERT_EQUALS( task2->value(), 1 );

    Thread::sleep_ms(100);
    TS_ASSERT_EQUALS( task3->value(), 0 );
    task1->kill();
    Thread::sleep_ms(100);
    TS_ASSERT_EQUALS( task1->value(), 3 );
    TS_ASSERT_EQUALS( task2->value(), 1 );
    TS_ASSERT_EQUALS( task3->value(), 1 );

    task2->kill();
    task3->kill();

    queue.join_all();
  }


  void test_threadpool_spawing_new_tasks()
  {
    boost::shared_ptr<TestTask> task1 (new TestTask);
    boost::shared_ptr<TestTask> task2 (new TestTask);
    boost::shared_ptr<TestTask> task3 (new TestTask);
    TS_ASSERT_EQUALS( task1->value(), 0 );
    TS_ASSERT_EQUALS( task2->value(), 0 );
    TS_ASSERT_EQUALS( task3->value(), 0 );

    set_debug_level(VerboseDebugMessage);
    FifoWorkQueue queue(8);
    queue.add_task(task1);
    queue.add_task(task2);
    queue.add_task(task3);

    // Give the tasks a chance to start up...
    Thread::sleep_ms(100);

    set_debug_level(VerboseDebugMessage);

    TS_ASSERT_EQUALS( task1->value(), 1 );
    TS_ASSERT_EQUALS( task2->value(), 1 );

    Thread::sleep_ms(100);
    task1->kill();
    Thread::sleep_ms(100);
    TS_ASSERT_EQUALS( task1->value(), 3 );
    TS_ASSERT_EQUALS( task2->value(), 1 );
    TS_ASSERT_EQUALS( task3->value(), 1 );


    boost::shared_ptr<TestTask> task4 (new TestTask);
    boost::shared_ptr<TestTask> task5 (new TestTask);
    boost::shared_ptr<TestTask> task6 (new TestTask);
    TS_ASSERT_EQUALS( task4->value(), 0 );
    TS_ASSERT_EQUALS( task5->value(), 0 );
    TS_ASSERT_EQUALS( task6->value(), 0 );
    queue.add_task(task4);
    queue.add_task(task5);
    queue.add_task(task6);
    Thread::sleep_ms(100);
    TS_ASSERT_EQUALS( task4->value(), 1 );
    TS_ASSERT_EQUALS( task5->value(), 1 );
    TS_ASSERT_EQUALS( task6->value(), 1 );

    task2->kill();
    task3->kill();

    task4->kill();
    task5->kill();
    task6->kill();
    queue.join_all();
  }


  void test_ordered_workqueue()
  {
    boost::shared_ptr<TestTask> task1 (new TestTask);
    boost::shared_ptr<TestTask> task2 (new TestTask);
    boost::shared_ptr<TestTask> task3 (new TestTask);
    TS_ASSERT_EQUALS( task1->value(), 0 );
    TS_ASSERT_EQUALS( task2->value(), 0 );
    TS_ASSERT_EQUALS( task3->value(), 0 );

    set_debug_level(VerboseDebugMessage);
    OrderedWorkQueue queue(8);
    queue.add_task(task1, 2);
    queue.add_task(task2, 1);

    set_debug_level(VerboseDebugMessage);

    TS_ASSERT_EQUALS( task1->value(), 0 );
    TS_ASSERT_EQUALS( task2->value(), 0 );
    TS_ASSERT_EQUALS( task3->value(), 0 );

    queue.add_task(task3, 0);
    Thread::sleep_ms(100);
    TS_ASSERT_EQUALS( task1->value(), 1 );
    TS_ASSERT_EQUALS( task2->value(), 1 );
    TS_ASSERT_EQUALS( task3->value(), 1 );

    boost::shared_ptr<TestTask> task4 (new TestTask);
    boost::shared_ptr<TestTask> task5 (new TestTask);
    boost::shared_ptr<TestTask> task6 (new TestTask);
    TS_ASSERT_EQUALS( task4->value(), 0 );
    TS_ASSERT_EQUALS( task5->value(), 0 );
    TS_ASSERT_EQUALS( task6->value(), 0 );
    queue.add_task(task4, 3);
    queue.add_task(task5, 5);
    Thread::sleep_ms(100);
    TS_ASSERT_EQUALS( task4->value(), 1 );
    TS_ASSERT_EQUALS( task5->value(), 0 );
    TS_ASSERT_EQUALS( task6->value(), 0 );

    queue.add_task(task6, 4);
    Thread::sleep_ms(100);
    TS_ASSERT_EQUALS( task4->value(), 1 );
    TS_ASSERT_EQUALS( task5->value(), 1 );
    TS_ASSERT_EQUALS( task6->value(), 1 );

    task1->kill();
    task2->kill();
    task3->kill();
    task4->kill();
    task5->kill();
    task6->kill();

    queue.join_all();
  }

};
