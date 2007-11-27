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
    FifoWorkQueue *taskgen = new FifoWorkQueue;
    taskgen->add_task(task1);
    taskgen->add_task(task2);
    taskgen->add_task(task3);

    boost::shared_ptr<WorkQueue> taskgen_ptr(taskgen);
    ThreadPool pool(taskgen_ptr, 4);
    
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
    
    pool.join_all();
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
    FifoWorkQueue *taskgen = new FifoWorkQueue;
    taskgen->add_task(task1);
    taskgen->add_task(task2);
    taskgen->add_task(task3);

    set_debug_level(VerboseDebugMessage);
    boost::shared_ptr<WorkQueue> taskgen_ptr(taskgen);
    ThreadPool pool(taskgen_ptr, 2);
    
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
    
    pool.join_all();
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
    FifoWorkQueue *taskgen = new FifoWorkQueue;
    taskgen->add_task(task1);
    taskgen->add_task(task2);
    taskgen->add_task(task3);

    set_debug_level(VerboseDebugMessage);
    boost::shared_ptr<WorkQueue> taskgen_ptr(taskgen);
    ThreadPool pool(taskgen_ptr, 8);
    
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
    taskgen->add_task(task4);
    taskgen->add_task(task5);
    taskgen->add_task(task6);
    Thread::sleep_ms(100);
    TS_ASSERT_EQUALS( task4->value(), 1 );
    TS_ASSERT_EQUALS( task5->value(), 1 );
    TS_ASSERT_EQUALS( task6->value(), 1 );

    task2->kill();
    task3->kill();

    task4->kill();
    task5->kill();
    task6->kill();
    pool.join_all();
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
    OrderedWorkQueue *taskgen = new OrderedWorkQueue;
    taskgen->add_task(task1, 2);
    taskgen->add_task(task2, 1);
    
    set_debug_level(VerboseDebugMessage);
    boost::shared_ptr<WorkQueue> taskgen_ptr(taskgen);
    ThreadPool pool(taskgen_ptr, 8);
    
    TS_ASSERT_EQUALS( task1->value(), 0 );
    TS_ASSERT_EQUALS( task2->value(), 0 );
    TS_ASSERT_EQUALS( task3->value(), 0 );
    
    taskgen->add_task(task3, 0);
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
    taskgen->add_task(task4, 3);
    taskgen->add_task(task5, 5);
    Thread::sleep_ms(100);
    TS_ASSERT_EQUALS( task4->value(), 1 );
    TS_ASSERT_EQUALS( task5->value(), 0 );
    TS_ASSERT_EQUALS( task6->value(), 0 );

    taskgen->add_task(task6, 4);
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

    pool.join_all();
  }

};
