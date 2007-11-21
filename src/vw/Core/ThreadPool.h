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

/// \file Core/ThreadsPool.h
/// 
/// Note: All tasks need to be of the same type, but you can have a
/// common abstract base class if you want.
///
#ifndef __VW_CORE_THREADPOOL_H__
#define __VW_CORE_THREADPOOL_H__

#include <vw/Core/Thread.h>
#include <vw/Core/Debugging.h>
#include <vw/Core/Exception.h>

namespace vw {
  // ----------------------  --------------  ---------------------------
  // ----------------------       Task       ---------------------------
  // ----------------------  --------------  ---------------------------
  struct Task {
    virtual ~Task() {}
    virtual void operator()() = 0;
  };
    
  template <class TaskT>
  class GenericTask : public Task {
    TaskT &m_task;
  public:
    GenericTask(TaskT &task) : m_task(task) {}
    virtual ~GenericTask() {}
    virtual void operator()() { m_task(); }
  };

  // ----------------------  --------------  ---------------------------
  // ----------------------  Task Generator  ---------------------------
  // ----------------------  --------------  ---------------------------

  struct TaskGenerator {
    virtual ~TaskGenerator() {}
    virtual bool idle() = 0;
    virtual boost::shared_ptr<Task> request_task() = 0;
  };

  class FifoTaskGenerator : public TaskGenerator {
    std::list<boost::shared_ptr<Task> > m_queued_tasks;
  public:

    int size() { return m_queued_tasks.size(); }
    
    // Add a task that is being tracked by a shared pointer.
    void add_task(boost::shared_ptr<Task> task) { m_queued_tasks.push_back(task); }

    // must keep the task object around until the TaskGenerator is
    // done using it.
    template <class TaskT>
    void add_task(TaskT &task) {
      boost::shared_ptr<Task> task_wrapper( new GenericTask<TaskT>(task) );
      add_task(task_wrapper);
    }

    virtual bool idle() { return !m_queued_tasks.size(); }
    virtual boost::shared_ptr<Task> request_task() { 
      if (m_queued_tasks.size() == 0) 
        vw_throw(LogicErr() << "FifoTaskGenerator: a task was requested when none were available.");
      boost::shared_ptr<Task> task = m_queued_tasks.front();
      m_queued_tasks.pop_front();
      return task;
    }
  };

  // ----------------------  --------------  ---------------------------
  // ----------------------   Thread Pool    ---------------------------
  // ----------------------  --------------  ---------------------------

  class ThreadPool {
    
    // The worker thread class is cognizant of both the thread and the
    // threadpool object.  When a worker thread finishes its task it
    // notifies the threadpool, which farms out the next task to the
    // worker thread.
    class WorkerThread {
      ThreadPool &m_pool;
      boost::shared_ptr<Task> m_task;
    public:
      WorkerThread(ThreadPool& pool, boost::shared_ptr<Task> task) : m_pool(pool), m_task(task) {}
      ~WorkerThread() {}
      void operator()() { 
        (*m_task)();  
        m_pool.task_complete(this); 
      }
    };

    typedef std::pair<boost::shared_ptr<WorkerThread>,
                      boost::shared_ptr<Thread> > threadpool_pair; 

    int m_active_workers;
    int m_max_workers;
    Mutex m_mutex;
    std::list<threadpool_pair> m_running_threads;
    boost::shared_ptr<TaskGenerator> m_taskgen;
    Condition m_joined_event;

    // This is called whenever a worker thread finishes its task. If
    // there are more tasks available, the worker is given more work.
    // Otherwise, the worker terminates.
    void task_complete(WorkerThread* worker_ptr) {
      {
        Mutex::Lock lock(m_mutex);
        
        m_active_workers--;
        vw_out(DebugMessage) << "ThreadPool: terminated worker thread.  [ " << m_active_workers << " / " << m_max_workers << " now active ].\n"; 
      }
      
      this->check_for_tasks();
      
      {
        Mutex::Lock lock(m_mutex);        

        // Erase the worker thread from the list of active threads
        bool erased = false;
        for (std::list<threadpool_pair>::iterator it = m_running_threads.begin(); it != m_running_threads.end(); ++it) {
          if (&( *(it->first) ) == worker_ptr) {
            m_running_threads.erase(it);
            erased = true;
            break;
          }
        }
        // If execution reaches this point, something went wrong during
        // thread cleanup.  Sound the alarm!
        if (!erased)
          vw_throw(LogicErr() << "ThreadPool: tried to clean up a non-existent thread.  This is bad news!");
      }
      
      m_joined_event.notify_all();
    }

  public:
    
    // Constructor & Destructor
    ThreadPool(boost::shared_ptr<TaskGenerator> taskgen, int num_threads = Thread::default_num_threads()) : 
      m_max_workers(num_threads), m_active_workers(0), m_taskgen(taskgen) {
      this->check_for_tasks();
    }
    ~ThreadPool() { this->join_all(); }

    void set_task_generator(boost::shared_ptr<TaskGenerator> taskgen) {
      Mutex::Lock lock(m_mutex);
      m_taskgen = taskgen;
      this->check_for_tasks();
    }

    void check_for_tasks() {
      Mutex::Lock lock(m_mutex);

      // While there are available threads, farm out the tasks from
      // the task generator
      while (m_taskgen && m_active_workers < m_max_workers && !(m_taskgen->idle()) ) {
        boost::shared_ptr<WorkerThread> next_worker( new WorkerThread(*this, m_taskgen->request_task()) );
        boost::shared_ptr<Thread> thread(new Thread(next_worker));
        m_running_threads.push_back(threadpool_pair(next_worker, thread));
        m_active_workers++;
        vw_out(DebugMessage) << "ThreadPool: created worker thread.  [ " << m_active_workers << " / " << m_max_workers << " now active. ]\n"; 
      } 
    }

    /// Return the max number threads that can run concurrently at any
    /// given time using this threadpool.
    int max_threads() { 
      Mutex::Lock lock(m_mutex);
      return m_max_workers;
    }

    /// Return the max number threads that can run concurrently at any
    /// given time using this threadpool.
    int active_threads() { 
      Mutex::Lock lock(m_mutex);
      return m_active_workers;
    }

    void clear() { 
      Mutex::Lock lock(m_mutex);
      m_taskgen.reset();
    }

    // Join all threads and clear all pending tasks from the list.
    void join_all() {
      bool finished = false;
      
      // Wait for the threads to clean up the threadpool state and exit.
      while(!finished) {
        Mutex::Lock lock(m_mutex);
        if (m_running_threads.size() != 0) 
          m_joined_event.wait(lock);
        else 
          finished = true;
      }
    }
    
  };

} // namespace vw

#endif // __VW_CORE_THREADPOOL_H__
