// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Core/ThreadsPool.h
///
/// Note: All tasks need to be of the same type, but you can have a
/// common abstract base class if you want.
///
#ifndef __VW_CORE_THREADPOOL_H__
#define __VW_CORE_THREADPOOL_H__

#include <vector>
#include <list>

#include <vw/Core/Settings.h>
#include <vw/Core/Thread.h>
#include <vw/Core/Log.h>
#include <vw/Core/Exception.h>

// STL
#include <map>

namespace vw {
  // ----------------------  --------------  ---------------------------
  // ----------------------       Task       ---------------------------
  // ----------------------  --------------  ---------------------------
  class Task {
    Mutex m_task_mutex;
    Condition m_finished_event;
    volatile bool m_finished;

  public:
    Task() : m_finished(false) {}
    virtual ~Task() {}
    virtual void operator()() = 0;

    // Keep track of whether task is finished
    // WorkQueue is responsible for calling "signal_finished" after task is finished
    bool is_finished() {
      Mutex::Lock lock(m_task_mutex);
      return m_finished;
    }
    void join() {
      while (1) {
        Mutex::Lock lock(m_task_mutex);
        if (m_finished) return;
        m_finished_event.wait(lock);
      }
    }
    void signal_finished() {
      Mutex::Lock lock(m_task_mutex);
      m_finished = true;
      m_finished_event.notify_all();
    }
  };

  // ----------------------  --------------  ---------------------------
  // ----------------------  Task Generator  ---------------------------
  // ----------------------  --------------  ---------------------------

  // Work Queue Base Class
  class WorkQueue {

    // The worker thread class is the Task object that is spun out to
    // do the actual work of the WorkQueue.  When a worker thread
    // finishes its task it notifies the threadpool, which farms out
    // the next task to the worker thread.
    class WorkerThread {
      WorkQueue &m_queue;
      boost::shared_ptr<Task> m_task;
      int m_thread_id;
      bool &m_should_die;
    public:
      WorkerThread(WorkQueue& queue, boost::shared_ptr<Task> initial_task,
                   int thread_id, bool &should_die) :
        m_queue(queue), m_task(initial_task), m_thread_id(thread_id), m_should_die(should_die) {}
      ~WorkerThread() {}
      void operator()() {
        do {
          vw_out(DebugMessage, "thread") << "ThreadPool: running worker thread "
                                         << m_thread_id << "\n";
          (*m_task)();
          m_task->signal_finished();

          {
            // We lock m_queue_mutex to prevent WorkQueue::notify() from running
            // until we either sucessfully have grabbed the next task, or we have
            // completely terminated the worker.
            Mutex::Lock lock(m_queue.m_queue_mutex);
            m_task = m_queue.get_next_task();

            if (!m_task)
              m_queue.worker_thread_complete(m_thread_id);
          }
        } while ( m_task && !m_should_die );
      }
    };

    int m_active_workers, m_max_workers;
    Mutex m_queue_mutex;
    std::vector<boost::shared_ptr<Thread> > m_running_threads;
    std::list<int> m_available_thread_ids;
    Condition m_joined_event;
    bool m_should_die;

    // This is called whenever a worker thread finishes its task. If
    // there are more tasks available, the worker is given more work.
    // Otherwise, the worker terminates.
    //
    // *************************************************************
    // IMPORANT NOTE: The worker_thread_complete() method is called
    // from within child thread so that we can clean up the list of
    // available threads.  As such, one must be very careful about
    // what is done in this method, because for the duration of this
    // method, the child thread has access to it's own shared pointer.
    // This method has been carefully written so as not to
    // accidentally allow the thread to delete it's own shared
    // pointer, which could cause the thread to call its own join()
    // method.
    // *************************************************************
    void worker_thread_complete(int worker_id) {
      m_active_workers--;
      vw_out(DebugMessage, "thread") << "ThreadPool: terminating worker thread " << worker_id << ".  [ " << m_active_workers << " / " << m_max_workers << " now active ]\n";

      // Erase the worker thread from the list of active threads
      VW_ASSERT(worker_id >= 0 && worker_id < int(m_running_threads.size()),
                LogicErr() << "WorkQueue: request to terminate thread " << worker_id << ", which does not exist.");
      m_available_thread_ids.push_back(worker_id);

      // Notify any threads that are waiting for the join event.
      m_joined_event.notify_all();
    }

  public:
    WorkQueue(int num_threads = vw_settings().default_num_threads() )
      : m_active_workers(0), m_max_workers(num_threads), m_should_die(false) {
      m_running_threads.resize(num_threads);
      for (int i = 0; i < num_threads; ++i)
        m_available_thread_ids.push_back(i);
    }
    virtual ~WorkQueue() { this->join_all(); }

    /// Return a shared pointer to the next task.  If no tasks are
    /// available, return an empty shared pointer.
    virtual boost::shared_ptr<Task> get_next_task() = 0;

    // Notify can be called by a child class that inherits from
    // WorkQueue.  A call to notify will cause the WorkQueue to
    // re-examine the list of tasks it has available for execution.
    // If there are any idle slots for worker threads, it will spin
    // off WorkerThreads to execute these tasks.
    void notify() {
      Mutex::Lock lock(m_queue_mutex);

      // While there are available threads, farm out the tasks from
      // the task generator
      boost::shared_ptr<Task> task;
      while ( !m_available_thread_ids.empty() &&
              (task = this->get_next_task()) ) {
        int next_available_thread_id = m_available_thread_ids.front();
        m_available_thread_ids.pop_front();

        boost::shared_ptr<WorkerThread> next_worker( new WorkerThread(*this, task,
                                                                      next_available_thread_id,
                                                                      m_should_die) );
        boost::shared_ptr<Thread> thread(new Thread(next_worker));
        m_running_threads[next_available_thread_id] = thread;
        m_active_workers++;
        vw_out(DebugMessage, "thread") << "ThreadPool: creating worker thread " << next_available_thread_id << ".  [ " << m_active_workers << " / " << m_max_workers << " now active ]\n";
      }
    }

    /// Return the max number threads that can run concurrently at any
    /// given time using this threadpool.
    int max_threads() {
      Mutex::Lock lock(m_queue_mutex);
      return m_max_workers;
    }

    /// Return the max number threads that can run concurrently at any
    /// given time using this threadpool.
    int active_threads() {
      Mutex::Lock lock(m_queue_mutex);
      return m_active_workers;
    }

    // Join all currently running threads and wait for the task pool to be empty.
    void join_all() {
      bool finished = false;

      // Wait for the threads to clean up the threadpool state and exit.
      while(!finished) {
        Mutex::Lock lock(m_queue_mutex);
        if (m_active_workers != 0) {
          m_joined_event.wait(lock);
        } else {
          finished = true;
        }
      }
    }

    void kill_and_join() {
      m_should_die = true;
      this->join_all();
    }

  };



  /// A simple, first-in, first-out work queue.
  class FifoWorkQueue : public WorkQueue {
    std::list<boost::shared_ptr<Task> > m_queued_tasks;
    Mutex m_mutex;
  public:

    FifoWorkQueue(int num_threads = vw_settings().default_num_threads()) : WorkQueue(num_threads) {}

    size_t size() {
      Mutex::Lock lock(m_mutex);
      return m_queued_tasks.size();
    }

    // Add a task that is being tracked by a shared pointer.
    void add_task(boost::shared_ptr<Task> task) {
      {
        Mutex::Lock lock(m_mutex);
        m_queued_tasks.push_back(task);
      }
      this->notify();
    }

    virtual boost::shared_ptr<Task> get_next_task() {
      Mutex::Lock lock(m_mutex);
      if (m_queued_tasks.empty())
        return boost::shared_ptr<Task>();

      boost::shared_ptr<Task> task = m_queued_tasks.front();
      m_queued_tasks.pop_front();
      return task;
    }
  };

  /// A simple ordered work queue.  Tasks are each given an "index"
  /// and they are processed in order starting with the task at index
  /// 0.  The idle() method returns true unless the task with the next
  /// expected index is present in the work queue.
  class OrderedWorkQueue : public WorkQueue {
    std::map<int, boost::shared_ptr<Task> > m_queued_tasks;
    int m_next_index;
    Mutex m_mutex;
  public:

    OrderedWorkQueue(int num_threads = vw_settings().default_num_threads()) : WorkQueue(num_threads) {
      m_next_index = 0;
    }

    size_t size() {
      Mutex::Lock lock(m_mutex);
      return m_queued_tasks.size();
    }

    // Add a task that is being tracked by a shared pointer.
    void add_task(boost::shared_ptr<Task> task, int index) {
      {
        Mutex::Lock lock(m_mutex);
        m_queued_tasks[index] = task;
      }
      this->notify();
    }

    virtual boost::shared_ptr<Task> get_next_task() {
      Mutex::Lock lock(m_mutex);

      // If there are no tasks available, we return the NULL task.
      if (m_queued_tasks.empty())
        return boost::shared_ptr<Task>();

      std::map<int, boost::shared_ptr<Task> >::iterator iter = m_queued_tasks.begin();

      // If the next task does not have the expected index, we
      // return the NULL task.
      if ((*iter).first != m_next_index)
        return boost::shared_ptr<Task>();


      boost::shared_ptr<Task> task = (*iter).second;
      m_queued_tasks.erase(m_queued_tasks.begin());
      m_next_index++;
      return task;
    }
  };

} // namespace vw

#endif // __VW_CORE_THREADPOOL_H__
