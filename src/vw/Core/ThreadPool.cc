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

#include <vw/vw_config.h>
#include <vw/Core/Exception.h>
#include <vw/Core/Log.h>
#include <vw/Core/ThreadPool.h>

#include <ostream>

using namespace vw;

//----------------------------------------------------
// Task


bool Task::is_finished() {
  Mutex::Lock lock(m_task_mutex);
  return m_finished;
}
void Task::join() {
  while (1) {
    Mutex::Lock lock(m_task_mutex);
    if (m_finished) return;
    m_finished_event.wait(lock);
  }
}
void Task::signal_finished() {
  Mutex::Lock lock(m_task_mutex);
  m_finished = true;
  m_finished_event.notify_all();
}

//----------------------------------------------------
// WorkQueue


WorkQueue::WorkerThread::WorkerThread(WorkQueue& queue, boost::shared_ptr<Task> initial_task,
                                      int thread_id, bool &should_die) :
  m_queue(queue), m_task(initial_task), m_thread_id(thread_id), m_should_die(should_die) {}


void WorkQueue::WorkerThread::operator()() {
  do {
    VW_OUT(DebugMessage, "thread") << "ThreadPool: running worker thread "
                                   << m_thread_id << "\n";
    // Run the task and then signal that it is finished
    (*m_task)();
    m_task->signal_finished();

    {
      // We lock m_queue_mutex to prevent WorkQueue::notify() from running
      // until we either sucessfully have grabbed the next task, or we have
      // completely terminated the worker.
      Mutex::Lock lock(m_queue.m_queue_mutex);
      m_task = m_queue.get_next_task();

      if (!m_task) // No more tasks, notify parent queue that we are finished.
        m_queue.worker_thread_complete(m_thread_id);
    }
  } while ( m_task && !m_should_die ); // Quit if no task or when instructed
}



void WorkQueue::worker_thread_complete(int worker_id) {
  m_active_workers--;
  VW_OUT(DebugMessage, "thread") << "ThreadPool: terminating worker thread " << worker_id << ".  [ " << m_active_workers << " / " << m_max_workers << " now active ]\n";

  // Erase the worker thread from the list of active threads
  VW_ASSERT(worker_id >= 0 && worker_id < int(m_running_threads.size()),
            LogicErr() << "WorkQueue: request to terminate thread " << worker_id << ", which does not exist.");
  m_available_thread_ids.push_back(worker_id);

  // Notify any threads that are waiting for the join event.
  m_joined_event.notify_all();
}

WorkQueue::WorkQueue(int num_threads )
  : m_active_workers(0), m_max_workers(num_threads), m_should_die(false) {
  m_running_threads.resize(num_threads);
  for (int i = 0; i < num_threads; ++i)
    m_available_thread_ids.push_back(i);
}
WorkQueue::~WorkQueue() { this->join_all(); }

void WorkQueue::notify() {
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
    VW_OUT(DebugMessage, "thread") << "ThreadPool: creating worker thread " << next_available_thread_id << ".  [ " << m_active_workers << " / " << m_max_workers << " now active ]\n";
  }
}

int WorkQueue::max_threads() {
  // TODO: Can this be atomic?
  Mutex::Lock lock(m_queue_mutex);
  return m_max_workers;
}

int WorkQueue::active_threads() {
  // TODO: Can this be atomic?
  Mutex::Lock lock(m_queue_mutex);
  return m_active_workers;
}

// Join all currently running threads and wait for the task pool to be empty.
void WorkQueue::join_all() {
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

void WorkQueue::kill_and_join() {
  m_should_die = true;
  this->join_all();
}


//----------------------------------------------------
// FifoWorkQueue

FifoWorkQueue::FifoWorkQueue(int num_threads) : WorkQueue(num_threads) {}

size_t FifoWorkQueue::size() {
  Mutex::Lock lock(m_mutex);
  return m_queued_tasks.size();
}

// Add a task that is being tracked by a shared pointer.
void FifoWorkQueue::add_task(boost::shared_ptr<Task> task) {
  {
    Mutex::Lock lock(m_mutex);
    m_queued_tasks.push_back(task);
  }
  this->notify();
}

boost::shared_ptr<Task> FifoWorkQueue::get_next_task() {
  Mutex::Lock lock(m_mutex);
  if (m_queued_tasks.empty())
    return boost::shared_ptr<Task>();

  boost::shared_ptr<Task> task = m_queued_tasks.front();
  m_queued_tasks.pop_front();
  return task;
}

//----------------------------------------------------
// OrderedWorkQueue

OrderedWorkQueue::OrderedWorkQueue(int num_threads) : WorkQueue(num_threads) {
  m_next_index = 0;
}

size_t OrderedWorkQueue::size() {
  Mutex::Lock lock(m_mutex);
  return m_queued_tasks.size();
}

// Add a task that is being tracked by a shared pointer.
void OrderedWorkQueue::add_task(boost::shared_ptr<Task> task, int index) {
  {
    Mutex::Lock lock(m_mutex);
    m_queued_tasks[index] = task;
  }
  this->notify();
}

boost::shared_ptr<Task> OrderedWorkQueue::get_next_task() {
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
