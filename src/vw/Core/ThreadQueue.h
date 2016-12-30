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


/// \file vw/Core/ThreadQueue.h
///

#ifndef __VW_CORE_QUEUE_H__
#define __VW_CORE_QUEUE_H__

#include <boost/bind.hpp>

#include <queue>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Thread.h>
#include <vw/Core/Condition.h>

namespace vw {

/// A thread-safe queue. Uses a condition variable to avoid busy-waiting.
/// - Use this class to safely pass messages and data between threads.
template<typename T>
class ThreadQueue : private boost::noncopyable {
  private:
    std::queue<T> m_queue;
    mutable Mutex m_mutex;
    Condition m_cond;
  public:
    
    /// Push an object on to the queue.
    void push(T const& data) {
      {
        Mutex::Lock lock(m_mutex);
        m_queue.push(data);
      }

      // notify AFTER we release the lock, so we don't wake up a thread just to
      // have it run into a lock
      m_cond.notify_one();
    }

    /// Returns true if the queue is empty.
    bool empty() const {
      Mutex::Lock lock(m_mutex);
      return m_queue.empty();
    }

    /// Try to pop something off, return indicates success.
    /// - Failure indicates that the queue is empty.
    bool try_pop(T& data) {
      Mutex::Lock lock(m_mutex);
      if(m_queue.empty()) {
        return false;
      }

      data = m_queue.front();
      m_queue.pop();
      return true;
    }

    /// Returns the number of messages waiting in the queue.
    size_t size() const {
      Mutex::Lock lock(m_mutex);
      return m_queue.size();
    }


    /// Wait forever until data is available, then get it.
    void wait_pop(T& data) {
      Mutex::Lock lock(m_mutex);
      while(m_queue.empty()) {
        m_cond.wait(lock);
      }

      data = m_queue.front();
      m_queue.pop();
    }

    /// Wait for data with a timeout (in ms).
    bool timed_wait_pop(T& data, unsigned long duration) {
      Mutex::Lock lock(m_mutex);

      // Use the predicate to wait until the queue is non-empty (to catch spurious wakeups)
      if(!m_cond.timed_wait(lock,duration,!boost::bind(&std::queue<T>::empty, boost::cref(m_queue))))
        return false;

      data = m_queue.front();
      m_queue.pop();
      return true;
    }

    /// Clear the contents of the queue.
    void flush() {
      Mutex::Lock lock(m_mutex);
      while (!m_queue.empty()) {
        m_queue.pop();
      }
    }
}; // End class ThreadQueue

} // namespace vw


#endif
