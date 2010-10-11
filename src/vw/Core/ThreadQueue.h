// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Core/Queue.h
///
/// A thread-safe queue. Uses a condition variable to avoid busy-waiting.
///

#ifndef __VW_CORE_QUEUE_H__
#define __VW_CORE_QUEUE_H__

#include <boost/bind.hpp>

#include <queue>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Thread.h>

namespace vw {

template<typename T>
class ThreadQueue : private boost::noncopyable {
  private:
    std::queue<T> m_queue;
    mutable Mutex m_mutex;
    Condition m_cond;
  public:
    void push(T const& data) {
      {
        Mutex::Lock lock(m_mutex);
        m_queue.push(data);
      }

      // notify AFTER we release the lock, so we don't wake up a thread just to
      // have it run into a lock
      m_cond.notify_one();
    }

    bool empty() const {
      Mutex::Lock lock(m_mutex);
      return m_queue.empty();
    }

    // Try to pop something off, return indicates success
    bool try_pop(T& data) {
      Mutex::Lock lock(m_mutex);
      if(m_queue.empty()) {
        return false;
      }

      data = m_queue.front();
      m_queue.pop();
      return true;
    }

    // Returns the number of messages waiting in the queue.
    int size() const {
      Mutex::Lock lock(m_mutex);
      return m_queue.size();
    }


    // Wait for data forever
    void wait_pop(T& data) {
      Mutex::Lock lock(m_mutex);
      while(m_queue.empty()) {
        m_cond.wait(lock);
      }

      data = m_queue.front();
      m_queue.pop();
    }

    // Wait for data with a timeout (in ms)
    bool timed_wait_pop(T& data, unsigned long duration) {
      Mutex::Lock lock(m_mutex);

      // Use the predicate to wait until the queue is non-empty (to catch spurious wakeups)
      if(!m_cond.timed_wait(lock,duration,!boost::bind(&std::queue<T>::empty, boost::cref(m_queue))))
        return false;

      data = m_queue.front();
      m_queue.pop();
      return true;
    }

    void flush() {
      Mutex::Lock lock(m_mutex);
      while (!m_queue.empty()) {
        m_queue.pop();
      }
    }
};

} // namespace vw


#endif
