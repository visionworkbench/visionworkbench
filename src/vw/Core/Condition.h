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

#ifndef __VW_CORE_CONDITION_H__
#define __VW_CORE_CONDITION_H__

#include <vw/Core/FundamentalTypes.h>

#include <boost/thread/condition.hpp>

namespace vw {

  // --------------------------------------------------------------
  //                            CONDITION
  // --------------------------------------------------------------

  /// Wrapper around the boost::condition class
  class Condition : private boost::condition,
                    private boost::noncopyable
  {

  public:
    // construct/copy/destruct
    inline Condition() : boost::condition() {}

    /// Unblock (notify) a single thread blocking waiting for the condition.
    void notify_one() {
      boost::condition::notify_one();
    }

    /// Unblocks (notifys) all threads blocking waiting for the condition.
    void notify_all() {
      boost::condition::notify_all();
    }

    /// Unlock the mutex in the lock, then wait until notified.
    /// - When notified and woken up, the mutex is re-locked.
    template<typename LockT> void wait(LockT &lock) {
      boost::condition::wait(lock);
    }

    /// Identical to:  while(!pred()){wait(lock);}
    template<typename LockT, typename Pred>
    void wait(LockT &lock, Pred pred) {
      boost::condition::wait(lock,pred);
    }

    /// As wait(), but also wake up if the time expires.
    template<typename LockT>
    bool timed_wait(LockT &lock, unsigned long milliseconds) {
      boost::xtime xt;
#if BOOST_VERSION >= 105000
      boost::xtime_get(&xt, boost::TIME_UTC_);
#else
      boost::xtime_get(&xt, boost::TIME_UTC);
#endif
      while (milliseconds >= 1000) {
        xt.sec++;
        milliseconds -= 1000;
      }
      xt.nsec += static_cast<uint32>(1e6) * milliseconds;
      return boost::condition::timed_wait(lock, xt);
    }

    /// ?
    template<typename LockT, typename Pred>
    bool timed_wait(LockT &lock, unsigned long milliseconds, Pred pred) {
      boost::xtime xt;
#if BOOST_VERSION >= 105000
      boost::xtime_get(&xt, boost::TIME_UTC_);
#else
      boost::xtime_get(&xt, boost::TIME_UTC);
#endif
      while (milliseconds >= 1000) {
        xt.sec++;
        milliseconds -= 1000;
      }
      xt.nsec += static_cast<uint32>(1e6) * milliseconds;
      return boost::condition::timed_wait(lock, xt, pred);
    }
  };

}

#endif//__VW_CORE_CONDITION_H__
