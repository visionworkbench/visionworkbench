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


/// \file Core/Thread.h
///
/// A very simple, abstracted threading system.  Currently this is
/// simply built on top of Boost threads, but its abstracted so that
/// it's easy to port the Vision Workbench to environments where Boost
/// threads do not work or otherwise cannot be used, and because the
/// Boost threads design itself is expected to be overhauled at some
/// point.
///
/// To create a thread, pass in any object that has the operator()
/// defined.
///
/// The rest of the interface is straightforward.  These are
/// the things you will need to implement if you want to use a
/// different thread library.  The interface consists of:
///
/// * A Thread class, whose constructor takes a shared pointer to a
///   Task.  The constructor should invoke the Task's operator()
///   function in a new thread, and the destructor should join the
///   child thread.  The static yield() function should yield the
///   active thread, and the static sleep_ms() function should sleep
///   the active thread for the given number of milliseconds.  The
///   Thread class may be non-copyable.
///
/// * A Mutex class, implementing a simple mutual exclusion lock.
///   This should have no methods other than the default constructor
///   and destructor.  All locking and unlocking takes place through
///   a nested Mutex::Lock class, whose constructor takes a Mutex&
///   to lock and whose destructor unlocks it.  Both the Mutex and
///   Lock classes may be non-copyable.
///

#ifndef __VW_CORE_THREAD_H__
#define __VW_CORE_THREAD_H__

#include <vw/Core/FundamentalTypes.h>

#include <boost/thread/locks.hpp>
#include <boost/thread/recursive_mutex.hpp>
#include <boost/thread/shared_mutex.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/xtime.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/noncopyable.hpp>
#include <boost/version.hpp>

namespace vw {

  // --------------------------------------------------------------
  //                            MUTEX
  // --------------------------------------------------------------

  /// A wrapper around some boost::mutex classes.
  /// - See docs: http://www.boost.org/doc/libs/1_59_0/doc/html/thread/synchronization.html
  /// - There are three locks being managed here:
  ///   - lock         = Normal, exclusive access lock.
  ///   - lock_shared  = Non-exclusive lock.
  ///   - lock_upgrade = Similar to lock_shared, but with the ability to upgrade to lock.
  ///                    Having this intermediate lock step is necessary to secure exclusive 
  ///                    access (lock) in a timely manner.
  class Mutex : private boost::shared_mutex {

    friend class WriteLock;
    friend class ReadLock;

  public:
    inline Mutex() {}

    /// Block until you own the requested lock.
    void lock        () { boost::shared_mutex::lock();         }
    void lock_shared () { boost::shared_mutex::lock_shared();  }
    void lock_upgrade() { boost::shared_mutex::lock_upgrade(); }
    
    /// Non-blocking attempt to obtain ownership of a lock.
    bool try_lock        () { return boost::shared_mutex::try_lock();         }
    bool try_lock_shared () { return boost::shared_mutex::try_lock_shared();  }
    bool try_lock_upgrade() { return boost::shared_mutex::try_lock_upgrade(); }
    
    /// Release ownership of a lock.
    void unlock        () { boost::shared_mutex::unlock();         }
    void unlock_shared () { boost::shared_mutex::unlock_shared();  }
    void unlock_upgrade() { boost::shared_mutex::unlock_upgrade(); }
    
    /// Atomic operations to switch ownership from one type of lock to another.
    void unlock_upgrade_and_lock       () { boost::shared_mutex::unlock_upgrade_and_lock(); }
    void unlock_and_lock_upgrade       () { boost::shared_mutex::unlock_and_lock_upgrade(); }
    void unlock_upgrade_and_lock_shared() { boost::shared_mutex::unlock_upgrade_and_lock_shared(); }

    /// A unique scoped lock class, used to lock and unlock a Mutex (only one can own at a time).
    /// - Automatically locks the mutex when constructed, unlocks on destruction.
    class WriteLock : private boost::unique_lock<Mutex>, private boost::noncopyable {
    public:
      inline WriteLock( Mutex &mutex ) : boost::unique_lock<Mutex>( mutex ) {}
      void lock()          { boost::unique_lock<Mutex>::lock(); }
      void unlock()        { boost::unique_lock<Mutex>::unlock(); }
    };
    
    /// A scoped lock class, used to lock and unlock a Mutex (allows shared ownership).
    /// - Automatically locks the mutex when constructed, unlocks on destruction.
    class ReadLock : private boost::shared_lock<Mutex>, private boost::noncopyable {
    public:
      inline ReadLock( Mutex &mutex ) : boost::shared_lock<Mutex>( mutex ) {}
      void lock()          { boost::shared_lock<Mutex>::lock(); }
      void unlock()        { boost::shared_lock<Mutex>::unlock(); }
    };

    typedef class WriteLock Lock; /// By default, use the non-shared lock type.
  };// End class Mutex

  /// A wrapper around the boost::recursive_mutex class.
  class RecursiveMutex : private boost::recursive_mutex {

    friend class Lock;

  public:
    inline RecursiveMutex() {}

    void lock()   { boost::recursive_mutex::lock(); }
    void unlock() { boost::recursive_mutex::unlock(); }

    // A unique scoped lock class, used to lock and unlock a Mutex (only one can own at a time).
    class Lock : private boost::unique_lock<RecursiveMutex>, private boost::noncopyable {
    public:
      inline Lock( RecursiveMutex &mutex ) : boost::unique_lock<RecursiveMutex>( mutex ) {}
      void lock()   { boost::unique_lock<RecursiveMutex>::lock(); }
      void unlock() { boost::unique_lock<RecursiveMutex>::unlock(); }
    }; // End class Lock
    
  }; // End class RecursiveLock


  // --------------------------------------------------------------
  //                            THREAD
  // --------------------------------------------------------------

  /// A thread class, that runs a "Task", which is an object or
  /// function that has the operator() defined.  When the Thread object
  /// is destroyed it will join the child thread if it has not already terminated.
  class Thread : private boost::noncopyable {

    boost::thread m_thread;

    // For some reason, the boost thread library makes a copy of the
    // Task object before handing it off to the thread.  This is
    // annoying, because it means that the parent thread no longer has
    // direct access to the child thread's instance of the task object.
    // This helper allows the parent thread to retain direct access to
    // the child instance of the task.
    template <class TaskT>
    class TaskHelper {
      boost::shared_ptr<TaskT> m_task;
    public:
      TaskHelper(boost::shared_ptr<TaskT> task) : m_task(task) {}
      void operator() () { (*m_task)(); }
    };

  public:

    /// This variant of the constructor takes a Task that is copy
    /// constructable.  The thread made a copy of the task, and this
    /// instance is no longer directly accessibly from the parent thread.
    template<class TaskT>
    inline Thread( TaskT task ) : m_thread( task ) {}

    /// This variant of the constructor takes a shared pointer to a task.
    /// The thread makes a copy of the shared pointer task, allowing
    /// the parent to still access the task instance that is running in the thread.
    template<class TaskT>
    inline Thread( boost::shared_ptr<TaskT> task ) : m_thread( TaskHelper<TaskT>(task) ) {}

    /// Destroys the thread. User is expected to either call join() themselves,
    /// or let the thread run free. We don't call join() here because most users
    /// will call join(), and a second call to join() is undefined.
    inline ~Thread() { }

    /// The current thread of execution blocks until this thread
    /// finishes execution of the task and all resources are reclaimed.
    inline void join() { m_thread.join(); }


    // --------------
    // STATIC METHODS
    // --------------

    /// Returns a unique ID for the current thread.  The ID for a
    /// thread is not determined until the thread calls the id()
    /// function for the first time, so there is no guarantee that IDs
    /// will be assigned in the same order that threads are created.
    static vw::uint64 id();

    /// Cause the current thread to yield the remainder of its
    /// execution time to the kernel's scheduler.
    static inline void yield() { boost::thread::yield(); }

    /// Cause the current thread to sleep for a specified amount of
    /// time.  The thread will not be scheduled to run at all for the
    /// duration, so machine resources are free for other
    /// threads/processes.
    static inline void sleep_ms( uint32 milliseconds ) {
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
      boost::thread::sleep(xt);
    }
  }; // End class Thread


} // namespace vw

#endif // __VW_CORE_THREAD_H__
