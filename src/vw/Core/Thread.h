// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
/// * A RunOnce class, implementing run-once semantics.  This must be
///   a POD class with an initializer macro VW_RUNONCE_INIT.  It must
///   implement a single method, run( void (*func)() ), which runs the
///   given function exactly once no matter how many times it is
///   called.  (This behavior is only defined for RunOnce objects
///   that are statically allocated at global or namespace scope and
///   statically initialized to VW_RUNONCE_INIT.)

#ifndef __VW_CORE_THREAD_H__
#define __VW_CORE_THREAD_H__

#include <vw/Core/FundamentalTypes.h>

#include <boost/thread.hpp>
#include <boost/thread/condition.hpp>
#include <boost/thread/xtime.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>

namespace vw {

  // --------------------------------------------------------------
  //                            RUNONCE
  // --------------------------------------------------------------

#define VW_RUNONCE_INIT { BOOST_ONCE_INIT }

  // A special POD class to enable safe library initialization.  You
  // should only define these objects at global or namespace scope,
  // and statically initialize them to VW_RUNONCE_INIT.
  struct RunOnce {
    boost::once_flag m_flag;

    inline void run( void (*func)() ) {
      boost::call_once( func, m_flag );
    }
  };

  // --------------------------------------------------------------
  //                            MUTEX
  // --------------------------------------------------------------

  // A simple mutual exclusion class.
  class Mutex : private boost::mutex {

    friend class Lock;

  public:
    inline Mutex() {}

    void lock()   { boost::mutex::lock(); }
    void unlock() { boost::mutex::unlock(); }

    // A scoped lock class, used to lock and unlock a Mutex.
    class Lock : private boost::unique_lock<Mutex>,
                 private boost::noncopyable {

    public:
      inline Lock( Mutex &mutex ) : boost::unique_lock<Mutex>( mutex ) {}
      void lock()   { boost::unique_lock<Mutex>::lock(); }
      void unlock() { boost::unique_lock<Mutex>::unlock(); }
    };
  };

  // A simple mutual exclusion class.
  class RecursiveMutex : private boost::recursive_mutex {

    friend class Lock;

  public:
    inline RecursiveMutex() {}

    void lock()   { boost::recursive_mutex::lock(); }
    void unlock() { boost::recursive_mutex::unlock(); }

    // A scoped lock class, used to lock and unlock a Mutex.
    class Lock : private boost::unique_lock<RecursiveMutex>,
                 private boost::noncopyable {

    public:
      inline Lock( RecursiveMutex &mutex ) : boost::unique_lock<RecursiveMutex>( mutex ) {}
      void lock()   { boost::unique_lock<RecursiveMutex>::lock(); }
      void unlock() { boost::unique_lock<RecursiveMutex>::unlock(); }
    };
  };

  // --------------------------------------------------------------
  //                            CONDITION
  // --------------------------------------------------------------

  class Condition : private boost::condition,
                    private boost::noncopyable
  {

  public:
    // construct/copy/destruct
    inline Condition() : boost::condition() {}

    // notification
    void notify_one() {
      boost::condition::notify_one();
    }

    void notify_all() {
      boost::condition::notify_all();
    }

    // waiting
    template<typename LockT> void wait(LockT &lock) {
      boost::condition::wait(lock);
    }

    template<typename LockT, typename Pred>
    void wait(LockT &lock, Pred pred) {
      boost::condition::wait(lock,pred);
    }

    template<typename LockT>
    bool timed_wait(LockT &lock, unsigned long milliseconds) {
      boost::xtime xt;
      boost::xtime_get(&xt, boost::TIME_UTC);
      while (milliseconds >= 1000) {
        xt.sec++;
        milliseconds -= 1000;
      }
      xt.nsec += static_cast<uint32>(1e6) * milliseconds;
      return boost::condition::timed_wait(lock, xt);
    }

    template<typename LockT, typename Pred>
    bool timed_wait(LockT &lock, unsigned long milliseconds, Pred pred) {
      boost::xtime xt;
      boost::xtime_get(&xt, boost::TIME_UTC);
      while (milliseconds >= 1000) {
        xt.sec++;
        milliseconds -= 1000;
      }
      xt.nsec += static_cast<uint32>(1e6) * milliseconds;
      return boost::condition::timed_wait(lock, xt, pred);
    }
  };

  // --------------------------------------------------------------
  //                            THREAD
  // --------------------------------------------------------------

  // A thread class, that runs a "Task", which is an object or
  // function that has the operator() defined.  When the Thread object
  // is destroyed it will join the child thread if it has not already
  // terminated.
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
    /// instance is no longer directly accessibly from the parent
    /// thread.
    template<class TaskT>
    inline Thread( TaskT task ) : m_thread( task ) {}

    /// This variant of the constructor takes a shared point to a task.
    /// The thread made a copy of the shared pointer task, allowing
    /// the parent to still access the task instance that is running in
    /// the thread.
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
      boost::xtime_get(&xt, boost::TIME_UTC);
      while (milliseconds >= 1000) {
        xt.sec++;
        milliseconds -= 1000;
      }
      xt.nsec += static_cast<uint32>(1e6) * milliseconds;
      boost::thread::sleep(xt);
    }
  };


} // namespace vw

#endif // __VW_CORE_THREAD_H__
