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

/// \file Core/Threads.h
/// 
/// A very simple, abstracted threading system.  Currently this is
/// simply built on top of Boost threads, but its abstracted so that
/// it's easy to port the Vision Workbench to environments where Boost
/// threads do not work or otherwise cannot be used, and because the
/// Boost threads design itself is expected to be overhauled at some
/// point.
///
/// To create a thread, first subclass the Runnable abstract class.
/// It has only two methods: the run() method, which will be run in
/// the child thread, and the kill() method, which should signal the
/// running thread to terminate and clean up.  Note that the kill()
/// function itself will run in a separate thread from the run()
/// function, so it should signify this using an atomic or
/// synchronized mechanism.  Also note that the thread may have
/// already terminated when kill() is called.
///
/// The rest of the interface is sraightforward.  These are the 
/// the things you will need to implement if you want to use a 
/// different thread library.  The interface consists of:
///
/// * A Thread class, whose constructor takes a shared pointer to a
///   Runnable.  The constructor should invoke the Runnable's run()
///   function in a new thread, and the destructor should invoke 
///   the kill() function and join the child thread.  The static 
///   yeild() function should yield the active thread, and the static 
///   sleep_ms() functin should sleep the active thread for the given 
///   number of milliseconds.  The Thread class may be non-copyable.
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

#ifndef __VW_CORE_THREADS_H__
#define __VW_CORE_THREADS_H__

#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>

namespace vw {

  // A thread class, that runs a Runnable.  The Runnable is run() in
  // the new thread.  When the Thread object is destroyed it will call
  // the runnable's kill() function and join the child thread if it
  // has not already terminated.
  class Thread {
  public:

    // An abstract runnable base class.
    class Runnable {
    public:
      virtual ~Runnable() {}
      virtual void run() = 0;  // Runs in the child thread
      virtual void kill() = 0; // Runs in the calling thread
    };

  private:
    boost::shared_ptr<Runnable> m_runnable;
    boost::thread m_thread;

    // Ensure non-copyable semantics
    Thread( Thread const& );
    Thread& operator=( Thread const& );

    // A wrapper so we can can construct a boost::thread
    class RunFunctor {
      Runnable& m_runnable;
    public:
      inline RunFunctor( Runnable& runnable ) : m_runnable( runnable ) {}
      inline void operator()() const { m_runnable.run(); }
    };

  public:
    inline Thread( boost::shared_ptr<Runnable> runnable ) : m_runnable( runnable ), m_thread( RunFunctor( *runnable ) ) {}
    inline ~Thread() { m_runnable->kill(); m_thread.join(); }
    static inline void yield() { boost::thread::yield(); }
    static inline void sleep_ms( unsigned long milliseconds ) {
      static long const nanoseconds_per_second = 1000L*1000L*1000L;
      boost::xtime xt;
      boost::xtime_get(&xt, boost::TIME_UTC);
      xt.nsec+=1000*1000*milliseconds;
      while (xt.nsec > nanoseconds_per_second) {
        xt.nsec -= nanoseconds_per_second;
        xt.sec++;
      }
      boost::thread::sleep(xt);
    }
    static int default_num_threads();
    static void set_default_num_threads(int);
  };


  // A simple mutual exclusion class.
  class Mutex : private boost::mutex {

    // Ensure non-copyable semantics
    Mutex( Mutex const& );
    Mutex& operator=( Mutex const& );

  public:
    inline Mutex() {}

    // A scoped lock class, used to lock and unlock a Mutex.
    class Lock : private boost::mutex::scoped_lock {

      // Ensure non-copyable semantics
      Lock( Lock const& );
      Lock& operator=( Lock const& );

    public:
      inline Lock( Mutex &mutex ) : boost::mutex::scoped_lock( mutex ) {}
    };
    friend class Lock;

  };


  // A special POD class to enable safe library initialization.  You 
  // should only define these objects at global or namespace scope, 
  // and statically initialize them to VW_RUNONCE_INIT.
  struct RunOnce {
    boost::once_flag m_flag;

    inline void run( void (*func)() ) {
      boost::call_once( func, m_flag );
    }
  };

#define VW_RUNONCE_INIT { BOOST_ONCE_INIT }

} // namespace vw

#endif // __VW_CORE_THREADS_H__
