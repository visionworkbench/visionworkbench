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

/// \file Core/Threads.cc

#if defined( VW_HAVE_PKG_PTHREADS )

namespace vw {
  namespace detail {

    class ThreadImpl {
      pthread_t m_thread;

      static void *run_func( void *runnable_ptr ) {
        ((vw::Runnable*)runnable_ptr)->run();
      }

    public:
      ThreadImpl( Runnable &runnable ) : m_runnable( runnable ) {
        int result = pthread_create( &m_thread, 0, &run_func, &m_runnable );
        if( result != 0 ) vw_throw( SystemErr() << "Error creating a pthread (" << result << ")." );
      }

      ~ThreadImpl() {
        m_runnable->kill();
        int result = pthread_join( &m_thread, 0 );
        if( result != 0 ) vw_throw( SystemErr() << "Error joining a pthread (" << result << ")." );
      }
    };

    class MutexImpl {
      pthread_mutex_t m_mutex;

    public:
      MutexImpl() {
        int result = pthread_mutex_init( &m_mutex, 0 );
        if( result != 0 ) vw_throw( SystemErr() << "Unable to initialize a pthread mutex (" << result << ")." );
      }

      ~MutexImpl() {
        int result = pthread_mutex_destroy( &m_mutex );
        if( result != 0 ) vw_throw( SystemErr() << "Error desroying a pthread mutex (" << result << ")." );
      }

      void lock() {
        int result = pthread_mutex_lock( &m_mutex );
        if( result != 0 ) vw_throw( SystemErr() << "Error locking a pthread mutex (" << result << ")." );
      }
    
      void unlock() {
        int result = pthread_mutex_unlock( &m_mutex );
        if( result != 0 ) vw_throw( SystemErr() << "Error locking a pthread mutex (" << result << ")." );
      }
      
    };

  } // namespace detail
} // namespace vw

#elif defined( WIN32 ) || defined( _WIN32 ) || defined( __WIN32__ )

#error No threading support yet under Windows!

#else

#error No threading support on this platform!

#endif

vw::Thread::Thread( Runnable &runnable ) : m_thread( new vw::detail::ThreadImpl( runnable ) ) {}
vw::Mutex::Mutex() : m_mutex( new vw::detail::MutexImpl() ) {}
vw::Mutex::lock() { m_mutex->lock(); }
vw::Mutex::unlock() { m_mutex->unlock(); }
