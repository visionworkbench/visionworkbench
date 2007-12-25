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

#include <vw/Core/Thread.h>
#include <vw/config.h>

namespace {
  static int num_threads = VW_NUM_THREADS;
}

namespace vw { 
namespace thread {

  // These static variables are used to generate unique identifiers
  // for threads as identifiers are requested using the static
  // Thread::id() method.
  static int vw_thread_next_available_id = 0;
  static boost::thread_specific_ptr<int> vw_thread_id_ptr;
  static Mutex vw_thread_id_mutex; 

}} // namespace vw::thread


int vw::Thread::id() { 
  
  // If the thread ID has not been accessed before, we initialize
  // it with a unique ID.
  if (thread::vw_thread_id_ptr.get() == 0) {
    Mutex::Lock lock(thread::vw_thread_id_mutex);
    thread::vw_thread_id_ptr.reset(new int(thread::vw_thread_next_available_id++));
  }
  
  // Then we return the result.
  int* result = thread::vw_thread_id_ptr.get();   
  return *result;
}


int vw::Thread::default_num_threads() {
  return num_threads;
}

void vw::Thread::set_default_num_threads(int threads) {
  num_threads = threads;
}
