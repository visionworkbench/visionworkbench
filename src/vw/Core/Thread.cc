// __BEGIN_LICENSE__
//  Copyright (c) 2006-2026, United States Government as represented by the
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
#include <vw/Core/Thread.h>
#include <vw/Core/FundamentalTypes.h>

namespace vw {
namespace thread {

  // These static variables are used to generate unique identifiers
  // for threads as identifiers are requested using the static
  // Thread::id() method.
  static vw::uint64 vw_thread_next_available_id = 0;

  // Both the mutex and the thread-local storage need to be available to all
  // destructors (in case they call Thread::id()). If the object being
  // destructed is static and defined in a different file, the destruction
  // order is undefined. That causes a destruction race. To prevent this race,
  // we use the construct-on-first-use idiom. See
  // http://www.parashift.com/c++-faq-lite/ctors.html#faq-10.15 for more info.
  static Mutex& vw_thread_id_mutex() {
    static Mutex* m = new Mutex();
    return *m;
  }

}} // namespace vw::thread

vw::uint64 vw::Thread::id() {

  // Assign a unique ID the first time this thread asks for one, kept in
  // thread-local storage. The global counter is bumped under a lock.
  static thread_local bool assigned = false;
  static thread_local vw::uint64 the_id = 0;
  if (!assigned) {
    Mutex::Lock lock(thread::vw_thread_id_mutex());
    the_id = thread::vw_thread_next_available_id++;
    assigned = true;
  }
  return the_id;
}
