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


#include <vw/Plate/detail/Seed.h>
#include <vw/Core/Thread.h>
#include <vw/Core/Log.h>
#include <vw/Core/Exception.h>
#include <vw/Core/Stopwatch.h>
#include <unistd.h>

namespace {
  using namespace vw;
  RunOnce seed_once = VW_RUNONCE_INIT;

  void hidden_seed_random() {
    uint64 seed = Stopwatch::microtime();
    unsigned int seed32 = (unsigned int)(seed & 0x00000000ffffffff);
    unsigned int pid = (unsigned int)::getpid();
    seed32 ^= pid;
    vw_out(DebugMessage) << "Seeding RNG with " << seed32 << std::endl;
    ::srandom(seed32);
  }
}


namespace vw {
namespace platefile {

  void plate_seed_random() {
    seed_once.run( hidden_seed_random );
  }

}}
