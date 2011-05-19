// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
