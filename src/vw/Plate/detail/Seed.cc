// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/detail/Seed.h>
#include <vw/Core/Thread.h>
#include <vw/Core/Stopwatch.h>
#include <sys/types.h>

namespace {
  vw::RunOnce seed_once = VW_RUNONCE_INIT;

  void hidden_seed_random() {
    srandom(vw::Stopwatch::microtime()*::getpid());
  }
}


namespace vw {
namespace platefile {

  void plate_seed_random() {
    seed_once.run( hidden_seed_random );
  }

}}
