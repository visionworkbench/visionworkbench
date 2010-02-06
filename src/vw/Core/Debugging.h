// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Core/Debugging.h
/// 
/// Types and functions to assist in debugging code.
///
#ifndef __VW_CORE_DEBUGGING_H__
#define __VW_CORE_DEBUGGING_H__

#include <vw/Core/Log.h>

#ifndef WIN32
#include <sys/time.h>
#endif

#include <ostream>
#include <string>

namespace vw {

/// The master compile-time debugging level flag.  The default value
/// for VW_DEBUG_LEVEL is guessed based on whether or not NDEBUG
/// is defined if the user has not specified it explicitly.
#ifndef VW_DEBUG_LEVEL
#ifdef NDEBUG
#define VW_DEBUG_LEVEL 0
#else
#define VW_DEBUG_LEVEL 1
#endif
#endif

/// A quick macro for selectively disabling code in non-debug builds.
#if VW_DEBUG_LEVEL == 0
#define VW_DEBUG(x)
#else
#define VW_DEBUG(x) x
#endif

  // *******************************************************************
  // Timing types and functions
  // *******************************************************************

  class Timer {
    std::string m_desc;
    MessageLevel m_level;
    std::string m_log_namespace;
#ifdef WIN32
    __int64 m_begin;
#else
    timeval m_begin;
#endif
  public:
    Timer( std::string const& desc, MessageLevel level=InfoMessage, std::string const& log_namespace = "console" );

    ~Timer();
  };


} // namespace vw

#endif  // __VW_CORE_DEBUGGING_H__
