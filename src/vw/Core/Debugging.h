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
#include <vw/Core/Features.h>

#ifndef WIN32
#include <sys/time.h>
#endif

#define _VW_STRINGIFY(x) #x
#define VW_STRINGIFY(x) _VW_STRINGIFY(x)

#if defined(_WIN32)
# define VW_CURRENT_FUNCTION __FUNCSIG__
#elif defined(__GNUC__)
# define VW_CURRENT_FUNCTION __PRETTY_FUNCTION__
#else
# define VW_CURRENT_FUNCTION (__FILE__ ":" VW_STRINGIFY(__LINE__))
#endif


#include <ostream>
#include <string>

namespace vw {


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
