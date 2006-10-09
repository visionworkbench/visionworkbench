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

/// \file Core/Debugging.h
/// 
/// Types and functions to assist in debugging code.
///
#ifndef __VW_CORE_DEBUGGING_H__
#define __VW_CORE_DEBUGGING_H__

#ifdef WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#include <ostream>
#include <string>

namespace vw {

/// The master compile-time debugging level flag.  The default value
/// for __VW_DEBUG_LEVEL__ is guessed based on whether or not NDEBUG
/// is defined if the user has not specified it explicitly.
#ifndef __VW_DEBUG_LEVEL__
#ifdef NDEBUG
#define __VW_DEBUG_LEVEL__ 0
#else
#define __VW_DEBUG_LEVEL__ 1
#endif
#endif

/// A quick macro for selectively disabling code in non-debug builds.
#if __VW_DEBUG_LEVEL__ == 0
#define VW_DEBUG(x)
#else
#define VW_DEBUG(x) x
#endif

  // *******************************************************************
  // Debugging output types and functions
  // *******************************************************************

  enum MessageLevel {
    ErrorMessage = 0,
    WarningMessage = 10,
    InfoMessage = 20,
    DebugMessage = 30,
    VerboseDebugMessage = 40
  };

  std::ostream& vw_out( MessageLevel level );
  void set_debug_level( MessageLevel level );
  void set_output_stream( std::ostream& stream );


  // *******************************************************************
  // Timing types and functions
  // *******************************************************************

  class Timer {
    std::string m_desc;
    MessageLevel m_level;
#ifdef WIN32
    LARGE_INTEGER m_begin;
#else
    timeval m_begin;
#endif
  public:
    Timer( std::string const& desc, MessageLevel level=InfoMessage )
      : m_desc(desc), m_level(level) {
#ifdef WIN32
      QueryPerformanceCounter( &m_begin );
#else
      gettimeofday( &m_begin, 0 );
#endif
    }

    ~Timer() {
#ifdef WIN32
      LARGE_INTEGER end, freq;
      QueryPerformanceCounter( &end );
      QueryPerformanceFrequency( &freq );
      double duration = (end.QuadPart - m_begin.QuadPart)/(double)freq.QuadPart;
#else
      timeval end;
      gettimeofday( &end, 0 );
      double duration = end.tv_sec - m_begin.tv_sec;
      duration += (end.tv_usec - m_begin.tv_usec)/1.0e6;
#endif
      vw_out(m_level) << m_desc << ": " << duration << std::endl;
    }
  };


} // namespace vw

#endif  // __VW_CORE_DEBUGGING_H__
