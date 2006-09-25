// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
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

#include <sys/time.h>
#include <ostream>

namespace vw {

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

  std::ostream& print( MessageLevel level );
  void set_debug_level( MessageLevel level );
  void set_output_stream( std::ostream& stream );


  // *******************************************************************
  // Timing types and functions
  // *******************************************************************

  class Timer {
    std::string m_desc;
    MessageLevel m_level;
    timeval m_begin;
  public:
    Timer( std::string const& desc, MessageLevel level=InfoMessage )
      : m_desc(desc), m_level(level) {
      gettimeofday( &m_begin, 0 );
    }

    ~Timer() {
      timeval end;
      gettimeofday( &end, 0 );
      double duration = end.tv_sec - m_begin.tv_sec;
      duration += (end.tv_usec - m_begin.tv_usec)/1.0e6;
      print(m_level) << m_desc << ": " << duration << std::endl;
    }
  };


} // namespace vw

#endif  // __VW_CORE_DEBUGGING_H__
