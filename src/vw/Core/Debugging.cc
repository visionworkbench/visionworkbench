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

/// \file Core/Debugging.cc
/// 
/// Types and functions to assist in debugging code.
///
#include <vw/Core/Debugging.h>

#include <iostream>

#ifdef WIN32
#include <Windows.h>
#endif

vw::Timer::Timer( std::string const& desc, MessageLevel level, std::string const& log_namespace )
  : m_desc(desc), m_level(level), m_log_namespace(log_namespace) {
#ifdef WIN32
  LARGE_INTEGER begin;
  QueryPerformanceCounter( &begin );
  m_begin= begin.QuadPart;
#else
  gettimeofday( &m_begin, 0 );
#endif
}

vw::Timer::~Timer() {
#ifdef WIN32
  LARGE_INTEGER end, freq;
  QueryPerformanceCounter( &end );
  QueryPerformanceFrequency( &freq );
  double duration = (end.QuadPart - m_begin)/(double)freq.QuadPart;
#else
  timeval end;
  gettimeofday( &end, 0 );
  double duration = end.tv_sec - m_begin.tv_sec;
  duration += (end.tv_usec - m_begin.tv_usec)/1.0e6;
#endif
  vw_out(m_level, m_log_namespace) << m_desc << ": " << duration << std::endl;
}
