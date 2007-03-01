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

namespace {
  // A null output stream buffer that silently ignores any data.
  class NullStreamBuf : public ::std::basic_streambuf<char> {
    virtual NullStreamBuf::int_type overflow( NullStreamBuf::int_type c ) {
      return NullStreamBuf::traits_type::not_eof( c );
    }
  };

  static NullStreamBuf g_the_nullbuf;
  static std::ostream g_the_nullstream( &g_the_nullbuf );
  static std::ostream g_the_ostream( std::clog.rdbuf() );
  static int g_the_level = vw::InfoMessage;
}

std::ostream& vw::vw_out( int level ) {
  if( level > g_the_level ) return g_the_nullstream;
  else return g_the_ostream;
}

void vw::set_debug_level( int level ) {
  g_the_level = level;
}

void vw::set_output_stream( std::ostream& stream ) {
  g_the_ostream.rdbuf( stream.rdbuf() );
}



vw::Timer::Timer( std::string const& desc, MessageLevel level )
  : m_desc(desc), m_level(level) {
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
  vw_out(m_level) << m_desc << ": " << duration << std::endl;
}
