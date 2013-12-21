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


/// \file Core/Debugging.cc
///
/// Types and functions to assist in debugging code.
///
#include <vw/Core/Debugging.h>
#include <vw/Core/Log.h>

#ifdef WIN32
#include <Windows.h>
#else
#include <sys/time.h>
#endif

#include <ostream>

#include <boost/numeric/conversion/cast.hpp>

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
  double duration = boost::numeric_cast<double>(end.tv_sec - m_begin.tv_sec);
  duration += boost::numeric_cast<double>(end.tv_usec - m_begin.tv_usec)/1.0e6;
#endif
  VW_OUT(m_level, m_log_namespace) << m_desc << ": " << duration << std::endl;
}
