// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_CORE_TERMINATION_HANDLER_H__
#define __VW_CORE_TERMINATION_HANDLER_H__

#include <boost/function.hpp>

namespace vw {

  // Handlers for terminate().
  typedef boost::function<void()> TerminateHandler;

  // Some built-in handlers to choose from.
  extern const TerminateHandler DefaultTerminate;
  extern const TerminateHandler PrettyTerminate;
  extern const TerminateHandler BacktraceTerminate;

  // Sets the application-wide terminate() handler. Pass zero to reinstate
  // the default. If you want the program to hang so you can attach gdb, your
  // TerminateHandler should not call abort!  The return value will tell you
  // whether the new handler was actually installed. This is useful if
  // PrettyTerminate or BacktraceTerminate are not supported in your build.
  bool set_terminate_handler( TerminateHandler const &th );

} // namespace vw


#endif
