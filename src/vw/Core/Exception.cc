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

/// \file Core/Exception.cc
/// 
/// Error / exception handling facilities.
///
/// See Core/Exception.h for documentation.
///
#include <vw/Core/Exception.h>

#include <cstdlib>

namespace {

  /// The default exception handler type and object, which throws 
  /// the exceptions its given unless VW_NO_EXCEPTIONS is 1, in 
  /// which case it prints the message and calls abort().
  static class DefaultExceptionHandler : public vw::ExceptionHandler {
  public:
    virtual void handle( vw::Exception const& e ) const {
#if defined(VW_NO_EXCEPTIONS) && (VW_NO_EXCEPTIONS==1)
      vw::vw_out(vw::ErrorMessage) << "Fatal error: " << e.what() << std::endl;
      std::abort();
#else
      e.default_throw();
#endif
    }
    virtual ~DefaultExceptionHandler() {}
  } _vw_default_exception_handler;

  /// The application-wide exception handler pointer.
  static vw::ExceptionHandler const* _vw_exception_handler = &_vw_default_exception_handler;

};

void vw::set_exception_handler( vw::ExceptionHandler const* eh ) {
  if( eh ) _vw_exception_handler = eh;
  else _vw_exception_handler = &_vw_default_exception_handler;
}

void vw::vw_throw( vw::Exception const& e ) {
  _vw_exception_handler->handle( e );
}
