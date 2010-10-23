// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Core/Exception.cc
///
/// Error / exception handling facilities.
///
/// See Core/Exception.h for documentation.
///
#include <vw/Core/Exception.h>
#include <vw/Core/Features.h>

#include <cstdlib>

namespace {

  /// The default exception handler type and object, which throws
  /// the exceptions its given unless VW_ENABLE_EXCEPTIONS is 0, in
  /// which case it prints the message and calls abort().
  static class DefaultExceptionHandler : public vw::ExceptionHandler {
  public:
    virtual void handle( vw::Exception const& e ) const VW_NORETURN {
#if defined(VW_ENABLE_EXCEPTIONS) && (VW_ENABLE_EXCEPTIONS==1)
      e.default_throw();
#else
      vw::vw_out(vw::ErrorMessage) << "Fatal error: " << e.what() << std::endl;
#endif
      std::abort();
    }
    virtual ~DefaultExceptionHandler() VW_NOTHROW {}
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
  // We cannot return.
  std::abort();
}
