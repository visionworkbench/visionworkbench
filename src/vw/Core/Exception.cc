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


/// \file Core/Exception.cc
///
/// Error / exception handling facilities.
///
/// See Core/Exception.h for documentation.
///
#include <vw/vw_config.h>
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
#if defined(VW_ENABLE_EXCEPTIONS)
      e.default_throw();
#else
      vw::VW_OUT(vw::ErrorMessage) << "Fatal error: " << e.what() << std::endl;
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
