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


/// \file vw/Core/Exception.h
///
/// Exception classes and related functions and macros.
///
/// The Vision Workbench is intended in part to be used in flight
/// systems, experimental multiprocessor systems, or other
/// environments where exceptions may not be fully supported.  As a
/// result, the use of exceptions within the Vision Workbench is
/// tightly controlled.  In particular, the exception usage rules were
/// designed to minimize the impact on platforms that do not support
/// exceptions at all.  There is a standard Vision Workbench
/// "exception" class hierarchy which is used to describe errors and
/// can be used even on platforms that do not support the C++
/// exception system.
///
/// The vw::Exception class serves as a base class for all VWB error
/// types.  It is designed to make it easy to throw exceptions with
/// meaningful error messages.  For example, this invocation:
///
///  <TT>vw_throw( vw::Exception() << "Unable to open file \"" << filename << "\"!" );</TT>
///
/// might generate a message like this:
///
///  <TT>terminate called after throwing an instance of 'vw::Exception'</TT>
///  <TT>     what():  Unable to open file "somefile.foo"! </TT>
///
/// A variety of standard derived exception types are provided; in the
/// above example, the exception should probably have been of type
/// vw::IOErr.  Also, two macros, VW_ASSERT(condition,exception) and
/// VW_DEBUG_ASSERT(condition,exception), are provided, with the usual
/// assertion semantics.  The only difference is that the debug
/// assertions will be disabled for increased performance in release
/// builds when VW_DEBUG_LEVEL is defined to zero (which happens by
/// default when NDEBUG is defined).
///
/// Note that in the example the exception was thrown by calling the
/// vw_throw() function rather than by using the C++ throw statement.
/// On platforms that do support C++ exceptions the default behavior
/// for vw_throw() is to throw the exception in the usual way.
/// However, the user can provide their own error-handling mechanism
/// if they choose.  For example, the default behavior when exceptions
/// are disabled is to print the error text to stderr and call
/// abort().
///
/// In general the only allowed usage of exceptions within the Vision
/// Workbench is throwing them using vw_throw().  In particular, try
/// and catch blocks are generally prohibited, so exceptions can only
/// be used to report fatal errors that the library is unable to
/// recover from by itself.  Other uses of exceptions are allowed only
/// under a few special circumstances.  If a part of the Vision
/// Workbench depends on a third-party library that fundamentally
/// relies on exceptions, then that part of the Vision Workbench may
/// use exceptions as well.  However, in that case that entire portion
/// of the Vision Workbench must be disabled when exceptions are not
/// supported.  Similarly, if a part of the Vision Workbench that
/// provides a high-level service cannot reasonably be written without
/// the full use of exceptions, then this portion may also be disabled
/// on platforms without exceptions.  In both of these cases it must
/// be clearly documented that these features are not available on
/// platforms that do not support exceptions.  Finally, it is legal to
/// catch an exception within the library for the sole purpose of
/// re-throwing an exception with a more informative data type or
/// error message.  This purely cosmetic usage must be conditionally
/// compiled like this:
///  #if defined(VW_ENABLE_EXCEPTIONS) && (VW_ENABLE_EXCEPTIONS==1) )
/// Obviously this functionality will be disabled on platforms
/// that do not support exceptions.
///
/// Exceptions are enabled or disabled based on the value of the
/// VW_ENABLE_EXCEPTIONS macro defined in vw/config.h.  This value can be
/// set by passing the command line options --enable-exeptions (the
/// default) or --disable-exceptions to the configure script prior to
/// buliding the Vision Workbench.  This option also sets an automake
/// variable called ENABLE_EXCEPTIONS which may be used by the build
/// system to conditionally compile entire source files.
///
/// In either case the default behavior of vw_throw() may be
/// overridden by calling set_exception_handler(), passing it a
/// pointer to a user-defined object derived from ExceptionHandler.
/// The user specifies the error-handling behavior by overriding the
/// abstract method handle().  When exceptions have not been disabled,
/// the Exception class and its children define a virtual method
/// default_throw() which the handler may call to have the exception
/// throw itself in a type-aware manner.
///
#ifndef __VW_CORE_EXCEPTION_H__
#define __VW_CORE_EXCEPTION_H__

#include <vw/Core/Features.h>
#include <vw/vw_config.h>

#include <string>
#include <sstream>
#include <ostream>

#if defined(VW_ENABLE_EXCEPTIONS)
#include <exception>
#define VW_IF_EXCEPTIONS(x) x
#else
#define VW_IF_EXCEPTIONS(x)
#endif

namespace vw {

  /// The core exception class.
  struct Exception VW_IF_EXCEPTIONS( : public std::exception )
  {

    /// The default constructor generates exceptions with empty error
    /// message text.  This is the cleanest approach if you intend to
    /// use streaming (via operator <<) to generate your message.
    Exception() VW_NOTHROW {}

    virtual ~Exception() VW_NOTHROW {}

    /// Copy Constructor
    Exception( Exception const& e ) VW_NOTHROW
      VW_IF_EXCEPTIONS( : std::exception(e) ) {
      m_desc << e.m_desc.str();
    }

    /// Assignment operator copies the error string.
    Exception& operator=( Exception const& e ) VW_NOTHROW {
      m_desc.str( e.m_desc.str() );
      return *this;
    }

    /// Returns a the error message text for display to the user.  The
    /// returned pointer must be used immediately; other operations on
    /// the exception may invalidate it.  If you need the data for
    /// later, you must save it to a local buffer of your own.
    virtual const char* what() const VW_NOTHROW {
      m_what_buf = m_desc.str();
      return m_what_buf.c_str();
    }

    /// Returns the error message text as a std::string.
    std::string desc() const { return m_desc.str(); }

    /// Returns a string version of this exception's type.
    virtual std::string name() const { return "Exception"; }

    void set( std::string const& s ) { m_desc.str(s); }
    void reset() { m_desc.str(""); }

    VW_IF_EXCEPTIONS( virtual void default_throw() const { throw *this; } )

  protected:
      virtual std::ostringstream& stream() {return m_desc;}

  private:
    // The error message text.
    std::ostringstream m_desc;

    // A buffer for storing the full exception description returned by
    // the what() method, which must generate its return value from
    // the current value of m_desc.  The what() method provides no
    // mechanism for freeing the returned string, and so we handle
    // allocation of that memory here, internally to the exception.
    mutable std::string m_what_buf;
  };

  // Use this macro to construct new exception types that do not add
  // additional functionality.  If you can think of a clean way to do
  // this using templates instead of the preprocessor, please do.  For
  // now, we're stuck with this.
  //
  // Some functions need to return the *this pointer with the correct
  // subclass type, and these are defined in the macro below rather
  // than the base exception class above.   These are:
  //
  // Exception::operator=():
  // The assignment operator must return an instance of the subclass.
  //
  // Exception::operator<<():
  // The streaming operator (<<) makes it possible to quickly
  // generate error message text.  This is currently implemented
  // by simply forwarding invocations of this method to an
  // internal ostringstream.

  /// Macro for quickly creating a hierarchy of exceptions, all of
  /// which share the same functionality.
  #define VW_DEFINE_EXCEPTION(exception_type,base)                             \
    struct exception_type : public base {                                      \
      VW_EXCEPTION_API(exception_type)                                         \
    }

  // This sentinel catches users who forget to use VW_EXCEPTION_API
  #define _VW_EXCEPTION_SENTINEL(e) _YouForgot_VW_EXCEPTION_API_OnException_ ## e

  // Macro for creating a hierarchy of exceptions that may have additional
  // functions or data. When using this macro, you must include the
  // VW_EXCEPTION_API macro inside the braces
  #define VW_DEFINE_EXCEPTION_EXT(exception_type,base)                         \
    struct _VW_EXCEPTION_SENTINEL(exception_type) : public base {              \
      virtual std::string name() const = 0;                                    \
    };                                                                         \
    struct exception_type : public _VW_EXCEPTION_SENTINEL(exception_type)

  #define VW_EXCEPTION_API(exception_type)                                     \
    virtual std::string name() const { return #exception_type; }               \
    VW_IF_EXCEPTIONS( virtual void default_throw() const { throw *this; } )    \
    template <class T>                                                         \
    exception_type& operator<<( T const& t ) { stream() << t; return *this; }

  /// Invalid function argument exception
  VW_DEFINE_EXCEPTION(ArgumentErr, Exception);

  /// Incorrect program logic exception
  VW_DEFINE_EXCEPTION(LogicErr, Exception);

  /// Invalid program input exception
  VW_DEFINE_EXCEPTION(InputErr, Exception);

  /// IO failure exception
  VW_DEFINE_EXCEPTION(IOErr, Exception);

  /// Arithmetic failure exception
  VW_DEFINE_EXCEPTION(MathErr, Exception);

  /// Unexpected null pointer exception
  VW_DEFINE_EXCEPTION(NullPtrErr, Exception);

  /// Invalid type exception
  VW_DEFINE_EXCEPTION(TypeErr, Exception);

  /// Not found exception
  VW_DEFINE_EXCEPTION(NotFoundErr, Exception);

  /// Unimplemented functionality exception
  VW_DEFINE_EXCEPTION(NoImplErr, Exception);

  /// Operation aborted partway through (e.g. with ProgressCallback
  /// returning Abort)
  VW_DEFINE_EXCEPTION(Aborted, Exception);

  /// The abstract exception handler base class, which users
  /// can subclass to install an alternative exception handler.
  class ExceptionHandler {
  public:
    virtual void handle( Exception const& e ) const VW_NORETURN = 0;
    virtual ~ExceptionHandler() VW_NOTHROW {}
  };

  /// Sets the application-wide exception handler.  Pass zero
  /// as an argument to reinstate the default handler.  The
  /// default behavior is to throw the exception unless the
  /// VW_ENABLE_EXCEPTIONS macro in vw/config.h was defined to 0
  /// at build time, in which case the default behavior is to
  /// print the error message at the ErrorMessage level and
  /// to call abort().
  void set_exception_handler( ExceptionHandler const* eh );

  /// Throws an exception via the Vision Workbench error
  /// handling mechanism, which may not actually involvle
  /// throwing an exception in the usual C++ sense.
  void vw_throw( Exception const& e ) VW_NORETURN;

} // namespace vw

/// The VW_ASSERT macro throws the given exception if the given
/// condition is not met.  The VW_DEBUG_ASSERT macro does the same
/// thing, but is disabled if VW_DEBUG_LEVEL is zero.  The default
/// value for VW_DEBUG_LEVEL is defined in Debugging.h.
#define VW_ASSERT(cond,excep) do { if(!(cond)) vw::vw_throw( excep ); } while(0)
#define VW_LINE_ASSERT(cond) do { if(!(cond)) vw::vw_throw( vw::LogicErr() << "Assertion failed (" << __FILE__ << ":" << __LINE__ << "): " << #cond); } while(0)
#if VW_DEBUG_LEVEL == 0
#define VW_DEBUG_ASSERT(cond,excep) do {} while(0)
#else
// Duplicate the definition to avoid extra macro expansion in recusion
#define VW_DEBUG_ASSERT(cond,excep) do { if(!(cond)) vw::vw_throw( excep ); } while(0)
#endif

#endif // __VW_CORE_EXCEPTION_H__
