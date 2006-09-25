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

/// \file Exception.h
/// 
/// Base exception classes and related macros.
/// 
/// The vw::Exception class serves as a base class for all VWB error
/// types.  It is designed to make it easy to throw exceptions with 
/// meaningful error messages.  For example, the this invocation:
///
///   <TT>throw vw::Exception() << "Unable to open file \"" << filename << "\"!";</TT>
///
/// might generate a message like this:
///
///  <TT>terminate called after throwing an instance of 'vw::Exception'</TT>
///
///  <TT>     what():  Unable to open file "somefile.foo"! </TT>
///
/// A variety of standard derived exception types are provided; in the
/// above example, the exception should probably have been of type
/// vw::IOErr.  Also, two macros, VW_ASSERT(condition,exception) and 
/// VW_DEBUG_ASSERT(condition,exception), are provided, with the usual 
/// assertion semantics.  The only difference is that the debug assertions 
/// will be disabled for increased performance in release builds when 
/// __VW_DEBUG_LEVEL__ is defined to zero (which happens by default when 
/// NDEBUG is defined).
///
/// Note that the base Exception class is not particularly lightweight.
/// If you are throwing an exception that you intend to catch as part
/// of the normal flow of your program, you should create your own 
/// (typically empty) class for this purpose.  If you do that, however, 
/// you must be absolutely certain that you catch it everywhere it may 
/// be thrown.  The user should only ever encounter exceptions derived 
/// from vw::Exception.
///
#ifndef __VW_CORE_EXCEPTION_H__
#define __VW_CORE_EXCEPTION_H__

#include <exception>
#include <string>
#include <sstream>
#include <ostream>

namespace vw {

  /// The core exception class.  
  struct Exception : public std::exception {

    /// The default constructor generates exceptions with empty error
    /// message text.  This is the cleanest approach if you intend to
    /// use streaming (via operator <<) to generate your message.
    Exception() throw() {}

    /// Generates exceptions with the given error message text.
    Exception( std::string const& s ) throw() { m_desc << s; }

    virtual ~Exception() throw() {}

    /// Copy Constructor
    Exception( Exception const& e ) throw() {
      m_desc << e.m_desc.str();
    }

    /// Assignment operator copies the error string.
    Exception& operator=( Exception const& e ) throw() {
      m_desc.str( e.m_desc.str() );
      return *this;
    }

    /// Returns a the error message text for display to the user.  The
    /// returned pointer must be used immediately; other operations on
    /// the exception may invalidate it.  If you need the data for
    /// later, you must save it to a local buffer of your own.
    virtual const char* what() const throw() {
      m_what_buf = m_desc.str();
      return m_what_buf.c_str();
    }

    /// Returns the error message text as a std::string.
    std::string desc() const { return m_desc.str(); }

  protected:
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
  //
  // Exception::set():
  // Sets the error message text to the provided string, returning a
  // reference to the exception for use with the streaming operator
  // (<<) if desired.
  //
  // Exception::reset():
  // Resets (i.e. clears) the error message text, returning a
  // reference to the exception for use with the streaming operator
  // (<<) if desired

  /// Macro for quickly creating a hierarchy of exceptions, all of
  /// which share the same functionality.
  #define VW_DEFINE_EXCEPTION(name,base)                              \
  struct name : public base {                                         \
    name() throw() : base() {}                                        \
    name(std::string const& s) throw() : base(s) {}                   \
    name( name const& e ) throw() : base( e ) {}                      \
    virtual ~name() throw() {}                                        \
                                                                      \
    inline name& operator=( name const& e ) throw() {                 \
      base::operator=( e );                                           \
      return *this;                                                   \
    }                                                                 \
                                                                      \
    template <class T>                                                \
    name& operator<<( T const& t ) { m_desc << t; return *this; }     \
                                                                      \
    name& set( std::string const& s ) { m_desc.str(s);  return *this; } \
                                                                      \
    name& reset() { m_desc.str("");  return *this; }                  \
  }

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

} // namespace vw

/// The VW_ASSERT macro throws the given exception if the given 
/// condition is not met.  The VW_DEBUG_ASSERT macro does the 
/// same thing, but is disabled if __VW_DEBUG_LEVEL__ is zero.
/// The default value for __VW_DEBUG_LEVEL__ is guessed based 
/// on whether or not NDEBUG is defined.
#ifndef __VW_DEBUG_LEVEL__
#ifdef NDEBUG
#define __VW_DEBUG_LEVEL__ 0
#else
#define __VW_DEBUG_LEVEL__ 1
#endif
#endif

#define VW_ASSERT(cond,excep) do { if(!(cond)) throw (excep); } while(0)
#if __VW_DEBUG_LEVEL__ == 0
#define VW_DEBUG_ASSERT(cond,excep) do {} while(0)
#else
// Duplicate the definition to avoid extra macro expansion in recusion
#define VW_DEBUG_ASSERT(cond,excep) do { if(!(cond)) throw (excep); } while(0)
#endif

#endif // __VW_CORE_EXCEPTION_H__
