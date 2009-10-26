// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <cstdio>
#include <stdexcept>
#include <cstdlib>

#include <typeinfo>

#include <vw/config.h>

#ifdef VW_HAVE_CXXABI_H
#include <cxxabi.h>
#endif

#ifdef VW_HAVE_EXECINFO_H
// for backtrace
#include <execinfo.h>
#endif

#ifdef VW_HAVE_DLFCN_H
#include <dlfcn.h>
#endif

#ifdef VW_HAVE_UNISTD_H
#include <unistd.h>
#endif

#include "TerminationHandler.h"


#if defined(VW_HAVE___CXA_DEMANGLE) && defined(VW_HAVE_CXXABI_H)
#  if defined(VW_HAVE___CXA_CURRENT_EXCEPTION_TYPE)
#    define USE_INTROSPECTION_EXCEPTION
#  endif
#  if defined(VW_HAVE_BACKTRACE) && defined(VW_HAVE_DLADDR) && defined(VW_HAVE_DLFCN_H) && defined(VW_HAVE_EXECINFO_H)
#    define USE_INTROSPECTION_BACKTRACE
#  endif
#endif

namespace {

// The original default terminate() handler. Should be std::abort, but isn't
// always
void (*real_default)()       = 0;

// Currently installed handler
vw::TerminateHandler current = 0;

void vw_terminate()
{
  static bool terminating = false;
  if (terminating) {
      fputs("terminate called recursively!\n", stderr);
      abort();
  }
  terminating = true;

  if (current)
    current();

  fprintf(stderr, "terminate() called! You can attach a debugger.");
#ifdef VW_HAVE_GETPID
  fprintf(stderr, "If you'd like to connect gdb to this process, use:\n");
  fprintf(stderr, "gdb ignored %d\n", getpid());
#endif
  fprintf(stderr, "or press 'enter' to abort...\n");
  // It's probably not safe to getchar here, but, given where we are,
  // we can't really make it worse...
  getchar();
  abort();
}

// This is a bit of a hack, but it lets us distinguish between
// DefaultTerminate and an attempt to unset the terminate handler
void noop() {}

// A replacement for the standard terminate_handler which prints
// more information about the terminating exception (if any) on
// stderr.
void print_exception_info()
{
#ifdef USE_INTROSPECTION_EXCEPTION
  // Make sure there was an exception; terminate is also called for an
  // attempt to rethrow when there is no suitable exception.
  std::type_info *t = abi::__cxa_current_exception_type();
  if (t)
  {
    // Note that "name" is the mangled name.
    char const *name = t->name();
    {
      int status = -1;
      char *dem  = abi::__cxa_demangle(name, 0, 0, &status);

      fputs("terminate called after throwing an instance of '", stderr);
      if (status == 0)
        fputs(dem, stderr);
      else
        fputs(name, stderr);
      fputs("'\n", stderr);

      if (status == 0)
        free(dem);
    }

#if defined(VW_ENABLE_EXCEPTIONS) && (VW_ENABLE_EXCEPTIONS == 1)
    // If the exception is derived from std::exception, we can
    // give more information.
      try { __throw_exception_again; }
      catch (std::exception const &exc)
        {
          char const *w = exc.what();
          fputs("  what():  ", stderr);
          fputs(w, stderr);
          fputs("\n", stderr);
        }
      catch (...) { }
#endif // VW_ENABLE_EXCEPTIONS
    }
  else
    fputs("terminate called without an active exception\n", stderr);
#else
    fputs("terminate() called, but no exception support\n", stderr);
#endif  //USE_INTROSPECTION_EXCEPTION
}

// Adapted by Randy Sargent 5/29/2006
// Much comes from glibc verbose termination handler
void do_backtrace()
{
#ifdef USE_INTROSPECTION_BACKTRACE
  print_exception_info();

  void *backtrace_buffer[200];
  int backtrace_size = backtrace(backtrace_buffer, 200);

  if (backtrace_size)
  {
    fprintf(stderr, "Stack backtrace:\n");
    for (int i = 0; i < backtrace_size; ++i)
    {
      Dl_info info;
      fprintf(stderr, "  ");
      if (dladdr(backtrace_buffer[i], &info))
      {
        if (info.dli_fname)
          fprintf(stderr, "%s", info.dli_fname);

        if (info.dli_sname)
        {
          int status = -1;
          char *dem  = abi::__cxa_demangle(info.dli_sname, 0, 0, &status);

          fprintf(stderr, "(%s+0x%lx)", status == 0 ? dem : info.dli_sname,
                  ((long)backtrace_buffer[i])-((long)info.dli_saddr));

          if (status == 0)
            free(dem);
        }
      }
      fprintf(stderr, " [0x%lx]\n", (long)backtrace_buffer[i]);
    }
  }
#else
    fputs("terminate() called, but no backtrace support\n", stderr);
#endif
}

} // anonymous namespace

const vw::TerminateHandler vw::DefaultTerminate = noop;

const vw::TerminateHandler vw::PrettyTerminate =
#ifdef USE_INTROSPECTION_EXCEPTION
  print_exception_info;
#else
  0;
#endif


const vw::TerminateHandler vw::BacktraceTerminate =
#ifdef USE_INTROSPECTION_BACKTRACE
  do_backtrace;
#else
  0;
#endif

bool vw::set_terminate_handler( vw::TerminateHandler const &th )
{
  if (real_default == 0)
  {
    real_default = std::set_terminate(std::abort);
    std::set_terminate(real_default);
  }

  if( th ) {
    current = th;
    std::set_terminate(vw_terminate);
    return true;
  }
  else
  {
    current = 0;
    std::set_terminate(real_default);
    return false;
  }
}
