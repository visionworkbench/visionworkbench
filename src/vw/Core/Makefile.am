# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_CORE

include_HEADERS = \
  Cache.h \
  CompoundTypes.h \
  ConfigParser.h \
  Debugging.h \
  Exception.h \
  Features.h \
  Functors.h \
  FundamentalTypes.h \
  Log.h \
  ProgressCallback.h \
  Settings.h \
  Stopwatch.h \
  System.h \
  Thread.h \
  ThreadPool.h \
  ThreadQueue.h \
  TypeDeduction.h \
  VarArray.h

libvwCore_la_SOURCES = \
  Cache.cc \
  ConfigParser.cc \
  Debugging.cc \
  Exception.cc \
  Log.cc \
  ProgressCallback.cc \
  Settings.cc \
  Stopwatch.cc \
  System.cc \
  Thread.cc

libvwCore_la_LIBADD = @MODULE_CORE_LIBS@

lib_LTLIBRARIES = libvwCore.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ @PKG_BOOST_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ $(LIBTOOL_LIB_FLAGS)

SUBDIRS = . tests

includedir = $(prefix)/include/vw/Core

include $(top_srcdir)/config/rules.mak
