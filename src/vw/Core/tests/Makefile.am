# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_CORE

TestCache_SOURCES            = TestCache.cxx
TestCompoundTypes_SOURCES    = TestCompoundTypes.cxx
TestExceptions_SOURCES       = TestExceptions.cxx
TestFunctors_SOURCES         = TestFunctors.cxx
TestFundamentalTypes_SOURCES = TestFundamentalTypes.cxx
TestLog_SOURCES              = TestLog.cxx
TestSettings_SOURCES         = TestSettings.cxx
TestThreadPool_SOURCES       = TestThreadPool.cxx
TestThreadQueue_SOURCES      = TestThreadQueue.cxx
TestThread_SOURCES           = TestThread.cxx
TestTypeDeduction_SOURCES    = TestTypeDeduction.cxx

TESTS = \
  TestCache \
  TestCompoundTypes \
  TestExceptions \
  TestFunctors \
  TestFundamentalTypes \
  TestLog \
  TestSettings \
  TestThread \
  TestThreadPool \
  TestThreadQueue \
  TestTypeDeduction

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_CORE_LIBS@

check_PROGRAMS = $(TESTS)
CLEANFILES = log.txt

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests.am
