# __BEGIN_LICENSE__
# Copyright (C) 2006, 2007 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_CORE

TestCompoundTypes_SOURCES    = TestCompoundTypes.cxx
TestFunctors_SOURCES         = TestFunctors.cxx
TestExceptions_SOURCES       = TestExceptions.cxx
TestThread_SOURCES           = TestThread.cxx
TestThreadPool_SOURCES       = TestThreadPool.cxx
TestSettings_SOURCES         = TestSettings.cxx
TestCache_SOURCES            = TestCache.cxx
TestLog_SOURCES              = TestLog.cxx
TestFundamentalTypes_SOURCES = TestFundamentalTypes.cxx

TESTS = TestCompoundTypes TestFunctors TestExceptions TestThread TestThreadPool TestSettings TestCache TestLog TestFundamentalTypes

include $(top_srcdir)/config/instantiate.am

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
