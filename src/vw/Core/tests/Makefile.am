# __BEGIN_LICENSE__
# Copyright (C) 2006, 2007 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_CORE

nodist_TestCompoundTypes_SOURCES = TestCompoundTypes.cxx
nodist_TestFunctors_SOURCES = TestFunctors.cxx
nodist_TestExceptions_SOURCES = TestExceptions.cxx
nodist_TestCache_SOURCES = TestCache.cxx
nodist_TestThread_SOURCES = TestThread.cxx
nodist_TestThreadPool_SOURCES = TestThreadPool.cxx
nodist_TestSettings_SOURCES = TestSettings.cxx
nodist_TestLog_SOURCES = TestLog.cxx
nodist_TestFundamentalTypes_SOURCES = TestFundamentalTypes.cxx

TESTS = TestCompoundTypes TestFunctors TestExceptions TestThread TestThreadPool TestSettings TestCache TestLog TestFundamentalTypes

include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ $(TEST_CPPFLAGS)
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_CORE_LIBS@

check_PROGRAMS = $(TESTS)
EXTRA_DIST = $(TESTS:%=%.h)
CLEANFILES = $(TESTS:%=%.cxx) log.txt

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests_old.am
