# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_INTERESTPOINT

TestMatcher_SOURCES  = TestMatcher.oldtest
TestIntegral_SOURCES = TestIntegral.oldtest
TestBoxFilter_SOURCES = TestBoxFilter.oldtest

TESTS = TestMatcher TestIntegral TestBoxFilter

include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ $(TEST_CPPFLAGS)
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_INTERESTPOINT_LIBS@

check_PROGRAMS = $(TESTS)
EXTRA_DIST = sub.png noisy_gradient_60.png
CLEANFILES = $(TESTS:%=%.cxx)

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests_old.am
