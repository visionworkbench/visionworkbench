# __BEGIN_LICENSE__
# Copyright (C) 2006, 2007 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_INTERESTPOINT

nodist_TestMatcher_SOURCES = TestMatcher.cxx 
#nodist_TestIntegral_SOURCES = TestIntegral.cxx

TESTS = TestMatcher #TestIntegral

include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ -I$(top_srcdir)/thirdparty/cxxtest
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_INTERESTPOINT_LIBS@ @PKG_VW_LIBS@

check_PROGRAMS = $(TESTS)
EXTRA_DIST = $(TESTS:%=%.h) sub.png noisy_gradient_60.png
CLEANFILES = $(TESTS:%=%.cxx)

include $(top_srcdir)/config/rules.mak
