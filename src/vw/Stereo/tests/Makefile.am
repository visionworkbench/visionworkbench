# __BEGIN_LICENSE__
# Copyright (C) 2006, 2007 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_STEREO

nodist_TestStereoModel_SOURCES = TestStereoModel.cxx
nodist_TestDisparity_SOURCES = TestDisparity.cxx
# nodist_TestCorrelator_SOURCES = TestCorrelator.cxx

TESTS = TestStereoModel TestDisparity # TestCorrelator

include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ $(TEST_CPPFLAGS)
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_STEREO_LIBS@

check_PROGRAMS = $(TESTS)
EXTRA_DIST = $(TESTS:%=%.h)
CLEANFILES = $(TESTS:%=%.cxx)

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests_old.am
