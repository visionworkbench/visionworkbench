# __BEGIN_LICENSE__
# Copyright (C) 2006-2011 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_STEREO

TestStereoModel_SOURCES  = TestStereoModel.cxx
TestDisparity_SOURCES    = TestDisparity.cxx
TestCorrelator_SOURCES   = TestCorrelator.cxx
TestSubPixel_SOURCES     = TestSubPixel.cxx

TESTS = TestStereoModel TestDisparity TestCorrelator TestSubPixel

#include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ $(TEST_CPPFLAGS)
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_STEREO_LIBS@

check_PROGRAMS = $(TESTS)

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests.am
