# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_STEREO

TestAlgorithms_SOURCES      = TestAlgorithms.cxx
TestPreFilter_SOURCES       = TestPreFilter.cxx
TestCostFunctions_SOURCES   = TestCostFunctions.cxx
TestCorrelation_SOURCES     = TestCorrelation.cxx
TestCorrelationView_SOURCES = TestCorrelationView.cxx
TestPyramidCorrelationView_SOURCES = TestPyramidCorrelationView.cxx

TESTS = TestAlgorithms TestPreFilter TestCostFunctions			\
  TestCorrelation TestCorrelationView TestPyramidCorrelationView

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
