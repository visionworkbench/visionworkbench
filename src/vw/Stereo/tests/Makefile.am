# __BEGIN_LICENSE__
# Copyright (C) 2006-2011 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_STEREO

TestAlgorithms_SOURCES    = TestAlgorithms.cxx
TestCorrelate_SOURCES     = TestCorrelate.cxx
TestCorrelationView_SOURCES = TestCorrelationView.cxx
TestCorrelation_SOURCES   = TestCorrelation.cxx
TestCostFunctions_SOURCES = TestCostFunctions.cxx
TestDisparity_SOURCES     = TestDisparity.cxx
TestPreFilter_SOURCES     = TestPreFilter.cxx
TestPyramidCorrelationView_SOURCES = TestPyramidCorrelationView.cxx
TestStereoModel_SOURCES   = TestStereoModel.cxx
TestSubPixel_SOURCES      = TestSubPixel.cxx

TESTS = TestAlgorithms TestCorrelate TestCorrelationView	\
	TestCostFunctions TestDisparity TestPreFilter		\
	TestPyramidCorrelationView TestStereoModel TestSubPixel

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
