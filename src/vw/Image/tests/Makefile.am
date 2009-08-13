# __BEGIN_LICENSE__
# Copyright (C) 2006, 2007 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_IMAGE

nodist_TestPixelMath_SOURCES = TestPixelMath.cxx
nodist_TestPixelTypes_SOURCES = TestPixelTypes.cxx
nodist_TestImageView_SOURCES = TestImageView.cxx
nodist_TestImageViewRef_SOURCES = TestImageViewRef.cxx
nodist_TestManipulation_SOURCES = TestManipulation.cxx
nodist_TestPerPixelViews_SOURCES = TestPerPixelViews.cxx
nodist_TestPerPixelAccessorViews_SOURCES = TestPerPixelAccessorViews.cxx
nodist_TestImageMath_SOURCES = TestImageMath.cxx
nodist_TestEdgeExtension_SOURCES = TestEdgeExtension.cxx
nodist_TestInterpolation_SOURCES = TestInterpolation.cxx
nodist_TestConvolution_SOURCES = TestConvolution.cxx
nodist_TestFilter_SOURCES = TestFilter.cxx
nodist_TestTransform_SOURCES = TestTransform.cxx
nodist_TestAlgorithms_SOURCES = TestAlgorithms.cxx
nodist_TestStatistics_SOURCES = TestStatistics.cxx

TESTS = TestPixelMath TestPixelTypes TestImageView TestImageViewRef		\
	TestManipulation TestPerPixelViews TestPerPixelAccessorViews		\
	TestImageMath TestEdgeExtension TestInterpolation TestConvolution	\
	TestFilter TestTransform TestAlgorithms TestStatistics

include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ $(TEST_CPPFLAGS)
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_IMAGE_LIBS@

check_PROGRAMS = $(TESTS)
EXTRA_DIST = $(TESTS:%=%.h)
CLEANFILES = $(TESTS:%=%.cxx)

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests_old.am
