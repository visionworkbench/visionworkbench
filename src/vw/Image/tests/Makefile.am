# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_IMAGE

TestPixelMath_SOURCES             = TestPixelMath.oldtest
TestPixelTypes_SOURCES            = TestPixelTypes.oldtest
TestImageView_SOURCES             = TestImageView.oldtest
TestImageViewRef_SOURCES          = TestImageViewRef.oldtest
TestManipulation_SOURCES          = TestManipulation.oldtest
TestPerPixelViews_SOURCES         = TestPerPixelViews.oldtest
TestPerPixelAccessorViews_SOURCES = TestPerPixelAccessorViews.oldtest
TestImageMath_SOURCES             = TestImageMath.oldtest
TestEdgeExtension_SOURCES         = TestEdgeExtension.oldtest
TestInterpolation_SOURCES         = TestInterpolation.oldtest
TestConvolution_SOURCES           = TestConvolution.oldtest
TestFilter_SOURCES                = TestFilter.oldtest
TestTransform_SOURCES             = TestTransform.oldtest
TestAlgorithms_SOURCES            = TestAlgorithms.oldtest
TestStatistics_SOURCES            = TestStatistics.oldtest

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
CLEANFILES = $(TESTS:%=%.cxx)

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests_old.am
