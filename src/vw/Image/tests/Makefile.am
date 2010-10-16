# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_IMAGE

TestAlgorithms_SOURCES            = TestAlgorithms.cxx
TestConvolution_SOURCES           = TestConvolution.cxx
TestEdgeExtension_SOURCES         = TestEdgeExtension.cxx
TestFilter_SOURCES                = TestFilter.cxx
TestImageMath_SOURCES             = TestImageMath.cxx
TestImageResource_SOURCES         = TestImageResource.cxx
TestImageView_SOURCES             = TestImageView.cxx
TestImageViewRef_SOURCES          = TestImageViewRef.cxx
TestInterpolation_SOURCES         = TestInterpolation.cxx
TestManipulation_SOURCES          = TestManipulation.cxx
TestPerPixelAccessorViews_SOURCES = TestPerPixelAccessorViews.cxx
TestPerPixelViews_SOURCES         = TestPerPixelViews.cxx
TestPixelMath_SOURCES             = TestPixelMath.cxx
TestPixelTypes_SOURCES            = TestPixelTypes.cxx
TestStatistics_SOURCES            = TestStatistics.cxx
TestTransform_SOURCES             = TestTransform.cxx
TestMaskedPixelMath_SOURCES       = TestMaskedPixelMath.cxx
TestMaskViews_SOURCES             = TestMaskViews.cxx
TestUtilityViews_SOURCES          = TestUtilityViews.cxx
TestMaskedPixelMath2_SOURCES      = TestMaskedPixelMath2.cxx
TestMaskedImageMath_SOURCES       = TestMaskedImageMath.cxx

TESTS = TestAlgorithms TestConvolution TestEdgeExtension TestFilter        \
        TestImageMath TestImageResource TestImageView TestImageViewRef     \
        TestInterpolation TestManipulation TestPerPixelAccessorViews       \
        TestPerPixelViews TestPixelMath TestPixelTypes TestStatistics      \
        TestTransform TestMaskedPixelMath TestMaskViews TestUtilityViews   \
        TestMaskedPixelMath2 TestMaskedImageMath

#include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ $(TEST_CPPFLAGS)
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_IMAGE_LIBS@

check_PROGRAMS = $(TESTS)

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests.am
