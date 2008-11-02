# __BEGIN_LICENSE__
# 
# Copyright (C) 2006 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration
# (NASA).  All Rights Reserved.
# 
# Copyright 2006 Carnegie Mellon University. All rights reserved.
# 
# This software is distributed under the NASA Open Source Agreement
# (NOSA), version 1.3.  The NOSA has been approved by the Open Source
# Initiative.  See the file COPYING at the top of the distribution
# directory tree for the complete NOSA document.
# 
# THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
# KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
# LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
# SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
# A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
# THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
# DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
# 
# __END_LICENSE__

########################################################################
# sources
########################################################################

if MAKE_MODULE_IMAGE

TESTS = TestPixelMath TestPixelTypes TestImageView TestImageViewRef	\
	TestManipulation TestPerPixelViews TestPerPixelAccessorViews			\
	TestImageMath TestEdgeExtension TestInterpolation TestConvolution	\
	TestFilter TestTransform TestAlgorithms TestStatistics

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

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ -I$(top_srcdir)/thirdparty/cxxtest
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_IMAGE_LIBS@ @PKG_MATH_LIBS@ @PKG_CORE_LIBS@

check_PROGRAMS = $(TESTS)
EXTRA_DIST = $(TESTS:%=%.h)
CLEANFILES = $(TESTS:%=%.cxx)

include $(top_srcdir)/config/rules.mak
