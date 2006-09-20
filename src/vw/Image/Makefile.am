# __BEGIN_LICENSE__
#
# Copyright (C) 2006 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration
# (NASA).  All Rights Reserved.
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

if MAKE_MODULE_VW

include_HEADERS = PixelIterator.h ImageViewBase.h PixelAccessors.h		\
	ImageView.h ImageViewRef.h PixelTypeInfo.h PixelMath.h PixelTypes.h	\
	GenericImageBuffer.h Manipulation.h PerPixelViews.h ImageMath.h			\
	EdgeExtend.h Interpolation.h Transform.h

libvwImage_la_SOURCES = GenericImageBuffer.cc

lib_LTLIBRARIES = libvwImage.la

endif

########################################################################
# tests
########################################################################

if MAKE_MODULE_VW

TESTS_H = tests/TestPixelMath.h tests/TestPixelTypes.h	\
	tests/TestImageView.h tests/TestImageViewRef.h				\
	tests/TestManipulation.h tests/TestPerPixelViews.h		\
	tests/TestImageMath.h tests/TestEdgeExtend.h					\
	tests/TestInterpolation.h tests/TestTransform.h

TESTS_EXTRAS = 

endif

########################################################################
# general
########################################################################

includedir = $(prefix)/include/vw/Image

AM_CPPFLAGS = @MODULE_VW_CPPFLAGS@
AM_LDFLAGS  = @MODULE_VW_LDFLAGS@

EXTRA_DIST = $(TESTS_H) $(TESTS_EXTRAS) 

include $(top_srcdir)/config/rules.mak
