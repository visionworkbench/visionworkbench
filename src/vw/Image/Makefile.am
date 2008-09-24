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

include_HEADERS = PixelIterator.h ImageViewBase.h PixelAccessors.h	\
	ImageView.h ImageIO.h ImageViewRef.h PixelTypeInfo.h		\
	PixelMath.h PixelTypes.h PixelMask.h ImageResource.h		\
	ImageResourceView.h Manipulation.h PerPixelViews.h		\
	PerPixelAccessorViews.h ImageMath.h EdgeExtension.h		\
	Interpolation.h Convolution.h Filter.h Filter.tcc Transform.h	\
	Statistics.h Algorithms.h EdgeExtend.h BlockProcessor.h		\
	BlockRasterize.h BlockCacheView.h UtilityViews.h Palette.h

libvwImage_la_SOURCES = PixelTypeInfo.cc ImageResource.cc Filter.cc
libvwImage_la_LIBADD = @MODULE_IMAGE_LIBS@

lib_LTLIBRARIES = libvwImage.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ -version-info @LIBTOOL_VERSION@ -no-undefined

SUBDIRS = . tests

includedir = $(prefix)/include/vw/Image

include $(top_srcdir)/config/rules.mak
