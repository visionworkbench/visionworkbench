# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_IMAGE

include_HEADERS = PixelIterator.h ImageViewBase.h PixelAccessors.h	\
	ImageView.h ImageIO.h ImageViewRef.h PixelTypeInfo.h		\
	PixelMath.h PixelTypes.h PixelMask.h ImageResource.h		\
	ImageResourceView.h ViewImageResource.h Manipulation.h		\
	PerPixelViews.h PerPixelAccessorViews.h ImageMath.h		\
	EdgeExtension.h Interpolation.h Convolution.h Filter.h		\
	Filter.tcc Transform.h Statistics.h Algorithms.h EdgeExtend.h	\
	BlockProcessor.h BlockRasterize.h UtilityViews.h Palette.h	\
	MaskViews.h SparseImageCheck.h

libvwImage_la_SOURCES = PixelTypeInfo.cc ImageResource.cc Filter.cc     \
	Interpolation.cc

libvwImage_la_LIBADD = @MODULE_IMAGE_LIBS@

lib_LTLIBRARIES = libvwImage.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ -version-info @LIBTOOL_VERSION@

SUBDIRS = . tests

includedir = $(prefix)/include/vw/Image

include $(top_srcdir)/config/rules.mak
