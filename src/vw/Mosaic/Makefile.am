# __BEGIN_LICENSE__
# Copyright (C) 2006, 2007 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_MOSAIC

include_HEADERS = QuadTreeGenerator.h ImageComposite.h	\
	KMLQuadTreeConfig.h TMSQuadTreeConfig.h		\
	UniviewQuadTreeConfig.h GMapQuadTreeConfig.h	\
	CelestiaQuadTreeConfig.h ToastQuadTreeConfig.h

libvwMosaic_la_SOURCES = QuadTreeGenerator.cc KMLQuadTreeConfig.cc	\
	TMSQuadTreeConfig.cc UniviewQuadTreeConfig.cc			\
	GMapQuadTreeConfig.cc CelestiaQuadTreeConfig.cc

libvwMosaic_la_LIBADD = @MODULE_MOSAIC_LIBS@

lib_LTLIBRARIES = libvwMosaic.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ -version-info @LIBTOOL_VERSION@

SUBDIRS = . tests

includedir = $(prefix)/include/vw/Mosaic

include $(top_srcdir)/config/rules.mak
