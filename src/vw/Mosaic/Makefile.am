# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_MOSAIC

include_HEADERS = \
  CelestiaQuadTreeConfig.h \
  GigapanQuadTreeConfig.h \
  GMapQuadTreeConfig.h \
  ImageComposite.h \
  KMLQuadTreeConfig.h \
  QuadTreeConfig.h \
  QuadTreeGenerator.h \
  TMSQuadTreeConfig.h \
  ToastQuadTreeConfig.h \
  UniviewQuadTreeConfig.h

libvwMosaic_la_SOURCES = \
  CelestiaQuadTreeConfig.cc \
  GigapanQuadTreeConfig.cc \
  GMapQuadTreeConfig.cc \
  KMLQuadTreeConfig.cc \
  QuadTreeConfig.cc \
  QuadTreeGenerator.cc \
  TMSQuadTreeConfig.cc \
  UniviewQuadTreeConfig.cc

libvwMosaic_la_LIBADD = @MODULE_MOSAIC_LIBS@

lib_LTLIBRARIES = libvwMosaic.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ $(LIBTOOL_LIB_FLAGS)

SUBDIRS = . tests

includedir = $(prefix)/include/vw/Mosaic

include $(top_srcdir)/config/rules.mak
