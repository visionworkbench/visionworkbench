# __BEGIN_LICENSE__
# Copyright (C) 2006, 2007 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_CARTOGRAPHY

if HAVE_PKG_GDAL
gdal_headers = GeoReferenceResourceGDAL.h
gdal_sources = GeoReferenceResourceGDAL.cc
endif

if MAKE_MODULE_CAMERA
camerabbox_headers = CameraBBox.h
camerabbox_sources = CameraBBox.cc
endif

include_HEADERS = FileIO.h GeoReferenceBase.h GeoReference.h		\
                  GeoTransform.h Datum.h PointImageManipulation.h	\
                  OrthoImageView.h GeoReferenceResourcePDS.h		\
                  Projection.h ToastTransform.h $(gdal_headers)		\
                  $(camerabbox_headers)

libvwCartography_la_SOURCES = Datum.cc GeoReference.cc GeoTransform.cc	\
                  GeoReferenceResourcePDS.cc ToastTransform.cc		\
                  $(gdal_sources) $(camerabbox_sources)

libvwCartography_la_LIBADD = @MODULE_CARTOGRAPHY_LIBS@

lib_LTLIBRARIES = libvwCartography.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ -version-info @LIBTOOL_VERSION@

SUBDIRS = . tests

includedir = $(prefix)/include/vw/Cartography

include $(top_srcdir)/config/rules.mak
