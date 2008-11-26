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
                  $(gdal_headers) $(camerabbox_headers)

libvwCartography_la_SOURCES = Datum.cc GeoReference.cc GeoTransform.cc	\
                  GeoReferenceResourcePDS.cc $(gdal_sources)		\
                  $(camerabbox_sources)

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
