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
gdal_headers = DiskImageResourceGeoReferenceHelperGDAL.h
gdal_sources = DiskImageResourceGeoReferenceHelperGDAL.cc
endif

if HAVE_PKG_JPEG2K
jpeg2k_headers = DiskImageResourceGeoReferenceHelperJP2.h
jpeg2k_sources = DiskImageResourceGeoReferenceHelperJP2.cc
endif

include_HEADERS = FileIO.h GeoReferenceBase.h GeoReference.h			\
                  GeoTransform.h Datum.h PointImageManipulation.h	\
                  OrthoImageView.h																\
                  DiskImageResourceGeoReferenceHelperPDS.h				\
                  $(gdal_headers) $(jpeg2k_headers)

libvwCartography_la_SOURCES = Datum.cc GeoReference.cc GeoTransform.cc	\
                  DiskImageResourceGeoReferenceHelperPDS.cc		\
                  $(gdal_sources) $(jpeg2k_sources)

libvwCartography_la_LDFLAGS= -version-info @LIBTOOL_VERSION@
libvwCartography_la_LIBADD = @MODULE_CARTOGRAPHY_LDFLAGS@

lib_LTLIBRARIES = libvwCartography.la

if ENABLE_EXCEPTIONS
bin_PROGRAMS = shaded_relief
shaded_relief_SOURCES = shaded_relief.cc
shaded_relief_LDFLAGS = @PKG_CARTOGRAPHY_LDFLAGS@ @PKG_BOOST_PROGRAM_OPTIONS_LDFLAGS@
endif

endif

########################################################################
# general
########################################################################

SUBDIRS = . tests

includedir = $(prefix)/include/vw/Cartography

AM_CPPFLAGS = @MODULE_CARTOGRAPHY_CPPFLAGS@
AM_LDFLAGS  = 

include $(top_srcdir)/config/rules.mak
