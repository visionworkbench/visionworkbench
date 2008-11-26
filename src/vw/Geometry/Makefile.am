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

if MAKE_MODULE_GEOMETRY

# Convex can use either PPL or Apron
if HAVE_PKG_PPL
if HAVE_PKG_GMP
convex_headers = Convex.h
convex_sources = Convex.cc
endif
endif
if HAVE_PKG_APRON
if HAVE_PKG_GMP
if HAVE_PKG_MPFR
convex_headers = Convex.h
convex_sources = Convex.cc
endif
endif
endif

include_HEADERS = Shape.h SpatialTree.h PointListIO.h Sphere.h Box.h	\
	$(convex_headers)

libvwGeometry_la_SOURCES = SpatialTree.cc $(convex_sources)
libvwGeometry_la_LIBADD = @MODULE_GEOMETRY_LIBS@

lib_LTLIBRARIES = libvwGeometry.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ -version-info @LIBTOOL_VERSION@

SUBDIRS = . tests

includedir = $(prefix)/include/vw/Geometry

include $(top_srcdir)/config/rules.mak
