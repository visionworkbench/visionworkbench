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

if MAKE_MODULE_HDR

include_HEADERS = CameraCurve.h GlobalToneMap.h LDRtoHDR.h	\
	LocalToneMap.h PolynomialCameraCurve.h

libvwHDR_la_SOURCES = GlobalToneMap.cc LocalToneMap.cc CameraCurve.cc	\
	PolynomialCameraCurve.cc

libvwHDR_la_LDFLAGS= -version-info @LIBTOOL_VERSION@
libvwHDR_la_LIBADD = @MODULE_HDR_LDFLAGS@

lib_LTLIBRARIES = libvwHDR.la

if ENABLE_EXCEPTIONS
bin_PROGRAMS = hdr_tonemap hdr_merge
hdr_merge_SOURCES = hdr_merge.cc 
hdr_merge_LDFLAGS = @PKG_HDR_LDFLAGS@
hdr_tonemap_SOURCES = hdr_tonemap.cc 
hdr_tonemap_LDFLAGS = @PKG_HDR_LDFLAGS@ 
endif

endif

########################################################################
# general
########################################################################

SUBDIRS = . tests

includedir = $(prefix)/include/vw/HDR

AM_CPPFLAGS = @MODULE_HDR_CPPFLAGS@
AM_LDFLAGS  = 

include $(top_srcdir)/config/rules.mak

