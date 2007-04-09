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

if MAKE_MODULE_CORE

include_HEADERS = Exception.h FundamentalTypes.h TypeDeduction.h	\
	Functors.h CompoundTypes.h Debugging.h Cache.h ProgressCallback.h

libvwCore_la_SOURCES = Debugging.cc Exception.cc Cache.cc ProgressCallback.cc
libvwCore_la_LDFLAGS= -version-info @LIBTOOL_VERSION@
libvwCore_la_LIBADD = @MODULE_CORE_LDFLAGS@

lib_LTLIBRARIES = libvwCore.la

endif

########################################################################
# general
########################################################################

SUBDIRS = . tests

includedir = $(prefix)/include/vw/Core

AM_CPPFLAGS = @MODULE_CORE_CPPFLAGS@
AM_LDFLAGS  =

include $(top_srcdir)/config/rules.mak
