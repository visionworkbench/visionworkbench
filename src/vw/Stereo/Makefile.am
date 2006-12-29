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

if MAKE_MODULE_STEREO

include_HEADERS = Correlate.h ReferenceCorrelator.h		\
	OptimizedCorrelator.h MultiresolutionCorrelator.h	\
	DisparityMap.h StereoModel.h SoftwareRenderer.h		\
	PyramidCorrelator.h BlockCorrelator.h

libvwStereo_la_SOURCES = StereoModel.cc OptimizedCorrelator.cc SoftwareRenderer.cc

libvwStereo_la_LDFLAGS= -version-info @LIBTOOL_VERSION@
libvwStereo_la_LIBADD = @MODULE_STEREO_LDFLAGS@

lib_LTLIBRARIES = libvwStereo.la

if ENABLE_EXCEPTIONS
bin_PROGRAMS = correlate
correlate_SOURCES = correlate.cc
correlate_LDFLAGS = @PKG_STEREO_LDFLAGS@ @PKG_BOOST_PROGRAM_OPTIONS_LDFLAGS@
endif

endif

########################################################################
# general
########################################################################

SUBDIRS = . tests

includedir = $(prefix)/include/vw/Stereo

AM_CPPFLAGS = @MODULE_STEREO_CPPFLAGS@
AM_LDFLAGS  = 

include $(top_srcdir)/config/rules.mak
