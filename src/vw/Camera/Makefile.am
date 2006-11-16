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

if MAKE_MODULE_CAMERA

if HAVE_PKG_LAPACK
opm_headers = OrbitingPushbroomModel.h
opm_sources = OrbitingPushbroomModel.cc
endif

include_HEADERS = CameraModel.h CAHVModel.h CAHVORModel.h LinearPushbroomModel.h \
                  CameraTransform.h Exif.h ExifData.h $(opm_headers)

libvwCamera_la_SOURCES = CAHVModel.cc CAHVORModel.cc LinearPushbroomModel.cc \
                         Exif.cc ExifData.cc $(opm_sources)

libvwCamera_la_LIBADD = @MODULE_CAMERA_LDFLAGS@

lib_LTLIBRARIES = libvwCamera.la

bin_PROGRAMS = print_exif fisheye
print_exif_SOURCES = print_exif.cc 
print_exif_LDFLAGS = @PKG_CAMERA_LDFLAGS@
fisheye_SOURCES = fisheye.cc 
fisheye_LDFLAGS = @PKG_CAMERA_LDFLAGS@ 

endif

########################################################################
# general
########################################################################

SUBDIRS = . tests

includedir = $(prefix)/include/vw/Camera

AM_CPPFLAGS = @MODULE_CAMERA_CPPFLAGS@
AM_LDFLAGS  = 

include $(top_srcdir)/config/rules.mak

