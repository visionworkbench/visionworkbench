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
core_subdirs = Core
core_headers = Core.h
endif

if MAKE_MODULE_MATH
math_subdirs = Math
math_headers = Math.h
endif

if MAKE_MODULE_IMAGE
image_subdirs = Image
image_headers = Image.h
endif

if MAKE_MODULE_FILEIO
fileio_subdirs = FileIO
fileio_headers = FileIO.h
endif

if MAKE_MODULE_CAMERA
camera_subdirs = Camera
camera_headers = Camera.h
endif

if MAKE_MODULE_MOSAIC
mosaic_subdirs = Mosaic
mosaic_headers = Mosaic.h
endif

if MAKE_MODULE_STEREO
stereo_subdirs = Stereo
stereo_headers = Stereo.h
endif

if MAKE_MODULE_CARTOGRAPHY
cartography_subdirs = Cartography
cartography_headers = Cartography.h
endif

include_HEADERS = $(core_headers) $(math_headers) $(image_headers)	\
	$(fileio_headers) $(camera_headers) $(mosaic_headers) $(stereo_headers) \
  $(cartography_headers) vw.h

SUBDIRS = $(core_subdirs) $(math_subdirs) $(image_subdirs)	\
	$(fileio_subdirs) . $(camera_subdirs) $(mosaic_subdirs) $(stereo_subdirs) \
  $(cartography_subdirs) tests tools

lib_LTLIBRARIES = libvw.la
libvw_la_SOURCES =
libvw_la_LIBADD = Core/libvwCore.la Math/libvwMath.la Image/libvwImage.la FileIO/libvwFileIO.la

########################################################################
# general
########################################################################

includedir = $(prefix)/include/vw

include $(top_srcdir)/config/rules.mak
