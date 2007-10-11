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

StandardShaders.cc:
	perl StandardShaders.pl

if MAKE_MODULE_GPU

include_HEADERS = Setup.h Utilities.h StandardShaders.h   \
  GPUProgram.h TexAlloc.h	\
  TexObj.h GPUImage.h  Filter.h 				\
  Manipulation.h Transform.h  GenericShaders.h						\
  ImageMath.h Statistics.h  Algorithms.h						\
  Interpolation.h EdgeExtension.h 

libvwGPU_la_SOURCES = Setup.cc Utilities.cc StandardShaders.cc 	\
  GPUProgram.cc TexAlloc.cc				\
  TexObj.cc GPUImage.cc Filter.cc 	\
  Manipulation.cc Transform.cc GenericShaders.cc			\
  ImageMath.cc Statistics.cc Algorithms.cc								

libvwGPU_la_LIBADD = @MODULE_GPU_LIBS@

lib_LTLIBRARIES = libvwGPU.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@

SUBDIRS = . tests

includedir = $(prefix)/include/vw/GPU

include $(top_srcdir)/config/rules.mak
