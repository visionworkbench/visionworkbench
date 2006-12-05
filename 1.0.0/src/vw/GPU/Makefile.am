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

if MAKE_MODULE_GPU

include_HEADERS = OpenGLManager.h GPUProgram.h TexAlloc.h TexBlock.h	\
	  TexObj.h GPUImage.h TokenReplacer.h Filter.h			\
	  ImageOperators.h SimpleTransforms.h Transform.h Display.h	\
	  GenericShaders.h ImageMath.h ImageStatistics.h Fourier.h	\
	  ImageAlgorithms.h Interpolation.h EdgeExtend.h Expressions.h	\
	  Timer.h

libvwGPU_la_SOURCES = OpenGLManager.cc GPUProgram.cc TexAlloc.cc	\
	TexBlock.cc TexObj.cc GPUImage.cc TokenReplacer.cc		\
	GenericShaders.cc ImageOperators.cc ImageMath.cc		\
	SimpleTransforms.cc Transform.cc Filter.cc ImageStatistics.cc	\
	Fourier.cc Display.cc ImageAlgorithms.cc Expressions.cc		\
	UtilityTimer.cc

libvwGPU_la_LDFLAGS= -version-info @LIBTOOL_VERSION@
libvwGPU_la_LIBADD = @MODULE_GPU_LDFLAGS@

lib_LTLIBRARIES = libvwGPU.la

bin_PROGRAMS = 

endif

########################################################################
# general
########################################################################

SUBDIRS = . tests

includedir = $(prefix)/include/vw/GPU

AM_CPPFLAGS = @MODULE_GPU_CPPFLAGS@
AM_LDFLAGS  = 

include $(top_srcdir)/config/rules.mak
