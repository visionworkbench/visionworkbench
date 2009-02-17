# __BEGIN_LICENSE__
# Copyright (C) 2006, 2007 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

StandardShaders.cc:
	perl StandardShaders.pl

CLEANFILES = StandardShaders.cc

if MAKE_MODULE_GPU

include_HEADERS = Setup.h Utilities.h StandardShaders.h GPUProgram.h	\
  TexAlloc.h TexObj.h GPUImage.h Filter.h Manipulation.h Transform.h	\
  GenericShaders.h ImageMath.h Statistics.h Algorithms.h		\
  Interpolation.h EdgeExtension.h

libvwGPU_la_SOURCES = Setup.cc Utilities.cc StandardShaders.cc	\
  GPUProgram.cc TexAlloc.cc TexObj.cc GPUImage.cc Filter.cc	\
  Manipulation.cc Transform.cc GenericShaders.cc ImageMath.cc	\
  Statistics.cc Algorithms.cc

libvwGPU_la_LIBADD = @MODULE_GPU_LIBS@

lib_LTLIBRARIES = libvwGPU.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ -version-info @LIBTOOL_VERSION@

SUBDIRS = . tests

includedir = $(prefix)/include/vw/GPU

include $(top_srcdir)/config/rules.mak
