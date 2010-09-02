# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_HDR

include_HEADERS = CameraCurve.h GlobalToneMap.h LDRtoHDR.h	\
	LocalToneMap.h

libvwHDR_la_SOURCES = LocalToneMap.cc CameraCurve.cc

libvwHDR_la_LIBADD = @MODULE_HDR_LIBS@

lib_LTLIBRARIES = libvwHDR.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ $(LIBTOOL_LIB_FLAGS)

SUBDIRS = . tests

includedir = $(prefix)/include/vw/HDR

include $(top_srcdir)/config/rules.mak

