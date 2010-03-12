# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_CORE

include_HEADERS = Exception.h FundamentalTypes.h TypeDeduction.h	\
	VarArray.h Functors.h CompoundTypes.h Debugging.h Thread.h	\
	ThreadPool.h Cache.h ProgressCallback.h Stopwatch.h Settings.h	\
	Log.h ConfigParser.h ThreadQueue.h Features.h

libvwCore_la_SOURCES = Debugging.cc Exception.cc Thread.cc Cache.cc	\
	ProgressCallback.cc Stopwatch.cc Settings.cc Log.cc		\
	ConfigParser.cc
libvwCore_la_LIBADD = @MODULE_CORE_LIBS@

lib_LTLIBRARIES = libvwCore.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ @PKG_BOOST_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ -version-info @LIBTOOL_VERSION@

SUBDIRS = . tests

includedir = $(prefix)/include/vw/Core

include $(top_srcdir)/config/rules.mak
