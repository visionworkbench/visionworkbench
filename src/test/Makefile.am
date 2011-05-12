# __BEGIN_LICENSE__
# Copyright (C) 2006-2011 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_CARTOGRAPHY
if HAVE_PKG_GDAL
geotif_generate_SOURCES = geotif-generate.cc
geotif_generate_CPPFLAGS = @VW_CPPFLAGS@ @PKG_CARTOGRAPHY_CPPFLAGS@
geotif_generate_LDADD = @PKG_CARTOGRAPHY_LIBS@
endif
endif

noinst_PROGRAMS = geotif-generate

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@

EXTRA_DIST = Helpers.h test_main.cc

includedir = $(prefix)/include/test
include $(top_srcdir)/config/rules.mak
