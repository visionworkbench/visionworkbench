# __BEGIN_LICENSE__
# Copyright (C) 2006-2011 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_STEREO

rewrite_headers = Algorithms.h Correlation.h CostFunctions.h  \
        CorrelationView.h PreFilter.h SubpixelView.h

rewrite_sources = Correlation.cc

include_HEADERS = $(rewrite_headers)

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ $(LIBTOOL_LIB_FLAGS)

SUBDIRS = .

includedir = $(prefix)/include/vw/Stereo/rewrite

include $(top_srcdir)/config/rules.mak
