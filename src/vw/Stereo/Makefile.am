# __BEGIN_LICENSE__
# Copyright (C) 2006-2011 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_STEREO

# Headers are referenced in the rewrite directory
# rewrite_headers =
rewrite_sources = rewrite/Correlation.cc

include_HEADERS = AffineMixtureComponent.h Correlate.h Correlate.tcc	\
        CorrelateResearch.h CorrelateResearch.tcc DisparityMap.h	\
        EMSubpixelCorrelatorView.h EMSubpixelCorrelatorView.hpp		\
        GammaMixtureComponent.h GaussianMixtureComponent.h		\
        MixtureComponent.h StereoModel.h StereoView.h SubpixelView.h \
        UniformMixtureComponent.h

libvwStereo_la_SOURCES = StereoModel.cc          \
        Correlate.cc  EMSubpixelCorrelatorView.cc \
        CorrelateResearch.cc $(rewrite_sources)

libvwStereo_la_LIBADD = @MODULE_STEREO_LIBS@

lib_LTLIBRARIES = libvwStereo.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ $(LIBTOOL_LIB_FLAGS)

SUBDIRS = . tests rewrite rewrite/tests

includedir = $(prefix)/include/vw/Stereo

include $(top_srcdir)/config/rules.mak
