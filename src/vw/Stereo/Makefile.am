# __BEGIN_LICENSE__
# Copyright (C) 2006-2011 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_STEREO

include_HEADERS = AffineMixtureComponent.h Algorithms.h Correlate.h	\
        Correlate.tcc Correlation.h CorrelationView.h			\
        CorrelateResearch.h CorrelateResearch.tcc CostFunctions.h	\
        DisparityMap.h EMSubpixelCorrelatorView.h			\
        EMSubpixelCorrelatorView.hpp GammaMixtureComponent.h		\
        GaussianMixtureComponent.h MixtureComponent.h PreFilter.h	\
        StereoModel.h StereoView.h SubpixelView.h			\
        UniformMixtureComponent.h

libvwStereo_la_SOURCES = StereoModel.cc Correlate.cc Correlation.cc	\
        EMSubpixelCorrelatorView.cc CorrelateResearch.cc

libvwStereo_la_LIBADD = @MODULE_STEREO_LIBS@

lib_LTLIBRARIES = libvwStereo.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ $(LIBTOOL_LIB_FLAGS)

SUBDIRS = . tests

includedir = $(prefix)/include/vw/Stereo

include $(top_srcdir)/config/rules.mak
