# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_STEREO

include_HEADERS = Correlate.h ReferenceCorrelator.h             \
        OptimizedCorrelator.h DisparityMap.h StereoModel.h      \
        PyramidCorrelator.h StereoView.h CorrelatorView.h       \
        SubpixelView.h EMSubpixelCorrelatorView.h               \
        MixtureComponent.h GammaMixtureComponent.h              \
        GaussianMixtureComponent.h                              \
        AffineMixtureComponent.h UniformMixtureComponent.h      \
        EMSubpixelCorrelatorView.hpp CorrelateResearch.h        \
        Correlate.tcc CorrelateResearch.tcc

libvwStereo_la_SOURCES = StereoModel.cc PyramidCorrelator.cc            \
        Correlate.cc OptimizedCorrelator.cc EMSubpixelCorrelatorView.cc \
        CorrelateResearch.cc

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
