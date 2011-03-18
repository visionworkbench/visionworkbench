# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_BUNDLEADJUSTMENT

if MAKE_MODULE_INTERESTPOINT
if MAKE_MODULE_STEREO
# These deal with loading of interest points
relation_headers = CameraRelation.h ControlNetworkLoader.h

relation_sources = CameraRelation.cc ControlNetworkLoader.cc
endif
endif

include_HEADERS = BundleAdjustReport.h ControlNetwork.h ModelBase.h         \
                  AdjustBase.h AdjustRef.h AdjustRobustRef.h AdjustSparse.h \
                  AdjustRobustSparse.h $(relation_headers)

libvwBundleAdjustment_la_SOURCES = BundleAdjustReport.cc ControlNetwork.cc  \
                  $(relation_sources)

libvwBundleAdjustment_la_LIBADD = @MODULE_BUNDLEADJUSTMENT_LIBS@

lib_LTLIBRARIES = libvwBundleAdjustment.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ $(LIBTOOL_LIB_FLAGS)

SUBDIRS = . tests

includedir = $(prefix)/include/vw/BundleAdjustment

include $(top_srcdir)/config/rules.mak
