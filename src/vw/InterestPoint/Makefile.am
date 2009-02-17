# __BEGIN_LICENSE__
# Copyright (C) 2006, 2007 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_INTERESTPOINT

include_HEADERS = Detector.h Descriptor.h Matcher.h Extrema.h	\
                  Localize.h InterestOperator.h WeightedHistogram.h	\
                  ImageOctave.h InterestData.h ImageOctaveHistory.h	\
                  InterestTraits.h MatrixIO.h VectorIO.h LearnPCA.h	\
		  SURF.h IntegralImage.h

libvwInterestPoint_la_SOURCES = InterestData.cc SURF.cc
libvwInterestPoint_la_LIBADD = @MODULE_INTERESTPOINT_LIBS@

lib_LTLIBRARIES = libvwInterestPoint.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ -version-info @LIBTOOL_VERSION@

SUBDIRS = . tests

includedir = $(prefix)/include/vw/InterestPoint

include $(top_srcdir)/config/rules.mak
