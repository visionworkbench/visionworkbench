# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_PHOTOMETRY

include_HEADERS = Albedo.h Camres.h Exposure.h Index.h Misc.h Outlier.h \
				  Reconstruct.h Reflectance.h Shadow.h ShapeFromShading.h \
				  Shape.h Weights.h

libvwPhotometry_la_SOURCES = Albedo.cc Camres.cc Exposure.cc Index.cc Misc.cc \
						  Outlier.cc Reflectance.cc Shadow.cc Shape.cc \
						  ShapeFromShading.cc Weights.cc

lib_LTLIBRARIES = libvwPhotometry.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ @PKG_BOOST_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ -version-info @LIBTOOL_VERSION@

SUBDIRS = .

includedir = $(prefix)/include/vw/Photometry

include $(top_srcdir)/config/rules.mak
