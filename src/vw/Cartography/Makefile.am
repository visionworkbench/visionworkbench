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
gdal_headers = GeoReferenceResourceGDAL.h
gdal_sources = GeoReferenceResourceGDAL.cc
endif

if MAKE_MODULE_CAMERA
camerabbox_headers = CameraBBox.h
camerabbox_sources = CameraBBox.cc
endif

if HAVE_PKG_PROTOBUF
protocol_headers = DatumDesc.pb.h GeoReferenceDesc.pb.h
protocol_sources = DatumDesc.pb.cc GeoReferenceDesc.pb.cc
endif

include_HEADERS = GeoReferenceBase.h GeoReference.h                     \
                  GeoTransform.h Datum.h SimplePointImageManipulation.h \
                  PointImageManipulation.h                              \
                  OrthoImageView.h GeoReferenceResourcePDS.h            \
                  Projection.h ToastTransform.h $(gdal_headers)         \
                  $(camerabbox_headers)

detail_HEADERS =  detail/BresenhamLine.h

libvwCartography_la_SOURCES = Datum.cc GeoReference.cc GeoTransform.cc  \
                  GeoReferenceResourcePDS.cc ToastTransform.cc          \
                  $(gdal_sources) $(camerabbox_sources)

nodist_libvwCartography_la_SOURCES = $(protocol_sources)

libvwCartography_la_LIBADD = @MODULE_CARTOGRAPHY_LIBS@

lib_LTLIBRARIES = libvwCartography.la

endif

########################################################################
# general
########################################################################

BUILT_SOURCES = $(protocol_sources)

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ $(LIBTOOL_LIB_FLAGS)

SUBDIRS = . tests
CLEANFILES = $(protocol_headers) $(protocol_sources)
EXTRA_DIST = DatumDesc.proto GeoReferenceDesc.proto $(protocol_headers) $(protocol_sources)

includedir = $(prefix)/include/vw/Cartography
detaildir = $(includedir)/detail

include $(top_srcdir)/config/rules.mak
if HAVE_PKG_PROTOBUF
include $(top_srcdir)/thirdparty/protobuf.mak
endif

