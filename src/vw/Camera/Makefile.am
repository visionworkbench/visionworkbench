# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_CAMERA

protocol_headers =
protocol_sources =
include_HEADERS =
libvwCamera_la_SOURCES =

if HAVE_PKG_PROTOBUF
protocol_headers += TsaiFile.pb.h
protocol_sources += TsaiFile.pb.cc
BUILT_SOURCES   = $(protocol_sources)
include_HEADERS += $(protocol_headers)
endif

if HAVE_PKG_LAPACK
include_HEADERS += Extrinsics.h LinescanModel.h LinearPushbroomModel.h   \
                   OrbitingPushbroomModel.h PinholeModelCalibrate.h      \
                   LensDistortion.h CameraGeometry.h
libvwCamera_la_SOURCES += Extrinsics.cc
endif

include_HEADERS += CameraModel.h CAHVModel.h CAHVORModel.h               \
                   CameraTransform.h Exif.h ExifData.h BayerFilter.h     \
                   PinholeModel.h

libvwCamera_la_SOURCES += PinholeModel.cc CAHVModel.cc CAHVORModel.cc    \
              Exif.cc ExifData.cc LensDistortion.cc $(protocol_sources)  \
              CameraModel.cc

libvwCamera_la_LIBADD = @MODULE_CAMERA_LIBS@

lib_LTLIBRARIES = libvwCamera.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ $(LIBTOOL_LIB_FLAGS)

SUBDIRS = . tests
CLEANFILES = $(protocol_headers) $(protocol_sources)
EXTRA_DIST = TsaiFile.proto

includedir = $(prefix)/include/vw/Camera

include $(top_srcdir)/config/rules.mak
if HAVE_PKG_PROTOBUF
include $(top_srcdir)/thirdparty/protobuf.mak
endif

