# __BEGIN_LICENSE__
# Copyright (C) 2006, 2007 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_CAMERA

if HAVE_PKG_LAPACK
lsm_headers = Extrinsics.h LinescanModel.h LinearPushbroomModel.h	\
              OrbitingPushbroomModel.h PinholeModel.h PinholeModelCalibrate.h \
              LensDistortion.h BundleAdjust.h
endif

include_HEADERS = CameraModel.h CAHVModel.h CAHVORModel.h		\
                  CameraTransform.h Exif.h ExifData.h BayerFilter.h	\
                  ControlNetwork.h $(lsm_headers)                       \
                  BundleAdjustReport.h 

libvwCamera_la_SOURCES = PinholeModel.cc CAHVModel.cc CAHVORModel.cc	\
                         Exif.cc ExifData.cc ControlNetwork.cc

libvwCamera_la_LIBADD = @MODULE_CAMERA_LIBS@

lib_LTLIBRARIES = libvwCamera.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ -version-info @LIBTOOL_VERSION@

SUBDIRS = . tests

includedir = $(prefix)/include/vw/Camera

include $(top_srcdir)/config/rules.mak

