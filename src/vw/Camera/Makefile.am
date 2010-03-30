# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_CAMERA

if HAVE_PKG_LAPACK
lsm_headers = Extrinsics.h LinescanModel.h LinearPushbroomModel.h       \
              OrbitingPushbroomModel.h PinholeModelCalibrate.h          \
              LensDistortion.h CameraGeometry.h


ba_headers = BundleAdjust.h BundleAdjustModelBase.h                     \
             BundleAdjustReport.h BundleAdjustmentBase.h                \
             BundleAdjustmentRef.h BundleAdjustmentRobustRef.h          \
             BundleAdjustmentSparse.h BundleAdjustmentRobustSparse.h    \
             BundleAdjustmentRobustSparseKGCP.h

ba_sources = BundleAdjustReport.cc
endif

if MAKE_MODULE_INTERESTPOINT
if MAKE_MODULE_STEREO
# These deal with loading of interest points
relation_headers = CameraRelation.h ControlNetworkLoader.h

relation_sources = CameraRelation.cc ControlNetworkLoader.cc
endif
endif

include_HEADERS = CameraModel.h CAHVModel.h CAHVORModel.h               \
             CameraTransform.h Exif.h ExifData.h BayerFilter.h          \
             ControlNetwork.h PinholeModel.h $(lsm_headers)             \
             $(ba_headers) $(relation_headers)

libvwCamera_la_SOURCES = PinholeModel.cc CAHVModel.cc CAHVORModel.cc    \
             Exif.cc ExifData.cc ControlNetwork.cc LensDistortion.cc    \
             $(ba_sources) $(relation_sources)

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

