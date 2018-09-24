# __BEGIN_LICENSE__
#  Copyright (c) 2006-2013, United States Government as represented by the
#  Administrator of the National Aeronautics and Space Administration. All
#  rights reserved.
#
#  The NASA Vision Workbench is licensed under the Apache License,
#  Version 2.0 (the "License"); you may not use this file except in
#  compliance with the License. You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_CAMERA

if HAVE_PKG_LAPACK
lapack_headers = \
	CameraGeometry.h \
	CameraSolve.h \
	Extrinsics.h \
	LensDistortion.h \
	LinescanModel.h

lapack_sources = \
	CameraGeometry.cc \
	LinescanModel.cc \
	Extrinsics.cc
endif

include_HEADERS =   \
  BayerFilter.h     \
  CAHVModel.h       \
  CAHVOREModel.h    \
  CAHVORModel.h     \
  CameraModel.h     \
  CameraTransform.h \
  CameraUtilities.h \
  ExifData.h        \
  Exif.h            \
  PinholeModel.h    \
  $(lapack_headers)

libvwCamera_la_SOURCES = \
  CAHVModel.cc       \
  CAHVOREModel.cc    \
  CAHVORModel.cc     \
  CameraModel.cc     \
  CameraUtilities.cc \
  Exif.cc            \
  ExifData.cc        \
  LensDistortion.cc  \
  PinholeModel.cc    \
  $(lapack_sources)

libvwCamera_la_LIBADD = @MODULE_CAMERA_LIBS@

lib_LTLIBRARIES = libvwCamera.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ $(LIBTOOL_LIB_FLAGS)

SUBDIRS = . tests

includedir = $(prefix)/include/vw/Camera

include $(top_srcdir)/config/rules.mak

