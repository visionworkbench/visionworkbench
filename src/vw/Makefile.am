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

if MAKE_MODULE_CORE
core_subdirs = Core
endif

if MAKE_MODULE_MATH
math_subdirs = Math
endif

if MAKE_MODULE_IMAGE
image_subdirs = Image
endif

if MAKE_MODULE_FILEIO
fileio_subdirs = FileIO
endif

if MAKE_MODULE_CAMERA
camera_subdirs = Camera
endif

if MAKE_MODULE_MOSAIC
mosaic_subdirs = Mosaic
endif

if MAKE_MODULE_STEREO
stereo_subdirs = Stereo
endif

if MAKE_MODULE_CARTOGRAPHY
cartography_subdirs = Cartography
endif

if MAKE_MODULE_INTERESTPOINT
interestpoint_subdirs = InterestPoint
endif

if MAKE_MODULE_HDR
hdr_subdirs = HDR
endif

if MAKE_MODULE_GEOMETRY
geometry_subdirs = Geometry
endif


if MAKE_MODULE_BUNDLEADJUSTMENT
bundleadjustment_subdirs = BundleAdjustment
endif

if MAKE_MODULE_TOOLS
tools_subdirs = tools
endif

if MAKE_MODULE_PYTHON
python_subdirs = Python
endif

include_HEADERS = config.h


SUBDIRS = \
  $(core_subdirs) \
  $(math_subdirs) \
  $(image_subdirs) \
  $(fileio_subdirs) \
  . \
  $(camera_subdirs) \
  $(stereo_subdirs) \
  $(cartography_subdirs) \
  $(mosaic_subdirs) \
  $(interestpoint_subdirs) \
  $(hdr_subdirs) \
  $(geometry_subdirs) \
  $(bundleadjustment_subdirs) \
  $(python_subdirs) \
  $(tools_subdirs)

lib_LTLIBRARIES = libvw.la
libvw_la_SOURCES =
libvw_la_LDFLAGS= $(LIBTOOL_LIB_FLAGS)
libvw_la_LIBADD = Core/libvwCore.la Math/libvwMath.la Image/libvwImage.la FileIO/libvwFileIO.la
# This makes sure we're using the c++ linker
nodist_EXTRA_libvw_la_SOURCES = dummy.cxx

CLEANFILES = config.h.pre.in~

########################################################################
# general
########################################################################

includedir = $(prefix)/include/vw

include $(top_srcdir)/config/rules.mak
DISTCLEANFILES += config.h
