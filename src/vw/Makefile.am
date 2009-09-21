# __BEGIN_LICENSE__
# Copyright (C) 2006, 2007 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_CORE
core_subdirs = Core
core_headers = Core.h
endif

if MAKE_MODULE_MATH
math_subdirs = Math
math_headers = Math.h
endif

if MAKE_MODULE_IMAGE
image_subdirs = Image
image_headers = Image.h
endif

if MAKE_MODULE_FILEIO
fileio_subdirs = FileIO
fileio_headers = FileIO.h
endif

if MAKE_MODULE_CAMERA
camera_subdirs = Camera
camera_headers = Camera.h
endif

if MAKE_MODULE_MOSAIC
mosaic_subdirs = Mosaic
mosaic_headers = Mosaic.h
endif

if MAKE_MODULE_STEREO
stereo_subdirs = Stereo
stereo_headers = Stereo.h
endif

if MAKE_MODULE_CARTOGRAPHY
cartography_subdirs = Cartography
cartography_headers = Cartography.h
endif

if MAKE_MODULE_GPU
gpu_subdirs = GPU
gpu_headers = GPU.h
endif

if MAKE_MODULE_INTERESTPOINT
interestpoint_subdirs = InterestPoint
interestpoint_headers = InterestPoint.h
endif

if MAKE_MODULE_HDR
hdr_subdirs = HDR
hdr_headers = HDR.h
endif

if MAKE_MODULE_GEOMETRY
geometry_subdirs = Geometry
geometry_headers = Geometry.h
endif

if MAKE_MODULE_TOOLS
tools_subdirs = tools
endif

if MAKE_MODULE_PLATE
plate_subdirs = plate
endif

if MAKE_MODULE_GUI
gui_subdirs = gui
endif

if MAKE_MODULE_PYTHON
python_subdirs = Python
endif

include_HEADERS = $(core_headers) $(math_headers) $(image_headers)	\
	$(fileio_headers) $(camera_headers) $(mosaic_headers)		\
	$(stereo_headers) $(cartography_headers) $(gpu_headers)		\
	$(interestpoint_headers) $(hdr_headers) $(geometry_headers)	\
	$(plate_headers) vw.h config.h

SUBDIRS = $(core_subdirs) $(math_subdirs) $(image_subdirs)		\
	$(fileio_subdirs) . $(camera_subdirs) $(mosaic_subdirs)		\
	$(stereo_subdirs) $(cartography_subdirs) $(gpu_subdirs)		\
	$(interestpoint_subdirs) $(hdr_subdirs) $(geometry_subdirs)	\
	$(plate_subdirs) $(tools_subdirs) $(gui_subdirs)		\
	$(python_subdirs)

lib_LTLIBRARIES = libvw.la
libvw_la_SOURCES =
libvw_la_LDFLAGS= -version-info @LIBTOOL_VERSION@
libvw_la_LIBADD = Core/libvwCore.la Math/libvwMath.la Image/libvwImage.la FileIO/libvwFileIO.la

CLEANFILES = config.h.pre.in~

########################################################################
# general
########################################################################

includedir = $(prefix)/include/vw

include $(top_srcdir)/config/rules.mak
