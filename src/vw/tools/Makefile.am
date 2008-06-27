# __BEGIN_LICENSE__
# 
# Copyright (C) 2006 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration
# (NASA).  All Rights Reserved.
# 
# Copyright 2006 Carnegie Mellon University. All rights reserved.
# 
# This software is distributed under the NASA Open Source Agreement
# (NOSA), version 1.3.  The NOSA has been approved by the Open Source
# Initiative.  See the file COPYING at the top of the distribution
# directory tree for the complete NOSA document.
# 
# THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
# KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
# LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
# SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
# A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
# THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
# DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
# 
# __END_LICENSE__

########################################################################
# sources
########################################################################

# For simplicity, these example applications depend on exceptions
if ENABLE_EXCEPTIONS

COMMON_LIBS = @PKG_VW_LIBS@ @PKG_BOOST_FILESYSTEM_LIBS@ \
	@PKG_BOOST_PROGRAM_OPTIONS_LIBS@ @PKG_BOOST_THREAD_LIBS@

# Command-line tools based on the Camera module
if MAKE_MODULE_CAMERA
camera_progs = print_exif bayer
# Prints out an image's EXIF headers
print_exif_SOURCES = print_exif.cc 
print_exif_LDADD = @PKG_CAMERA_LIBS@ $(COMMON_LIBS)
# Applies an inverse Bayer filter to a raw camera image
bayer_SOURCES = bayer.cc
bayer_LDADD = @PKG_CAMERA_LIBS@ $(COMMON_LIBS)
endif

# Command-line tools based on the Cartography module
if MAKE_MODULE_CARTOGRAPHY
# Converts a DEM to a simple shaded-relief map
shaded_relief_SOURCES = shaded_relief.cc
shaded_relief_LDADD = @PKG_CARTOGRAPHY_LIBS@ $(COMMON_LIBS)
# Converts a DEM to a simple color map
colorize_dem_SOURCES = colorize_dem.cc
colorize_dem_LDADD = @PKG_CARTOGRAPHY_LIBS@ $(COMMON_LIBS)
# Adds or adjusts an image's georeferencing information
georef_SOURCES = georef.cc
georef_LDADD = @PKG_CARTOGRAPHY_LIBS@ $(COMMON_LIBS)
if HAVE_PKG_GDAL
cart_gdal_progs = georef shaded_relief
endif
cartography_progs = colorize_dem $(cart_gdal_progs)
endif

# Command-line tools based on the HDR module
if MAKE_MODULE_HDR
hdr_progs = hdr_tonemap hdr_merge
# Merge multiple source images into an HDR image
hdr_merge_SOURCES = hdr_merge.cc 
hdr_merge_LDADD = @PKG_HDR_LIBS@ $(COMMON_LIBS)
# Map an HDR image into an 8-bit dynamic range
hdr_tonemap_SOURCES = hdr_tonemap.cc 
hdr_tonemap_LDADD = @PKG_HDR_LIBS@ $(COMMON_LIBS)
endif

# Command-line tools based on the InterestPoint module
if MAKE_MODULE_INTERESTPOINT 
if HAVE_PKG_LAPACK
interestpoint_progs = ipalign ipmatch ipfind learnpca
# 
ipfind_SOURCES = ipfind.cc
ipfind_LDADD = @PKG_INTERESTPOINT_LIBS@ $(COMMON_LIBS)
# 
ipalign_SOURCES = ipalign.cc
ipalign_LDADD = @PKG_INTERESTPOINT_LIBS@ $(COMMON_LIBS)
#
ipmatch_SOURCES = ipmatch.cc
ipmatch_LDADD = @PKG_INTERESTPOINT_LIBS@ $(COMMON_LIBS)
#
learnpca_SOURCES = learnpca.cc
learnpca_LDADD = @PKG_INTERESTPOINT_LIBS@ $(COMMON_LIBS)
#
endif
endif

# Command-line tools based on the Mosaic module
if MAKE_MODULE_MOSAIC
mosaic_progs = qtree blend
# Generates a quad-tree representation of an image
qtree_SOURCES = qtree.cc
qtree_LDADD = @PKG_MOSAIC_LIBS@ $(COMMON_LIBS)
# Blends multiple images into a composite
blend_SOURCES = blend.cc
blend_LDADD = @PKG_MOSAIC_LIBS@ $(COMMON_LIBS)
endif

# Command-line tools based on the Cartography and Mosaic modules
if MAKE_MODULE_CARTOGRAPHY
if MAKE_MODULE_MOSAIC
if HAVE_PKG_GDAL
cart_mos_progs = image2kml image2tms geoblend
# Generates a KML superoverlay from a georeferenced image
image2kml_SOURCES = image2kml.cc
image2kml_LDADD = @PKG_CARTOGRAPHY_LIBS@ @PKG_MOSAIC_LIBS@ $(COMMON_LIBS)
# Generates a TMS map heirarchy from a georeferenced image
image2tms_SOURCES = image2tms.cc
image2tms_LDADD = @PKG_CARTOGRAPHY_LIBS@ @PKG_MOSAIC_LIBS@ $(COMMON_LIBS)
# Creates a blended composite from multiple georeferenced images
geoblend_SOURCES = geoblend.cc
geoblend_LDADD = @PKG_CARTOGRAPHY_LIBS@ @PKG_MOSAIC_LIBS@ $(COMMON_LIBS)
endif
endif
endif

# Command-line tools based on the Stereo module
if MAKE_MODULE_STEREO
stereo_progs = correlate
# Apply the block stereo correlator to two images, producing a disparity map
correlate_SOURCES = correlate.cc
correlate_LDADD = @PKG_STEREO_LIBS@ $(COMMON_LIBS)
endif

# Command-line tools based on the GPU module
if MAKE_MODULE_GPU
gpu_progs = gpu_correlate
# A GPU-based version of the Stereo correlator
gpu_correlate_SOURCES = gpu_correlate.cc
gpu_correlate_LDADD = @PKG_GPU_LIBS@ $(COMMON_LIBS)
endif

bin_PROGRAMS = $(camera_progs) $(cartography_progs) $(hdr_progs)	\
               $(interestpoint_progs) $(mosaic_progs)			\
               $(cart_mos_progs) $(stereo_progs) $(gpu_progs)	

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@

includedir = $(prefix)/include/vw

include $(top_srcdir)/config/rules.mak
