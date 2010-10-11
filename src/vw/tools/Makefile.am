# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

# For simplicity, these example applications depend on exceptions
if ENABLE_EXCEPTIONS

include_HEADERS = Common.h
libvwTools_la_SOURCES = Common.cc
libvwTools_la_LIBADD = @PKG_TOOLS_LIBS@
libvwTools_la_LDFLAGS = @VW_LDFLAGS@ $(LIBTOOL_LIB_FLAGS)
lib_LTLIBRARIES = libvwTools.la

COMMON_LIBS = libvwTools.la

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
hillshade_SOURCES = hillshade.cc
hillshade_LDADD = @PKG_CARTOGRAPHY_LIBS@ $(COMMON_LIBS)
# Converts a DEM to a simple color map
colormap_SOURCES = colormap.cc
colormap_LDADD = @PKG_CARTOGRAPHY_LIBS@ $(COMMON_LIBS)
# Calculates aspect and gradient angle for DEM
slopemap_SOURCES = slopemap.cc
slopemap_LDADD = @PKG_CARTOGRAPHY_LIBS@ $(COMMON_LIBS)
# Adds or adjusts an image's georeferencing information
georef_SOURCES = georef.cc
georef_LDADD = @PKG_CARTOGRAPHY_LIBS@ $(COMMON_LIBS)
# Fades an image's alpha channel for blending
grassfirealpha_SOURCES = grassfirealpha.cc
grassfirealpha_LDADD = @PKG_CARTOGRAPHY_LIBS@ $(COMMON_LIBS)
if HAVE_PKG_GDAL
cart_gdal_progs = georef hillshade
endif
cartography_progs = colormap $(cart_gdal_progs) slopemap grassfirealpha
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
mosaic_progs = blend
# Blends multiple images into a composite
blend_SOURCES = blend.cc
blend_LDADD = @PKG_MOSAIC_LIBS@ $(COMMON_LIBS)
endif

# Command-line tools based on the Cartography and Mosaic modules
if MAKE_MODULE_CARTOGRAPHY
if MAKE_MODULE_MOSAIC
if HAVE_PKG_GDAL
cart_mos_progs = image2qtree geoblend image2toast
# Generates a quadtree of some sort from a georeferenced image
image2qtree_SOURCES = image2qtree.cc
image2qtree_LDADD = @PKG_CARTOGRAPHY_LIBS@ @PKG_MOSAIC_LIBS@ $(COMMON_LIBS)
# Creates a blended composite from multiple georeferenced images
geoblend_SOURCES = geoblend.cc
geoblend_LDADD = @PKG_CARTOGRAPHY_LIBS@ @PKG_MOSAIC_LIBS@ $(COMMON_LIBS)
# Creates a TOAST quadtree from multiple georeferenced images
image2toast_SOURCES = image2toast.cc
image2toast_LDADD = @PKG_CARTOGRAPHY_LIBS@ @PKG_MOSAIC_LIBS@ $(COMMON_LIBS)
endif
endif
endif

# Command-line tools based on the Stereo module
if MAKE_MODULE_STEREO
if MAKE_MODULE_INTERESTPOINT
stereo_progs = correlate
# Apply the block stereo correlator to two images, producing a disparity map
correlate_SOURCES = correlate.cc
correlate_LDADD = @PKG_STEREO_LIBS@ @PKG_INTERESTPOINT_LIBS@ $(COMMON_LIBS)
endif
endif

# Command-line tools based on the GPU module
if MAKE_MODULE_GPU
gpu_progs = gpu_correlate
# A GPU-based version of the Stereo correlator
gpu_correlate_SOURCES = gpu_correlate.cc
gpu_correlate_LDADD = @PKG_GPU_LIBS@ $(COMMON_LIBS)
endif

# Bundle adjustment test programs
if MAKE_MODULE_CAMERA
if MAKE_MODULE_STEREO
if MAKE_MODULE_INTERESTPOINT
if MAKE_MODULE_BUNDLEADJUSTMENT
batest_progs = ba_test make_ba_test_data
ba_test_SOURCES = ba_test.cc
ba_test_LDADD = @PKG_BUNDLEADJUSTMENT_LIBS@ $(COMMON_LIBS)
make_ba_test_data_SOURCES = make_ba_test_data.cc
make_ba_test_data_LDADD = @PKG_BUNDLEADJUSTMENT_LIBS@ $(COMMON_LIBS)
cnet_editor_progs = cnet_editor
cnet_editor_SOURCES = cnet_editor.cc
cnet_editor_LDADD = @PKG_BUNDLEADJUSTMENT_LIBS@ $(COMMON_LIBS)
endif
endif
endif
endif

# Contour generation program
if HAVE_PKG_CAIROMM
contourgen_progs = contourgen
contourgen_SOURCES = contour.cpp FitCurves.cpp contourgen.cc contour.h
contourgen_CPPFLAGS = @VW_CPPFLAGS@ @PKG_CAIROMM_CPPFLAGS@
contourgen_LDADD = $(COMMON_LIBS) @PKG_CAIROMM_LIBS@
endif

if MAKE_MODULE_MOSAIC
doc_generate_progs = doc-generate
doc_generate_SOURCES = doc-generate.cc
doc_generate_LDADD   = @PKG_VW_LIBS@ @PKG_MOSAIC_LIBS@
endif

bin_PROGRAMS = $(camera_progs) $(cartography_progs) $(hdr_progs) \
               $(interestpoint_progs) $(mosaic_progs)            \
               $(cart_mos_progs) $(stereo_progs) $(gpu_progs)    \
               $(cnet_editor_progs) $(contourgen_progs)

noinst_PROGRAMS      = $(doc_generate_progs) $(batest_progs)
dist_noinst_SCRIPTS  = ba_unit_test run_ba_tests

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@

EXTRA_DIST = run_ba_tests_lib ba_test.cfg.example ba_test_plan.txt.example

includedir = $(prefix)/include/vw/tools

include $(top_srcdir)/config/rules.mak
