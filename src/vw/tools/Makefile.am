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

# For simplicity, these example applications depend on exceptions
if ENABLE_EXCEPTIONS

include_HEADERS = Common.h hillshade.h
libvwTools_la_SOURCES = Common.cc
libvwTools_la_LIBADD = @PKG_TOOLS_LIBS@
libvwTools_la_LDFLAGS = @VW_LDFLAGS@ $(LIBTOOL_LIB_FLAGS)
lib_LTLIBRARIES = libvwTools.la

COMMON_LIBS = libvwTools.la
noinst_LIBRARIES =

# Command-line tools based on the Camera module
if MAKE_MODULE_CAMERA
camera_progs = undistort_image convert_pinhole_model

if MAKE_APP_LEGACY
camera_progs += print_exif bayer
endif

# Prints out an image's EXIF headers
print_exif_SOURCES = print_exif.cc
print_exif_LDADD = @PKG_CAMERA_LIBS@ $(COMMON_LIBS)


# Applies an inverse Bayer filter to a raw camera image
bayer_SOURCES = bayer.cc
bayer_LDADD = @PKG_CAMERA_LIBS@ $(COMMON_LIBS)

# Correct pinhole camera lens distortion
undistort_image_SOURCES = undistort_image.cc
undistort_image_LDADD = @PKG_VW_LIBS@ @PKG_CARTOGRAPHY_LIBS@ @PKG_CAMERA_LIBS@ $(COMMON_LIBS)

# Convert a pinhole model to a different type
convert_pinhole_model_SOURCES = convert_pinhole_model.cc
convert_pinhole_model_LDADD = @PKG_CAMERA_LIBS@ $(COMMON_LIBS)
endif

# Command-line tools based on the Cartography module
if MAKE_MODULE_CARTOGRAPHY
# Converts a DEM to a simple shaded-relief map
hillshade_SOURCES = hillshade.cc hillshade.h
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
# Adds or adjusts an image's georeferencing information
georefquery_SOURCES = georef_query.cc
georefquery_LDADD = @PKG_CARTOGRAPHY_LIBS@ $(COMMON_LIBS)
# Fades an image's alpha channel for blending
grassfirealpha_SOURCES = grassfirealpha.cc
grassfirealpha_LDADD = @PKG_CARTOGRAPHY_LIBS@ $(COMMON_LIBS)
if HAVE_PKG_GDAL
cart_gdal_progs = georef hillshade
endif
cartography_progs = colormap $(cart_gdal_progs) georefquery
if MAKE_APP_LEGACY
cartography_progs += slopemap grassfirealpha
endif
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
interestpoint_progs = ipalign ipmatch ipfind
if MAKE_APP_LEGACY
interestpoint_progs += learnpca
endif
#
ipfind_SOURCES = ipfind.cc
ipfind_LDADD = @PKG_INTERESTPOINT_LIBS@ $(COMMON_LIBS)
#
ipalign_SOURCES = ipalign.cc
ipalign_LDADD = @PKG_INTERESTPOINT_LIBS@ @PKG_CAMERA_LIBS@ $(COMMON_LIBS)
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
if MAKE_APP_LEGACY
mosaic_progs = blend
# Blends multiple images into a composite
blend_SOURCES = blend.cc
blend_LDADD = @PKG_MOSAIC_LIBS@ $(COMMON_LIBS)
endif
endif

# Command-line tools based on the Cartography and Mosaic modules
if MAKE_MODULE_CARTOGRAPHY
if MAKE_MODULE_MOSAIC
if HAVE_PKG_GDAL
cart_mos_progs = image2qtree
if MAKE_APP_LEGACY
cart_mos_progs += image2toast
endif

# Generates a quadtree of some sort from a georeferenced image
image2qtree_SOURCES = image2qtree.cc
image2qtree_ARCHIVES = libimage2qtree_grayau8.a libimage2qtree_grayau16.a \
        libimage2qtree_grayai16.a libimage2qtree_grayaf32.a               \
        libimage2qtree_rgbau8.a libimage2qtree_rgbau16.a                  \
        libimage2qtree_rgbai16.a libimage2qtree_rgbaf32.a
image2qtree_LDADD = $(image2qtree_ARCHIVES) @PKG_CARTOGRAPHY_LIBS@        \
        @PKG_MOSAIC_LIBS@ $(COMMON_LIBS)
# The many archive files are to help speed compile times
noinst_LIBRARIES += $(image2qtree_ARCHIVES)
libimage2qtree_grayau8_a_SOURCES = image2qtree_help.cc
libimage2qtree_grayau8_a_CPPFLAGS = $(AM_CPPFLAGS) -DPIXEL_TYPE=PixelGrayA -DCHANNEL_TYPE=uint8
libimage2qtree_grayau16_a_SOURCES = image2qtree_help.cc
libimage2qtree_grayau16_a_CPPFLAGS = $(AM_CPPFLAGS) -DPIXEL_TYPE=PixelGrayA -DCHANNEL_TYPE=uint16
libimage2qtree_grayai16_a_SOURCES = image2qtree_help.cc
libimage2qtree_grayai16_a_CPPFLAGS = $(AM_CPPFLAGS) -DPIXEL_TYPE=PixelGrayA -DCHANNEL_TYPE=int16
libimage2qtree_grayaf32_a_SOURCES = image2qtree_help.cc
libimage2qtree_grayaf32_a_CPPFLAGS = $(AM_CPPFLAGS) -DPIXEL_TYPE=PixelGrayA -DCHANNEL_TYPE=float32
libimage2qtree_rgbau8_a_SOURCES = image2qtree_help.cc
libimage2qtree_rgbau8_a_CPPFLAGS = $(AM_CPPFLAGS) -DPIXEL_TYPE=PixelRGBA -DCHANNEL_TYPE=uint8
libimage2qtree_rgbau16_a_SOURCES = image2qtree_help.cc
libimage2qtree_rgbau16_a_CPPFLAGS = $(AM_CPPFLAGS) -DPIXEL_TYPE=PixelRGBA -DCHANNEL_TYPE=uint16
libimage2qtree_rgbai16_a_SOURCES = image2qtree_help.cc
libimage2qtree_rgbai16_a_CPPFLAGS = $(AM_CPPFLAGS) -DPIXEL_TYPE=PixelRGBA -DCHANNEL_TYPE=int16
libimage2qtree_rgbaf32_a_SOURCES = image2qtree_help.cc
libimage2qtree_rgbaf32_a_CPPFLAGS = $(AM_CPPFLAGS) -DPIXEL_TYPE=PixelRGBA -DCHANNEL_TYPE=float32

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
correlate_LDADD = @PKG_VW_LIBS@ @PKG_STEREO_LIBS@ @PKG_INTERESTPOINT_LIBS@ @PKG_CARTOGRAPHY_LIBS@ $(COMMON_LIBS)
endif
endif

# Contour generation program
if HAVE_PKG_CAIROMM
if MAKE_APP_LEGACY
contourgen_progs = contourgen
contourgen_SOURCES = contour.cpp FitCurves.cpp contourgen.cc contour.h
contourgen_CPPFLAGS = @VW_CPPFLAGS@ @PKG_CAIROMM_CPPFLAGS@
contourgen_LDADD = $(COMMON_LIBS) @PKG_CAIROMM_LIBS@
endif
endif

if MAKE_MODULE_MOSAIC
if MAKE_APP_LEGACY
doc_generate_progs = doc-generate
doc_generate_SOURCES = doc-generate.cc
doc_generate_LDADD   = @PKG_VW_LIBS@ @PKG_MOSAIC_LIBS@
endif
endif

if MAKE_MODULE_FLOOD_DETECT
detect_water_progs = detect_water clean_sentinel1_borders
detect_water_SOURCES = detect_water.cc multispectral.h radar.h flood_common.h
detect_water_LDADD   = @PKG_VW_LIBS@ @PKG_CARTOGRAPHY_LIBS@ $(COMMON_LIBS)
clean_sentinel1_borders_SOURCES = clean_sentinel1_borders.cc flood_common.h
clean_sentinel1_borders_LDADD   = @PKG_VW_LIBS@ @PKG_CARTOGRAPHY_LIBS@ $(COMMON_LIBS)
endif

bin_PROGRAMS = $(camera_progs) $(cartography_progs) $(hdr_progs) \
               $(interestpoint_progs) $(mosaic_progs)            \
               $(cart_mos_progs) $(stereo_progs)     \
               $(contourgen_progs) $(detect_water_progs)

noinst_PROGRAMS      = $(doc_generate_progs)

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@

includedir = $(prefix)/include/vw/tools

include $(top_srcdir)/config/rules.mak
