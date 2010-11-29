# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_FILEIO

if HAVE_PKG_PNG
png_headers = DiskImageResourcePNG.h  MemoryImageResourcePNG.h  PngIO.h
png_sources = DiskImageResourcePNG.cc MemoryImageResourcePNG.cc PngIO.cc
endif

if HAVE_PKG_JPEG
jpeg_headers = DiskImageResourceJPEG.h  MemoryImageResourceJPEG.h  JpegIO.h
jpeg_sources = DiskImageResourceJPEG.cc MemoryImageResourceJPEG.cc JpegIO.cc
endif

if HAVE_PKG_TIFF
tiff_headers = DiskImageResourceTIFF.h
tiff_sources = DiskImageResourceTIFF.cc
endif

if HAVE_PKG_OPENEXR
if ENABLE_EXCEPTIONS
openexr_headers = DiskImageResourceOpenEXR.h
openexr_sources = DiskImageResourceOpenEXR.cc
endif
endif

if HAVE_PKG_HDF
hdf_headers = DiskImageResourceHDF.h
hdf_sources = DiskImageResourceHDF.cc
endif

if HAVE_PKG_GDAL
gdal_headers = DiskImageResourceGDAL.h
gdal_sources = DiskImageResourceGDAL.cc
endif

include_HEADERS = \
  DiskImageResource.h \
  DiskImageResourcePBM.h \
  DiskImageResourcePDS.h \
  DiskImageView.h \
  MemoryImageResource.h \
  KML.h \
  ScanlineIO.h \
  $(gdal_headers) \
  $(hdf_headers) \
  $(jpeg_headers) \
  $(openexr_headers) \
  $(png_headers) \
  $(tiff_headers)

libvwFileIO_la_SOURCES = \
  DiskImageResource.cc \
  DiskImageResourcePBM.cc \
  DiskImageResourcePDS.cc \
  KML.cc \
  MemoryImageResource.cc \
  ScanlineIO.cc \
  $(gdal_sources) \
  $(hdf_sources) \
  $(jpeg_sources) \
  $(openexr_sources) \
  $(png_sources) \
  $(tiff_sources)

libvwFileIO_la_LIBADD = @MODULE_FILEIO_LIBS@

lib_LTLIBRARIES = libvwFileIO.la

noinst_HEADERS = DiskImageResource_internal.h

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ $(LIBTOOL_LIB_FLAGS)

SUBDIRS = . tests

includedir = $(prefix)/include/vw/FileIO

include $(top_srcdir)/config/rules.mak
