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

if MAKE_MODULE_FILEIO

if HAVE_PKG_PNG
png_headers = DiskImageResourcePNG.h
png_sources = DiskImageResourcePNG.cc
endif

if HAVE_PKG_JPEG
jpeg_headers = DiskImageResourceJPEG.h
jpeg_sources = DiskImageResourceJPEG.cc
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

include_HEADERS = DiskImageResource.h DiskImageView.h JP2.h	\
	DiskImageResourcePDS.h $(png_headers) $(jpeg_headers)	\
	$(tiff_headers) $(openexr_headers) $(hdf_headers)	\
	$(gdal_headers)

libvwFileIO_la_SOURCES = DiskImageResource.cc DiskImageResourcePDS.cc	\
	$(png_sources) $(jpeg_sources) $(tiff_sources)			\
	$(openexr_sources) $(hdf_sources) $(gdal_sources)

libvwFileIO_la_LIBADD = @MODULE_FILEIO_LIBS@

lib_LTLIBRARIES = libvwFileIO.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ -version-info @LIBTOOL_VERSION@

SUBDIRS = . tests

includedir = $(prefix)/include/vw/FileIO

include $(top_srcdir)/config/rules.mak
